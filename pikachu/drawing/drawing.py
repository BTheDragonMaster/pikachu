#!/usr/bin/env python
from typing import Union, List, Tuple, Dict, Set, Generator, Optional, TYPE_CHECKING
import math
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from io import StringIO
import re

from pikachu.drawing.rings import Ring, RingOverlap, find_neighbouring_rings, rings_connected_by_bridge, \
    find_bridged_systems
from pikachu.math_functions import Vector, Polygon, Line
from pikachu.chem.chirality import get_chiral_permutations
from pikachu.chem.atom_properties import ATOM_PROPERTIES
from pikachu.chem.atom import Atom
from pikachu.chem.bond import Bond
from pikachu.chem.lone_pair import LonePair
from pikachu.chem.structure import Structure
from pikachu.errors import DrawingError
from pikachu.drawing.colours import RANDOM_PALETTE_2, get_hex

if TYPE_CHECKING:
    pass


class Options:
    def __init__(self):
        self.width: int = 500
        self.height: int = 500
        self.bond_thickness: int = 2
        self.bond_length: int = 15
        self.chiral_bond_width: float = self.bond_length * 0.1
        self.bond_length_squared: int = self.bond_length ** 2
        self.short_bond_length: float = 0.50
        self.double_bond_length: float = 0.80
        self.bond_spacing: float = 0.18 * self.bond_length
        self.isomeric: bool = True
        self.padding: int = 30
        self.font_size_large: int = 5
        self.font_size_small: int = 3
        self.kk_threshold: float = 0.1
        self.kk_inner_threshold: float = 0.1
        self.kk_max_iteration: int = 2000
        self.kk_max_inner_iteration: int = 50
        self.kk_max_energy: float = 1e9
        self.overlap_sensitivity: float = 0.10
        self.overlap_resolution_iterations: int = 5
        self.background_color: str = 'white'
        self.draw_hydrogens: bool = False
        self.finetune: bool = True
        self.strict_mode: bool = False
        self.svg_font: str = "verdana"
        self.svg_font_size: int = 8
        self.svg_font_size_small: int = 6
        self.svg_letter_spacing: int = -2


class KKLayout:

    def __init__(self, structure, atoms, center, start_atom, bond_length,
                 threshold=0.1, inner_threshold=0.1, max_iteration=2000,
                 max_inner_iteration=50, max_energy=1e9):
        self.structure = structure
        self.atoms = atoms
        self.center = center
        self.start_atom = start_atom
        self.edge_strength = bond_length
        self.threshold = threshold
        self.inner_threshold = inner_threshold
        self.max_iteration = max_iteration
        self.max_inner_iteration = max_inner_iteration
        self.max_energy = max_energy

        self.x_positions = {}
        self.y_positions = {}
        self.positioned = {}

        self.length_matrix = {}
        self.distance_matrix = {}
        self.spring_strengths = {}
        self.energy_matrix = {}

        self.energy_sums_x = {}
        self.energy_sums_y = {}

        self.initialise_matrices()
        self.get_kk_layout()

    def initialise_matrices(self):

        self.distance_matrix = self.get_subgraph_distance_matrix(self.atoms)
        length = len(self.atoms)

        radius = Polygon.find_polygon_radius(500, length)
        angle = Polygon.get_central_angle(length)

        a: float = 0.0

        for atom in self.atoms:
            if not atom.draw.positioned:
                self.x_positions[atom] = self.center.x + math.cos(a) * radius
                self.y_positions[atom] = self.center.y + math.sin(a) * radius
            else:
                self.x_positions[atom] = atom.draw.position.x
                self.y_positions[atom] = atom.draw.position.y

            self.positioned[atom] = atom.draw.positioned
            a += angle

        for atom_1 in self.atoms:
            self.length_matrix[atom_1] = {}
            self.spring_strengths[atom_1] = {}
            self.energy_matrix[atom_1] = {}

            self.energy_sums_x[atom_1] = None
            self.energy_sums_y[atom_1] = None

            for atom_2 in self.atoms:
                self.length_matrix[atom_1][atom_2] = self.edge_strength * self.distance_matrix[atom_1][atom_2]
                self.spring_strengths[atom_1][atom_2] = self.edge_strength * self.distance_matrix[atom_1][atom_2] ** -2.0
                self.energy_matrix[atom_1][atom_2] = None

        for atom_1 in self.atoms:
            ux = self.x_positions[atom_1]
            uy = self.y_positions[atom_1]
            d_ex = 0.0
            d_ey = 0.0

            for atom_2 in self.atoms:

                if atom_1 == atom_2:
                    continue

                vx = self.x_positions[atom_2]
                vy = self.y_positions[atom_2]

                denom = 1.0 / math.sqrt((ux - vx) ** 2 + (uy - vy) ** 2)

                self.energy_matrix[atom_1][atom_2] = (self.spring_strengths[atom_1][atom_2] * ((ux - vx) - self.length_matrix[atom_1][atom_2] * (ux - vx) * denom),
                                                      self.spring_strengths[atom_1][atom_2] * ((uy - vy) - self.length_matrix[atom_1][atom_2] * (uy - vy) * denom))

                self.energy_matrix[atom_2][atom_1] = self.energy_matrix[atom_1][atom_2]

                d_ex += self.energy_matrix[atom_1][atom_2][0]
                d_ey += self.energy_matrix[atom_1][atom_2][1]

            self.energy_sums_x[atom_1] = d_ex
            self.energy_sums_y[atom_1] = d_ey

    def get_kk_layout(self):
        iteration = 0

        max_energy = self.max_energy

        while max_energy > self.threshold and self.max_iteration > iteration:
            iteration += 1
            max_energy_atom, max_energy, d_ex, d_ey = self.highest_energy()
            delta = max_energy
            inner_iteration = 0

            while delta > self.inner_threshold and self.max_inner_iteration > inner_iteration:
                inner_iteration += 1
                self.update(max_energy_atom, d_ex, d_ey)
                delta, d_ex, d_ey = self.energy(max_energy_atom)

        for atom in self.atoms:
            atom.draw.position.x = self.x_positions[atom]
            atom.draw.position.y = self.y_positions[atom]
            atom.draw.positioned = True
            atom.draw.force_positioned = True

    def energy(self, atom):
        energy = [self.energy_sums_x[atom]**2 + self.energy_sums_y[atom]**2, self.energy_sums_x[atom], self.energy_sums_y[atom]]
        return energy

    def highest_energy(self):
        max_energy = 0.0
        max_energy_atom = None
        max_d_ex = 0.0
        max_d_ey = 0.0

        for atom in self.atoms:
            delta, d_ex, d_ey = self.energy(atom)

            if delta > max_energy and not self.positioned[atom]:
                max_energy = delta
                max_energy_atom = atom
                max_d_ex = d_ex
                max_d_ey = d_ey

        return max_energy_atom, max_energy, max_d_ex, max_d_ey

    def update(self, atom, d_ex, d_ey):
        dxx = 0.0
        dyy = 0.0
        dxy = 0.0

        ux = self.x_positions[atom]
        uy = self.y_positions[atom]

        lengths_array = self.length_matrix[atom]
        strengths_array = self.spring_strengths[atom]

        for atom_2 in self.atoms:
            if atom == atom_2:
                continue

            vx = self.x_positions[atom_2]
            vy = self.y_positions[atom_2]

            length = lengths_array[atom_2]
            strength = strengths_array[atom_2]

            squared_xdiff = (ux - vx) ** 2
            squared_ydiff = (uy - vy) ** 2

            denom = 1.0 / (squared_xdiff + squared_ydiff) ** 1.5

            dxx += strength * (1 - length * squared_ydiff * denom)
            dyy += strength * (1 - length * squared_xdiff * denom)
            dxy += strength * (length * (ux - vx) * (uy - vy) * denom)

        if dxx == 0:
            dxx = 0.1

        if dyy == 0:
            dyy = 0.1

        if dxy == 0:
            dxy = 0.1

        dy = (d_ex / dxx + d_ey / dxy) / (dxy / dxx - dyy / dxy)
        dx = -(dxy * dy + d_ex) / dxx

        self.x_positions[atom] += dx
        self.y_positions[atom] += dy

        d_ex = 0.0
        d_ey = 0.0

        ux = self.x_positions[atom]
        uy = self.y_positions[atom]

        for atom_2 in self.atoms:
            if atom == atom_2:
                continue

            vx = self.x_positions[atom_2]
            vy = self.y_positions[atom_2]

            previous_ex = self.energy_matrix[atom][atom_2][0]
            previous_ey = self.energy_matrix[atom][atom_2][1]

            denom = 1.0 / math.sqrt((ux - vx) ** 2 + (uy - vy) ** 2)
            dx = strengths_array[atom_2] * ((ux - vx) - lengths_array[atom_2] * (ux - vx) * denom)
            dy = strengths_array[atom_2] * ((uy - vy) - lengths_array[atom_2] * (uy - vy) * denom)

            self.energy_matrix[atom][atom_2] = [dx, dy]

            d_ex += dx
            d_ey += dy

            self.energy_sums_x[atom_2] += dx - previous_ex
            self.energy_sums_y[atom_2] += dy - previous_ey

        self.energy_sums_x[atom] = d_ex
        self.energy_sums_y[atom] = d_ey

    def get_subgraph_distance_matrix(self, atoms):
        adjacency_matrix = self.get_subgraph_adjacency_matrix(atoms)
        distance_matrix = {}

        for atom_1 in atoms:
            if atom_1 not in distance_matrix:
                distance_matrix[atom_1] = {}
            for atom_2 in atoms:
                distance_matrix[atom_1][atom_2] = float('inf')
                if adjacency_matrix[atom_1][atom_2] == 1:
                    distance_matrix[atom_1][atom_2] = 1

        for atom_1 in atoms:
            for atom_2 in atoms:
                for atom_3 in atoms:
                    if distance_matrix[atom_2][atom_3] > distance_matrix[atom_2][atom_1] + distance_matrix[atom_1][atom_3]:
                        distance_matrix[atom_2][atom_3] = distance_matrix[atom_2][atom_1] + distance_matrix[atom_1][atom_3]

        return distance_matrix

    def get_subgraph_adjacency_matrix(self, atoms):
        adjacency_matrix = {}

        for atom_1 in atoms:
            adjacency_matrix[atom_1] = {}
            for atom_2 in atoms:
                adjacency_matrix[atom_1][atom_2] = 0

        for atom_1 in atoms:
            for atom_2 in atoms:
                if atom_1 != atom_2:
                    if atom_1 in self.structure.bond_lookup[atom_2]:
                        adjacency_matrix[atom_1][atom_2] = 1

        return adjacency_matrix


class Drawer:
    """
    Class to draw a PIKAChU structure
    """
    def __init__(self, structure: Structure, options: Union[Options, None] = None,
                 coords_only: bool = False, multiple: bool = False, kekulise: bool = True) -> None:
        """
        Draw a PIKAChU structure from a structure instance

        Parameters
        ----------
        structure: Structure instance
        options: Options instance, use this to define custom drawing parameters
        coords_only: bool, if True, only assign positions to atoms; draw the structure otherwise
        multiple: bool, if True, multiple structures are drawn, False otherwise.
             Do not toggle this parameter. If you want to draw multiple structures, use the function
             draw_multiple instead
        kekulise: bool, if True, kekulise the structure before drawing, False otherwise
        """

        if options is None:
            self.options: Options = Options()
        else:
            self.options: Options = options
        if kekulise:
            self.structure = structure.kekulise()
        else:
            self.structure = structure

        self.multiple: bool = multiple

        # Used for tracking rings in the structure
        self.rings: List[Ring] = []
        self.ring_overlaps: List[RingOverlap] = []
        self.id_to_ring: Dict[int, Ring] = {}
        self.has_bridged_ring: bool = False
        self.ring_id_tracker: int = 0
        self.ring_overlap_id_tracker: int = 0

        # Used for tracking rings in the structure once bridged ring systems have been established
        self.original_rings: List[Ring] = []
        self.original_ring_overlaps: List[RingOverlap] = []

        self.drawn_atoms: List[Atom] = []
        self.drawn_bonds: List[Bond] = []

        self.total_overlap_score: float = 0.0

        # TODO: Check if we cannot use the atom dictionary in Structure instead
        self.atom_nr_to_atom: Dict[int, Atom] = {}

        # Used to keep track of bonds that are under steric restraint
        self.chiral_bonds: List[Bond] = []
        self.chiral_bond_to_orientation: Dict[Bond, Tuple[str, Atom]] = {}
        self.fixed_chiral_bonds: Set[Bond] = set()

        # Used for grouping atoms and bonds within an SVG file
        self.svg_groups: Dict[str, Dict[str, List[str]]] = {}
        self.annotation_to_colour: Optional[Dict[str, str]] = None
        self.annotation: Optional[str] = None
        self.structure_id: Optional[str] = None

        # self.svg_style: str = """<style> line {stroke: black; stroke_width: 1px;} polygon {fill: black;} </style>"""
        self.svg_style: str = ""

        self.draw(coords_only=coords_only)

    def set_annotation_for_grouping(self, annotation: str) -> None:
        """
        Sets structure-wide annotation to name svg groups. Especially useful when drawing multiple structures

        Parameters
        ----------
        annotation: str, label for the current structure

        """
        self.svg_groups = {}
        self.annotation = annotation

    def colour_by_annotation(self, value_to_colour: Optional[Dict[str, str]] = None, verbose: bool = False) -> None:
        """
        Colour atoms by their annotations

        Parameters
        ----------
        value_to_colour: dict of {annotation: colour, ->}, with annotation and colour both str, and colour a hex code.
            If dictionary is not given, PIKAChU will use a default palette containing 18 colours
        verbose: bool, if True, print warnings when there are insufficient colours in the default palette

        """
        assert self.annotation
        self.annotation_to_colour = {}

        annotation_values: Set[str] = set()
        for atom in self.structure.graph:
            if atom.annotations.has_annotation(self.annotation):
                annotation_values.add(atom.annotations.get_annotation(self.annotation))

        annotation_values: List[str] = list(annotation_values)
        annotation_values.sort()

        if value_to_colour is not None:
            for annotation_value in annotation_values:
                if annotation_value not in value_to_colour:
                    raise ValueError(f"Could not find colour for {annotation_value}. Please provide a colour for each \
                    annotation value in the structure.")

            for annotation_value, colour in value_to_colour.items():
                self.annotation_to_colour[annotation_value] = get_hex(colour)

        else:
            if len(RANDOM_PALETTE_2) < len(annotation_values):
                if verbose:
                    print("Warning: Not enough colours in the default palette to assign a unique colour to each \
                    annotation. Use the parameter 'value_to_colour' to provide your own hex codes instead.")
            self.annotation_to_colour = {}
            for i, annotation_value in enumerate(annotation_values):
                colour = RANDOM_PALETTE_2[i % len(RANDOM_PALETTE_2)]
                self.annotation_to_colour[annotation_value] = get_hex(colour)

        for atom in self.structure.graph:
            if atom.annotations.has_annotation(self.annotation):
                atom.draw.colour = self.annotation_to_colour[atom.annotations.get_annotation(self.annotation)]

    # TODO: Should this function move to Structure?

    def find_shortest_path(self, atom_1: Atom, atom_2: Atom, path_type: str = 'bond') -> List[Union[Bond, Atom]]:
        """
        Return the shortest path between two atoms, either as a list of bonds or as a list of atoms

        Parameters
        ----------
        atom_1: Atom instance, must be in structure
        atom_2: Atom instance, must be in structure
        path_type: str, 'bond' or 'atom'

        Returns
        -------
        list of atoms or bonds describing the shortest path between atom_1 and atom_2

        """
        distances: Dict[Atom, float] = {}
        previous_hop: Dict[Atom, Optional[Atom]] = {}
        unvisited: Set[Atom] = set()

        for atom in self.structure.graph:
            distances[atom] = float('inf')
            previous_hop[atom] = None
            unvisited.add(atom)

        distances[atom_1] = 0.0

        while unvisited:

            current_atom: Optional[Atom] = None
            minimum = float('inf')

            # Find the atom with the smallest distance value that has not yet been visited
            for atom in unvisited:
                dist = distances[atom]
                if dist < minimum:
                    current_atom = atom
                    minimum = dist

            unvisited.remove(current_atom)

            if current_atom == atom_2:
                break

            # If there exists a shorter path between the source atom and the neighbours, update distance
            for neighbour in self.structure.graph[current_atom]:
                if neighbour in unvisited:
                    alternative_distance: float = distances[current_atom] + 1.0

                    if alternative_distance < distances[neighbour]:
                        distances[neighbour] = alternative_distance
                        previous_hop[neighbour] = current_atom

        # Construct the path of atoms

        path_atoms: List[Atom] = []
        current_atom: Optional[Atom] = atom_2
        if previous_hop[current_atom] or current_atom == atom_1:
            while current_atom:
                path_atoms.insert(0, current_atom)
                current_atom = previous_hop[current_atom]

        if path_type == 'bond':
            path: List[Union[Bond, Atom]] = []
            for i in range(1, len(path_atoms)):
                atom_1 = path_atoms[i - 1]
                atom_2 = path_atoms[i]
                bond = self.structure.bond_lookup[atom_1][atom_2]
                path.append(bond)

            return path
        elif path_type == 'atom':
            return path_atoms
        else:
            raise ValueError("Path type must be 'bond' or 'atom'.")

    def _finetune_overlap_resolution(self) -> None:
        """
        Find clashing atoms, find the shortest path between those atoms and rotate a bond on that path to
            resolve the overlap
        """

        if self.total_overlap_score > self.options.overlap_sensitivity:
            clashing_atoms = self._find_clashing_atoms()

            best_bonds: List[Bond] = []
            for atom_1, atom_2 in clashing_atoms:
                shortest_path = self.find_shortest_path(atom_1, atom_2)
                rotatable_bonds: List[Bond] = []
                distances: List[float] = []

                for i, bond in enumerate(shortest_path):
                    distance_1: int = i
                    distance_2: int = len(shortest_path) - i

                    average_distance = len(shortest_path) / 2

                    distance_metric = abs(average_distance - distance_1) + abs(average_distance - distance_2)

                    if self.bond_is_rotatable(bond):
                        rotatable_bonds.append(bond)
                        distances.append(distance_metric)

                best_bond: Optional[Bond] = None
                optimal_distance: float = float('inf')
                for i, distance in enumerate(distances):
                    if distance < optimal_distance:
                        best_bond = rotatable_bonds[i]
                        optimal_distance = distance

                if best_bond is not None:
                    best_bonds.append(best_bond)

            best_bonds = list(set(best_bonds))

            for best_bond in best_bonds:
                if self.total_overlap_score > self.options.overlap_sensitivity:

                    subtree_size_1 = self.get_subgraph_size(best_bond.atom_1, {best_bond.atom_2})
                    subtree_size_2 = self.get_subgraph_size(best_bond.atom_2, {best_bond.atom_1})

                    if subtree_size_1 < subtree_size_2:
                        rotating_atom = best_bond.atom_1
                        parent_atom = best_bond.atom_2
                    else:
                        rotating_atom = best_bond.atom_2
                        parent_atom = best_bond.atom_1

                    overlap_score, _, _ = self.get_overlap_score()

                    scores: List[float] = [overlap_score]

                    # Attempt 12 rotations
                    for i in range(12):
                        self.rotate_subtree(rotating_atom, parent_atom, math.radians(30), parent_atom.draw.position)
                        new_overlap_score, _, _ = self.get_overlap_score()
                        scores.append(new_overlap_score)

                    assert len(scores) == 13

                    scores = scores[:12]

                    best_i = 0
                    best_score = scores[0]

                    for i, score in enumerate(scores):
                        if score < best_score:
                            best_score = score
                            best_i = i

                    self.total_overlap_score = best_score

                    self.rotate_subtree(rotating_atom, parent_atom, math.radians(30 * best_i + 1),
                                        parent_atom.draw.position)

    @staticmethod
    def _find_ring_neighbour(atom: Atom, bond: Bond) -> Atom:
        """
        Return an atom that is inside the same ring as atom and bond, neighbours atom,
            but does not neighbour bond

        Used for mirroring incorrectly positioned stereobonds in rings

        Parameters
        ----------
        atom: Atom instance, must be adjacent to bond and be an atom in a ring
        bond: Bond instance, must neighbour atom and be a bond in a ring

        Returns
        -------
        cyclic_neighbour: Atom instance that is in the same ring as atom, neighbours atom, but does not
            neighbour bond
        """
        # Define the set of rings that both atom 1 and atom 2 are in
        rings: Set[int] = set(bond.atom_1.draw.rings).intersection(set(bond.atom_2.draw.rings))

        cyclic_neighbour: Optional[Atom] = None

        for neighbour in atom.neighbours:
            if len(set(neighbour.draw.rings).intersection(rings)) > 0 and neighbour.draw.is_drawn \
                    and neighbour != bond.get_connected_atom(atom):
                cyclic_neighbour = neighbour
                break

        assert cyclic_neighbour

        return cyclic_neighbour

    def _find_ring_branch_to_flip(self, bond: Bond, neighbours_1: List[Atom],
                                  neighbours_2: List[Atom]) -> Tuple[Optional[Atom],
                                                                      Optional[Tuple[Atom]]]:
        """
        Returns a triplet of atoms, one central and two flanking, which can be mirrored to fix
            stereochemical representation of a stereobond in a ring

        Parameters
        ----------
        bond: Bond instance, stereobond
        neighbours_1: list of drawn Atom instances that neighbour bond.atom_1 but are not bond.atom_2
        neighbours_2: list of drawn Atom instances that neighbour bond.atom_2 but are not bond.atom_1

        Returns
        -------
        central_atom: Optional[Atom instance], atom that will be flipped
        flanking_atoms: Optional[Tuple[Atom, Atom]]. The two atoms about whose midline the central atom will be flipped

        """
        rings: Set[int] = set(bond.atom_1.draw.rings).intersection(set(bond.atom_2.draw.rings))

        # Keeps track of whether or not an incorrectly positioned stereobond in a ring can be resolved
        resolvable: bool = True

        # If one of the atoms only has one drawn neighbour, use that as the central atom to flip around
        if len(neighbours_1) == 1:
            central_atom = bond.atom_1
            flanking_atoms = (neighbours_1[0], bond.atom_2)

        elif len(neighbours_2) == 1:
            central_atom = bond.atom_2
            flanking_atoms = (neighbours_2[0], bond.atom_1)

        # Otherwise, determine which non-ring subtree is the smallest
        else:

            subtree_1_size: Optional[int] = None
            neighbour_1_in_cycle: bool = False
            neighbour_1: Optional[Atom] = None

            subtree_2_size: Optional[int] = None
            neighbour_2_in_cycle: bool = False
            neighbour_2: Optional[Atom] = None

            for neighbour in neighbours_1:
                if len(set(neighbour.draw.rings).intersection(rings)) == 0:

                    subtree_size = self.get_subgraph_size(neighbour, {bond.atom_1})
                    if neighbour_1 and not neighbour_1_in_cycle:
                        if subtree_size < subtree_1_size:
                            subtree_1_size = subtree_size
                            neighbour_1 = neighbour
                    else:
                        subtree_1_size = subtree_size
                        neighbour_1 = neighbour

                    # If a non-cyclic neighbour is selected, always choose that one
                    neighbour_1_in_cycle = False

                else:
                    # Neighbour will only be set as a cyclic one if there is no previous neighbour set
                    if not neighbour_1:
                        neighbour_1 = neighbour
                        neighbour_1_in_cycle = True

            for neighbour in neighbours_2:

                if len(set(neighbour.draw.rings).intersection(rings)) == 0:
                    # If a non-cyclic neighbour is selected, always choose that one
                    neighbour_2_in_cycle = False

                    subtree_size = self.get_subgraph_size(neighbour, {bond.atom_2})
                    if neighbour_2:
                        if subtree_size < subtree_2_size:
                            subtree_2_size = subtree_size
                            neighbour_2 = neighbour
                    else:
                        subtree_2_size = subtree_size
                        neighbour_2 = neighbour

                else:
                    # Neighbour will only be set as a cyclic one if there is no previous neighbour set
                    if not neighbour_2:
                        neighbour_2 = neighbour
                        neighbour_2_in_cycle = True

            assert neighbour_1 and neighbour_2

            # If both atoms have a neighbour that is not in a shared cycle, let the subtree size decide which
            # branch gets flipped.

            if not neighbour_1_in_cycle and not neighbour_2_in_cycle:
                if subtree_2_size > subtree_1_size:
                    central_atom = bond.atom_1
                    ring_atom = self._find_ring_neighbour(bond.atom_1, bond)
                    flanking_atoms = (bond.atom_2, ring_atom)

                else:
                    central_atom = bond.atom_2
                    ring_atom = self._find_ring_neighbour(bond.atom_2, bond)
                    flanking_atoms = (bond.atom_1, ring_atom)

            # Otherwise, prioritise branches that are not in a shared cycle
            elif neighbour_1_in_cycle and not neighbour_2_in_cycle:
                central_atom = bond.atom_2
                ring_atom = self._find_ring_neighbour(bond.atom_2, bond)
                flanking_atoms = (bond.atom_1, ring_atom)

            elif neighbour_2_in_cycle and not neighbour_1_in_cycle:
                central_atom = bond.atom_1
                ring_atom = self._find_ring_neighbour(bond.atom_1, bond)
                flanking_atoms = (bond.atom_2, ring_atom)

            # If both branches contain a shared cycle, the structure is not resolvable
            else:
                central_atom = None
                flanking_atoms = None
                resolvable = False

        if resolvable:
            return central_atom, flanking_atoms
        else:
            return None, None

    def _flip_stereobond_in_ring(self, bond: Bond) -> None:
        """
        Mirror an atom to correct the stereochemistry of a double bond in a ring

        Example:

           C =======                   //
          /                           //
         /              ->           //
        /                   ---------C

        Parameters
        ----------
        bond: Bond instance, a stereochemically restricted double bond

        """
        neighbours_1: List[Atom] = bond.atom_1.get_drawn_neighbours()[:]
        neighbours_2: List[Atom] = bond.atom_2.get_drawn_neighbours()[:]

        neighbours_1.remove(bond.atom_2)
        neighbours_2.remove(bond.atom_1)

        # An atom neighbouring the chiral double bond which will be flipped 'inside' the ring
        central_atom: Optional[Atom] = None

        # Atoms that are members of stereobonds have at least two neighbours: the atom at the other side
        # of the bond, and an atom that determines the cis/trans configuration.

        resolvable: bool = True
        flanking_atoms: List[Atom] = []

        # Check if the neighbouring atoms are adjacent to stereobonds, and if those
        # stereobonds have already been fixed

        neighbours_1_adjacent_to_stereobond = False
        neighbours_2_adjacent_to_stereobond = False

        neighbours_1_adjacent_to_fixed_stereobond = False
        neighbours_2_adjacent_to_fixed_stereobond = False

        for neighbour_1 in neighbours_1:
            if neighbour_1._adjacent_to_stereobond():
                neighbours_1_adjacent_to_stereobond = True
                for bond_1 in neighbour_1.bonds:
                    if bond_1.chiral and bond_1 in self.fixed_chiral_bonds:
                        neighbours_1_adjacent_to_fixed_stereobond = True

        for neighbour_2 in neighbours_2:
            if neighbour_2._adjacent_to_stereobond():
                neighbours_2_adjacent_to_stereobond = True
                for bond_2 in neighbour_2.bonds:
                    if bond_2.chiral and bond_2 in self.fixed_chiral_bonds:
                        neighbours_2_adjacent_to_fixed_stereobond = True

        if not neighbours_1_adjacent_to_stereobond and not neighbours_2_adjacent_to_stereobond:

            central_atom, flanking_atoms = self._find_ring_branch_to_flip(bond, neighbours_1, neighbours_2)
            if not central_atom:
                resolvable = False

        if neighbours_1_adjacent_to_stereobond and not neighbours_2_adjacent_to_stereobond:
            central_atom = bond.atom_2
            ring_neighbour = self._find_ring_neighbour(bond.atom_2, bond)
            flanking_atoms = [bond.atom_1, ring_neighbour]

        elif neighbours_2_adjacent_to_stereobond and not neighbours_1_adjacent_to_stereobond:
            central_atom = bond.atom_1
            ring_neighbour = self._find_ring_neighbour(bond.atom_1, bond)
            flanking_atoms = [bond.atom_2, ring_neighbour]

        elif neighbours_1_adjacent_to_stereobond and neighbours_2_adjacent_to_stereobond:
            if neighbours_1_adjacent_to_fixed_stereobond and not neighbours_2_adjacent_to_fixed_stereobond:
                central_atom = bond.atom_2
                ring_neighbour = self._find_ring_neighbour(bond.atom_2, bond)
                flanking_atoms = [bond.atom_1, ring_neighbour]

            elif neighbours_2_adjacent_to_fixed_stereobond and not neighbours_1_adjacent_to_fixed_stereobond:
                central_atom = bond.atom_1
                ring_neighbour = self._find_ring_neighbour(bond.atom_1, bond)
                flanking_atoms = [bond.atom_2, ring_neighbour]
            elif not neighbours_1_adjacent_to_fixed_stereobond and not neighbours_2_adjacent_to_fixed_stereobond:
                central_atom, flanking_atoms = self._find_ring_branch_to_flip(bond, neighbours_1, neighbours_2)
                if not central_atom:
                    resolvable = False
            else:
                resolvable = False

        if resolvable:
            self._flip_subtree(central_atom, flanking_atoms[0], flanking_atoms[1])

        else:
            if self.options.strict_mode:
                raise DrawingError('chiral bond ring')
            else:
                print("Warning! Cis/trans stereochemistry of cyclic system incorrectly drawn.")
            
    def _flip_subtree(self, root: Atom, atom_1: Atom, atom_2: Atom) -> None:
        """
        Mirror a subtree starting at root about an imaginary line between atom_1 and atom_2

        Parameters
        ----------
        root: Atom instance, the first atom to be flipped. Must be adjacent to atom_1 and atom_2
        atom_1: Atom instance, start point of imaginary mirror line
        atom_2: Atom instance, end point of imaginary mirror line
        """

        for atom in self.traverse_substructure(root, {atom_1, atom_2}):
            atom.draw.position.mirror_about_line(atom_1.draw.position, atom_2.draw.position)
            for anchored_ring in atom.draw.anchored_rings:
                if anchored_ring.center:
                    anchored_ring.center.mirror_about_line(atom_1.draw.position, atom_2.draw.position)

    def _find_clashing_atoms(self) -> List[Tuple[Atom, Atom]]:
        """
        Returns a list of tuples of non-adjacent atoms which are closer than a bond length together
        """
        clashing_atoms: List[Tuple[Atom, Atom]] = []
        for i, atom_1 in enumerate(self.drawn_atoms):
            for j in range(i + 1, len(self.drawn_atoms)):
                atom_2 = self.drawn_atoms[j]
                if not self.structure.bond_exists(atom_1, atom_2):
                    distance = Vector.subtract_vectors(atom_1.draw.position, atom_2.draw.position).get_squared_length()
                    if distance < 0.8 * self.options.bond_length_squared:
                        clashing_atoms.append((atom_1, atom_2))

        return clashing_atoms

    # TODO: Check if implicit hydrogens are excluded
    def _prioritise_chiral_bonds(self, chiral_center: Atom) -> List[Atom]:
        """
        Returns a list of atoms, indicating the order in which atoms adjacent to a chiral center should be prioritised
            for being depicted with a wedged bond

        Parameters
        ----------
        chiral_center: Atom instance, must be a chiral center in the structure
        """

        subtree_1_size = self.get_subgraph_size(chiral_center.neighbours[0], {chiral_center})
        subtree_2_size = self.get_subgraph_size(chiral_center.neighbours[1], {chiral_center})
        subtree_3_size = self.get_subgraph_size(chiral_center.neighbours[2], {chiral_center})

        sizes_and_atoms = [(subtree_1_size, chiral_center.neighbours[0]),
                           (subtree_2_size, chiral_center.neighbours[1]),
                           (subtree_3_size, chiral_center.neighbours[2])]

        if len(chiral_center.neighbours) == 4:
            subtree_4_size = self.get_subgraph_size(chiral_center.neighbours[3], {chiral_center})
            sizes_and_atoms.append((subtree_4_size, chiral_center.neighbours[3]))

        # Sort first by subtree size and then by atom priority
        sizes_and_atoms.sort(key=lambda x: (x[0], ATOM_PROPERTIES.element_to_atomic_nr[x[1].type]))

        options_h: List[Atom] = []
        options: List[Atom] = []
        backup_options_rings: List[Atom] = []
        backup_options_chiral_noring: List[Atom] = []
        backup_options_chiral_ring: List[Atom] = []
        backup_options_chiral_neighbour: List[Atom] = []
        backup_options_chiral_neighbour_ring: List[Atom] = []
        non_options: List[Atom] = []

        for neighbour in [size_and_atom[1] for size_and_atom in sizes_and_atoms]:
            # Priotitise explicit hydrogens
            if neighbour.type == 'H':
                options_h.append(neighbour)
            else:
                other_chiral_centre = False

                for atom in neighbour.neighbours:
                    if atom != chiral_center and atom.chiral:
                        other_chiral_centre = True
                        break

                in_ring = neighbour.in_ring(self.structure)

                if self.structure.bond_lookup[chiral_center][neighbour].type != 'single':
                    non_options.append(neighbour)

                # Then, prioritise atoms not in a ring and not adjacent to a chiral centre

                elif not other_chiral_centre and not in_ring and not neighbour.chiral:
                    options.append(neighbour)

                # Then, prioritise atoms not in a ring and adjacent to a chiral centre

                elif other_chiral_centre and not neighbour.chiral and not in_ring:
                    backup_options_chiral_noring.append(neighbour)

                # Finally, consider some backup options

                elif in_ring and not other_chiral_centre and not neighbour.chiral:
                    backup_options_rings.append(neighbour)

                elif in_ring and other_chiral_centre and not neighbour.chiral:
                    backup_options_chiral_ring.append(neighbour)

                elif not in_ring and neighbour.chiral:
                    backup_options_chiral_neighbour.append(neighbour)

                else:
                    backup_options_chiral_neighbour_ring.append(neighbour)

        priority = options_h + options + backup_options_chiral_noring + \
            backup_options_rings + backup_options_chiral_ring + backup_options_chiral_neighbour + \
            backup_options_chiral_neighbour_ring + non_options

        if chiral_center.type == 'C':
            assert len(priority) == 4

        return priority

    @staticmethod
    def _place_eclipsed_bond(hydrogen: Atom, angles_between_lines: List[float], atom_order: List[Atom],
                             wedge_atom: Atom) -> None:
        """
        Find the conventional position for an eclipsed bond and adjust atom order accordingly

        Parameters
        ----------
        hydrogen: Atom instance, hydrogen to be placed adjacent to chiral center
        angles_between_lines: list of float, with each float representing an angle between two bonds
        atom_order: list of Atom instances, with each atom a non-hydrogen neighbour of a chiral center
        wedge_atom: Atom instance, neighbour of chiral center whose incoming bond is drawn as a wedge
        """
        position: Optional[int] = None

        for i, angle in enumerate(angles_between_lines):

            if round(angle, 3) >= round(math.pi, 3):
                position = i

        if position is not None:
            atom_order.insert(position, hydrogen)
        else:
            wedge_index = atom_order.index(wedge_atom)
            atom_order.insert(wedge_index + 1, hydrogen)

    @staticmethod
    def _reorder(atom_order: List[Atom], wedge_atom: Atom) -> List[Atom]:
        """
        Re-order order in which atoms are drawn around chiral centre such that the list starts with wedge_atom

        Parameters
        ----------
        atom_order: list of Atom instances, neighbours of a chiral center
        wedge_atom: Atom instance, neighbour of chiral center whose incoming bond is drawn as a wedge

        Returns
        -------
        new_order: list of Atom instances, reordered

        """
        first_index = atom_order.index(wedge_atom)
        new_order = atom_order[first_index:] + atom_order[:first_index]
        return new_order

    def move_structure(self, x: float = 0.0, y: float = 0.0) -> None:
        """
        Move the entire structure by vector [x, y]

        Parameters
        ----------
        x: number of pixels to move horizontally
        y: number of pixels to move vertically
        """
        for atom in self.structure.graph:
            if atom.draw.is_drawn:
                atom.draw.position.x += x
                atom.draw.position.y += y
    
    def _determine_chirality(self, chiral_center: Atom) -> Tuple[Bond, str]:
        """
        Return bond to be drawn as wedge as well as the direction of the wedge

        Parameters
        ----------
        chiral_center: Atom instance, a chiral center

        Returns
        -------
        wedge_bond: Bond instance adjacent to the chiral center which is to be drawn as a wedge
        wedge: str, 'front' or 'back', identifies the type of wedge that should be drawn

        """

        # Get all angles of all drawn atoms to the chiral center

        angles_and_atoms = []

        for neighbour in chiral_center.neighbours:
            if neighbour.draw.is_drawn:
                angle = Vector.get_line_angle(chiral_center.draw.position, neighbour.draw.position)
                if angle < 0:
                    angle += 2 * math.pi
                angles_and_atoms.append((angle, neighbour))

        # sort them by angle

        angles_and_atoms.sort(key=lambda x: x[0])

        # get the angles between all adjacent atoms

        angles_between_lines = []
        atom_order = []

        for i, (angle_1, atom) in enumerate(angles_and_atoms):

            angle_2 = angles_and_atoms[i - 1][0]
            if angle_2 > angle_1:
                angle_between_lines = 2 * math.pi - angle_2 + angle_1
            else:
                angle_between_lines = angle_1 - angle_2

            atom_order.append(atom)

            angles_between_lines.append(angle_between_lines)

        # Determine which atoms are the best candidates for receiving an incoming chiral bond
        priority: List[Atom] = self._prioritise_chiral_bonds(chiral_center)

        # This means a hydrogen or lone pair that takes up a sp3 orbital was not drawn
        if len(atom_order) == 3:

            if chiral_center.has_neighbour('H'):
                hydrogen = chiral_center.get_neighbour('H')
                assert not hydrogen.draw.is_drawn
                eclipsed_element = hydrogen
                wedge_atom = priority[1]  # The atom after the implicit hydrogen

            elif chiral_center.lone_pairs:
                lone_pair = chiral_center.lone_pairs[0]
                eclipsed_element = lone_pair
                wedge_atom = priority[0]
            else:
                raise DrawingError("chiral center")

            # Eclipsed elements are NOT drawn, but their implied placement matters in correctly
            # drawing of a chiral center
            self._place_eclipsed_bond(eclipsed_element, angles_between_lines, atom_order, wedge_atom)

        else:
            wedge_atom = priority[0]

        assert len(atom_order) == 4

        atom_order = self._reorder(atom_order, wedge_atom)

        original_order: List[Union[Atom, LonePair]] = chiral_center.neighbours[:]

        # It is possible that a chiral center only has 3 bonds, if a lone pair occupies the fourth sp3 orbital
        if len(original_order) == 3:
            if chiral_center.lone_pairs:
                original_order.append(chiral_center.lone_pairs[0])
            else:
                raise DrawingError("chiral center")

        order_matches_chirality: bool = False

        # Get all 12 list permutations that match the chirality of the original order
        chiral_permutations = get_chiral_permutations(original_order)
        if tuple(atom_order) in chiral_permutations:
            order_matches_chirality = True

        if order_matches_chirality and chiral_center.chiral == 'counterclockwise':
            wedge = 'front'
        elif order_matches_chirality and chiral_center.chiral == 'clockwise':
            wedge = 'back'
        elif not order_matches_chirality and chiral_center.chiral == 'counterclockwise':
            wedge = 'back'
        else:
            wedge = 'front'

        wedge_bond: Bond = self.structure.bond_lookup[wedge_atom][chiral_center]

        return wedge_bond, wedge

    def set_chiral_bonds(self) -> None:
        """
        Find and store wedged bonds in the structure
        """
        for atom in self.structure.graph:
            if atom.chiral:
                bond, wedge = self._determine_chirality(atom)

                self.chiral_bonds.append(bond)
                self.chiral_bond_to_orientation[bond] = (wedge, atom)

    def flip_y_axis(self) -> None:
        """
        Flip a structure across its y-axis
        """
        for atom in self.structure.graph:
            if atom.draw.positioned:
                atom.draw.position.y = -atom.draw.position.y

    def convert_to_int(self) -> None:
        """
        Convert all atom positions to integer values
        """
        for atom in self.structure.graph:
            if atom.draw.positioned:
                atom.draw.position.x = int(atom.draw.position.x)
                atom.draw.position.y = int(atom.draw.position.y)

    def move_to_positive_coords(self) -> None:
        """
        Move the structure to a positive coordinate system
        """
        min_x = 100000000
        min_y = 100000000

        for atom in self.structure.graph:
            if atom.draw.positioned:
                if atom.draw.position.x < min_x:
                    min_x = atom.draw.position.x
                if atom.draw.position.y < min_y:
                    min_y = atom.draw.position.y

        x_translation = abs(min(0, min_x))
        y_translation = abs(min(0, min_y))

        self.move_structure(x_translation + self.options.padding + 1, y_translation + self.options.padding + 1)

    def _draw_text(self, text: Union[str, int], x: float, y: float, font_size: Optional[int] = None,
                   color: str = 'black') -> str:
        """
        Return svg element encoding text

        Parameters
        ----------
        text: text (str for atom types, int for ato counts) to be written in svg format
        x: float, horizontal position (centered) of the text
        y: float, vertical position (centered) of the text
        font_size: Optional[int], if None, use default font size

        Returns
        -------
        svg_text: str, svg text

        """
        if font_size is None:
            font_size = self.options.svg_font_size
        svg_text = f"""<text x="{x}" y="{y}" text-anchor="middle" font-family="{self.options.svg_font}" fill="{color}" font-size="{font_size}"><tspan y="{y}" dy="0.35em">{text}</tspan></text>"""

        return svg_text

    # TODO: Convert structures to a json dict first, then use that dict to draw in matplotlib or directly write to svg
    def draw_svg(self, annotation: Optional[str] = None, numbered_atoms: List[Atom] = None) -> str:
        """
        Returns a str containing svg elements needed for drawing the entire structure

        Parameters
        ----------
        annotation: str, annotation for this structure to group the svg by. Especially useful when
            drawing multiple structures
        numbered_atoms: list of Atom instances, atoms to be labelled with their indices in the drawing

        Returns
        -------
        svg_string: str, containing all svg elements needed for drawing the entire structure. Excludes svg header and
            footer. To include these, call the method 'write_svg' instead.

        """

        self.set_annotation_for_grouping(annotation)

        for ring in self.rings:
            self.set_ring_center(ring)

        for bond_nr, bond in self.structure.bonds.items():
            if bond.atom_1.draw.positioned and bond.atom_2.draw.positioned:
                line = Line(bond.atom_1.draw.position, bond.atom_2.draw.position, bond.atom_1, bond.atom_2)
                midpoint = line.get_midpoint()

                if bond.type == 'single':
                    if bond in self.chiral_bonds:

                        orientation, chiral_center = self.chiral_bond_to_orientation[bond]
                        self._draw_chiral_bond(orientation, chiral_center, line, midpoint)
                    else:
                        self._draw_halflines(line, midpoint)
                elif bond.type in {'double', 'aromatic'}:
                    aromatic = False

                    if bond.type == 'aromatic':
                        aromatic = True

                    if not self._is_terminal(bond.atom_1) and not self._is_terminal(bond.atom_2):
                        self._draw_halflines(line, midpoint)

                        common_ring_numbers = self._get_common_rings(bond.atom_1, bond.atom_2)

                        # One bond is as usual, the other is drawn a little to the side

                        if common_ring_numbers:

                            common_rings = []
                            for ring_nr in common_ring_numbers:
                                common_rings.append(self.get_ring(ring_nr))

                            common_rings.sort(key=lambda x: len(x.members))
                            common_ring = common_rings[0]
                            ring_centre = common_ring.center

                            # If the double bond is in a ring, draw the second line towards the centre of the ring
                            second_line = line.double_line_towards_center(ring_centre, self.options.bond_spacing,
                                                                          self.options.double_bond_length)
                            second_line_midpoint = second_line.get_midpoint()
                            self._draw_halflines_double(second_line, second_line_midpoint, aromatic)

                        else:

                            bond_neighbours = bond.atom_1.drawn_neighbours + bond.atom_2.drawn_neighbours
                            assert bond_neighbours

                            vectors = [atom.draw.position for atom in bond_neighbours]
                            # If the bond is not in a ring, draw the second line where you have more neighbouring bonds
                            gravitational_point = Vector.get_average(vectors)
                            second_line = line.double_line_towards_center(gravitational_point,
                                                                          self.options.bond_spacing,
                                                                          self.options.double_bond_length)
                            second_line_midpoint = second_line.get_midpoint()
                            self._draw_halflines_double(second_line, second_line_midpoint, aromatic)

                    else:
                        if self._is_terminal(bond.atom_1) and self._is_terminal(bond.atom_2):
                            dummy_1 = Vector(bond.atom_1.draw.position.x + 1, bond.atom_1.draw.position.y + 1)
                            dummy_2 = Vector(bond.atom_1.draw.position.x - 1, bond.atom_1.draw.position.y - 1)
                            double_bond_line_1 = line.double_line_towards_center(dummy_1,
                                                                                 self.options.bond_spacing / 2.0,
                                                                                 self.options.double_bond_length)
                            double_bond_line_1_midpoint = double_bond_line_1.get_midpoint()
                            double_bond_line_2 = line.double_line_towards_center(dummy_2,
                                                                                 self.options.bond_spacing / 2.0,
                                                                                 self.options.double_bond_length)
                            double_bond_line_2_midpoint = double_bond_line_2.get_midpoint()

                            self._draw_halflines_double(double_bond_line_1, double_bond_line_1_midpoint)
                            self._draw_halflines_double(double_bond_line_2, double_bond_line_2_midpoint, aromatic)

                        else:

                            if self._is_terminal(bond.atom_1):
                                terminal_atom = bond.atom_1
                                branched_atom = bond.atom_2
                            else:
                                terminal_atom = bond.atom_2
                                branched_atom = bond.atom_1

                            if len(branched_atom.drawn_neighbours) >= 3:
                                closest_two = self.get_sorted_distances_from_list(terminal_atom,
                                                                                  branched_atom.drawn_neighbours)
                                closest_atom_1 = closest_two[0][1]
                                closest_atom_2 = closest_two[1][1]

                                line = Line(terminal_atom.draw.position, branched_atom.draw.position, terminal_atom,
                                            branched_atom)

                                double_bond_line_1, double_bond_line_2 = line.get_perpendicular_lines(
                                    self.options.bond_spacing / 2.0)
                                terminal_atom_pos_1 = double_bond_line_1.get_atom_coords(terminal_atom)
                                terminal_atom_pos_2 = double_bond_line_2.get_atom_coords(terminal_atom)

                                closest_atom_to_pos_1 = terminal_atom_pos_1.get_closest_atom(closest_atom_1,
                                                                                             closest_atom_2)
                                closest_atom_to_pos_2 = terminal_atom_pos_2.get_closest_atom(closest_atom_1,
                                                                                             closest_atom_2)

                                bond_1_line = Line(branched_atom.draw.position, closest_atom_to_pos_1.draw.position,
                                                   branched_atom, closest_atom_to_pos_1)
                                bond_2_line = Line(branched_atom.draw.position, closest_atom_to_pos_2.draw.position,
                                                   branched_atom, closest_atom_to_pos_2)

                                double_bond_line_1_midpoint = double_bond_line_1.get_midpoint()
                                double_bond_line_2_midpoint = double_bond_line_2.get_midpoint()

                                intersection_1 = double_bond_line_1.find_intersection(bond_1_line)
                                intersection_2 = double_bond_line_2.find_intersection(bond_2_line)

                                if terminal_atom.draw.position.x > branched_atom.draw.position.x:
                                    # check for parallel lines
                                    if intersection_1 and intersection_1.x < 100000 and intersection_1.y < 100000:
                                        double_bond_line_1.point_1 = intersection_1
                                    if intersection_2 and intersection_2.x < 100000 and intersection_2.y < 100000:
                                        double_bond_line_2.point_1 = intersection_2

                                else:
                                    # check for parallel lines
                                    if intersection_1 and intersection_1.x < 100000 and intersection_1.y < 100000:
                                        double_bond_line_1.point_2 = intersection_1
                                    if intersection_2 and intersection_2.x < 100000 and intersection_2.y < 100000:
                                        double_bond_line_2.point_2 = intersection_2

                                self._draw_halflines(double_bond_line_1, double_bond_line_1_midpoint)
                                self._draw_halflines(double_bond_line_2, double_bond_line_2_midpoint, aromatic)

                            else:
                                self._draw_halflines(line, midpoint)

                                bond_neighbours = bond.atom_1.drawn_neighbours + bond.atom_2.drawn_neighbours
                                if bond_neighbours:
                                    vectors = [atom.draw.position for atom in bond_neighbours]
                                    gravitational_point = Vector.get_average(vectors)
                                    second_line = line.get_parallel_line(gravitational_point,
                                                                         self.options.bond_spacing)
                                    second_line_midpoint = second_line.get_midpoint()
                                    self._draw_halflines(second_line, second_line_midpoint, aromatic)
                                else:
                                    print("Shouldn't happen!")

                elif bond.type == 'triple':
                    self._draw_halflines(line, midpoint)
                    line_1, line_2 = line.get_parallel_lines(self.options.bond_spacing)
                    line_1_midpoint = line_1.get_midpoint()
                    line_2_midpoint = line_2.get_midpoint()
                    self._draw_halflines(line_1, line_1_midpoint)
                    self._draw_halflines(line_2, line_2_midpoint)

        for atom in self.structure.graph:
            if atom.draw.positioned:
                svg_text = ''
                svg_h_text = ''
                svg_charge_text = ''
                svg_h_count_text = ''

                if atom.type != 'C' or atom.draw.draw_explicit or atom.charge:
                    if atom.type == 'C' and not atom.charge:
                        svg_text = self._draw_text('.', atom.draw.position.x, atom.draw.position.y - 2,
                                                   color=atom.draw.colour)
                    else:
                        svg_text = self._draw_text(atom.type, atom.draw.position.x, atom.draw.position.y,
                                                   color=atom.draw.colour)

                        # TODO: Make this possible in svg writing
                        # text = self.set_r_group_indices_subscript(atom.type)

                orientation = self._get_hydrogen_text_orientation(atom)

                # Swap up-down orientation due to swapped svg coordinate system
                if orientation == 'H_below_atom':
                    orientation = 'H_above_atom'
                elif orientation == 'H_above_atom':
                    orientation = 'H_below_atom'

                if atom.type != 'C' or atom.draw.draw_explicit or atom.charge:

                    hydrogen_count = 0

                    for neighbour in atom.neighbours:
                        if neighbour.type == 'H' and not neighbour.draw.is_drawn:
                            hydrogen_count += 1

                    h_x = atom.draw.position.x
                    h_y = atom.draw.position.y
                    h_subscript_x = atom.draw.position.x
                    h_subscript_y = atom.draw.position.y + 3
                    charge_x = atom.draw.position.x + 6
                    charge_y = atom.draw.position.y - 3

                    if orientation == 'H_below_atom':
                        h_y = atom.draw.position.y + 7
                        h_subscript_x += 2
                        h_subscript_y += 10
                    elif orientation == 'H_above_atom':
                        h_y = atom.draw.position.y - 7
                        h_subscript_x += 2
                        h_subscript_y -= 4
                    elif orientation == 'H_before_atom':
                        if hydrogen_count > 1:
                            h_x -= 10
                            h_subscript_x -= 5
                        elif hydrogen_count == 1:
                            h_x -= 6
                    else:
                        h_x += len(atom.type) * 6
                        h_subscript_x += len(atom.type) * 6 + 5

                    if abs(atom.charge):
                        if hydrogen_count and orientation == 'H_after_atom':
                            if hydrogen_count > 1:
                                charge_x += len(atom.type) * 6 + 7
                            else:
                                charge_x += len(atom.type) * 6 + 2

                    h_pos = Vector(h_x, h_y)
                    h_subscript_pos = Vector(h_subscript_x, h_subscript_y)
                    charge_pos = Vector(charge_x, charge_y)

                    if hydrogen_count:
                        svg_h_text = self._draw_text('H', h_pos.x, h_pos.y,
                                                     color=atom.draw.colour)
                        if hydrogen_count > 1:
                            svg_h_count_text = self._draw_text(hydrogen_count, h_subscript_pos.x, h_subscript_pos.y,
                                                               font_size=self.options.svg_font_size_small,
                                                               color=atom.draw.colour)
                    if atom.charge:
                        if atom.charge > 0:
                            charge_symbol = '+'
                        else:
                            charge_symbol = '-'

                        if abs(atom.charge) > 1:
                            charge_text = f"{abs(atom.charge)}{charge_symbol}"
                        else:
                            charge_text = charge_symbol

                        svg_charge_text = self._draw_text(charge_text, charge_pos.x, charge_pos.y,
                                                          font_size=self.options.svg_font_size_small,
                                                          color=atom.draw.colour)

                if svg_text or svg_h_text or svg_charge_text or svg_h_count_text:
                    if self.structure_id:
                        text_group = f'<g id="atom_{atom}_{self.structure_id}_text">\n'
                    else:
                        text_group = f'<g id="atom_{atom}_text">\n'
                    if svg_text:
                        text_group += svg_text
                        text_group += '\n'
                    if svg_h_text:
                        text_group += svg_h_text
                        text_group += '\n'
                    if svg_charge_text:
                        text_group += svg_charge_text
                        text_group += '\n'
                    if svg_h_count_text:
                        text_group += svg_h_count_text
                        text_group += '\n'
                    text_group += '</g>'
                    self._add_svg_element(text_group, atom)
                    
        if numbered_atoms:
            self.show_atom_numbers(numbered_atoms)

        svg = self._assemble_svg()
        return svg

    def write_svg(self, out_file: str, annotation: Optional[str] = None, numbered_atoms: List[Atom] = None) -> None:
        """
        Write svg of structure to an output file

        Parameters
        ----------
        out_file: str, path to output file
        annotation: str, annotation for the structure
        numbered_atoms: list of Atom instances to be labelled with their indices in the drawing
        """

        self.flip_y_axis()
        self.move_to_positive_coords()
        self.convert_to_int()

        min_x = 100000000
        max_x = -100000000
        min_y = 100000000
        max_y = -100000000

        for atom in self.structure.graph:
            if atom.draw.positioned:
                if atom.draw.position.x < min_x:
                    min_x = atom.draw.position.x
                if atom.draw.position.y < min_y:
                    min_y = atom.draw.position.y
                if atom.draw.position.x > max_x:
                    max_x = atom.draw.position.x
                if atom.draw.position.y > max_y:
                    max_y = atom.draw.position.y

        width = max_x - min_x + 2 * self.options.padding
        height = max_y - min_y + 2 * self.options.padding

        x1 = max(0, min_x - self.options.padding)
        y1 = max(0, min_y - self.options.padding)
        x2 = max_x + self.options.padding
        y2 = max_y + self.options.padding

        svg_string = f"""<svg width="{width}" height="{height}" viewBox="{x1} {y1} {x2} {y2}" xmlns="http://www.w3.org/2000/svg">"""
        svg_string += self.svg_style
        svg_string += self.draw_svg(annotation=annotation, numbered_atoms=numbered_atoms)
        svg_string += "</svg>"

        with open(out_file, 'w') as out:
            out.write(svg_string)

    def _assemble_svg(self) -> str:
        """
        Assemble and return string containing svg elements from dictionary of svg elements

        Returns
        -------
        svg_string: str containing all svg elements required for visualising the structure

        """

        svg_string = ''

        for svg_group, atom_to_elements in self.svg_groups.items():
            svg_string += f"""<g id="{svg_group}">\n"""

            for atom, elements in atom_to_elements.items():
                svg_string += f"""<g id="{atom}">\n"""

                for element in elements:
                    svg_string += element
                    svg_string += '\n'
                svg_string += "</g>\n"
            svg_string += "</g>\n"

        return svg_string

    # TODO: This should be refactored. Better to always only do the positioning, and then specifically call 'draw'
    def draw(self, coords_only: bool = False) -> None:
        """
        Position the atoms of the structure in a plane. If coords_only is False, also draw the structure in matplotlib

        Parameters
        ----------
        coords_only: bool, only position the atoms if True, also draw them in matplotlib if False
        """

        if not self.options.draw_hydrogens:
            self.hide_hydrogens()

        self.get_atom_nr_to_atom()
        self.define_rings()

        if not self.multiple:
            self._process_structure()
            self.set_chiral_bonds()
            if not coords_only:
                self.plot_structure()
        else:
            self.restore_ring_information()

    @staticmethod
    def _get_hydrogen_text_orientation(atom: Atom) -> str:
        """
        Returns the optimal orientation of the text depicting hydrogen atoms to achieve minimal overlap

        Parameters
        ----------
        atom: Atom instance, atom that has neighbouring hydrogen atoms

        Returns
        -------
        orientation: str, 'H_before_atom', 'H_after_atom', 'H_above_atom' or 'H_below_atom'
        """
        four_positions = [Vector(atom.draw.position.x, atom.draw.position.y + 3),
                          Vector(atom.draw.position.x, atom.draw.position.y - 3),
                          Vector(atom.draw.position.x + 3, atom.draw.position.y),
                          Vector(atom.draw.position.x - 3, atom.draw.position.y)]

        positions_to_angles: List[List[float]] = [[], [], [], []]

        for neighbour in atom.drawn_neighbours:
            for i, position in enumerate(four_positions):
                angle = Vector.get_angle_between_vectors(position, neighbour.draw.position, atom.draw.position)
                positions_to_angles[i].append(angle)

        if not positions_to_angles[0]:
            orientation = 'H_after_atom'
        elif min(positions_to_angles[2]) > 1.57078 or min(positions_to_angles[3]) > 1.57078:
            if min(positions_to_angles[2]) >= min(positions_to_angles[3]):
                orientation = 'H_after_atom'

            else:
                orientation = 'H_before_atom'
        else:
            smallest_angles = [min(angles) for angles in positions_to_angles]

            angle: float = 0.0
            position = None

            for j, largest_angle in enumerate(smallest_angles):
                if largest_angle > angle:
                    angle = largest_angle
                    position = j

            if position == 0:
                orientation = 'H_above_atom'
            elif position == 1:
                orientation = 'H_below_atom'
            elif position == 2:
                orientation = 'H_after_atom'
            elif position == 3:
                orientation = 'H_before_atom'
            else:
                orientation = 'H_after_atom'

        return orientation

    @staticmethod
    def _get_common_rings(atom_1: Atom, atom_2: Atom) -> Optional[List[int]]:
        """
        Returns a list of ring identifiers if atom_1 and atom_2 share one or more rings, None otherwise

        Parameters
        ----------
        atom_1: Atom instance
        atom_2: Atom instance

        Returns
        -------
        joined_rings: list of int, with each int a ring identifier
        """
        if atom_1.draw.rings and atom_2.draw.rings:
            joined_rings = list(set(atom_1.draw.rings).intersection(set(atom_2.draw.rings)))
            return joined_rings

        return None

    @staticmethod
    def _plot_chiral_bond_front(polygon: List[Vector], ax: Axes, color: str = 'black') -> None:
        """
        Plot a forward-facing chiral bond, visualised as a solid-coloured wedge, in matplotlib

        Parameters
        ----------
        polygon: list of Vector instances, with each vector one corner of the wedge.
            Must have length 3 for the first half of the wedge, and length 4 for the second half of the wedge
        ax: matplotlib Axes instance to plot in
        color: str, denoting colour of wedge. Default: black
        """
        x: List[float] = []
        y: List[float] = []
        for point in polygon:
            x.append(point.x)
            y.append(point.y)

        ax.fill(x, y, color=color)

    def _plot_chiral_bond_back(self, lines: List[Line], ax: Axes, color: str = 'black') -> None:
        """
        Plot a backward-facing chiral bond, visualised as a series of thin lines, in matplotlib

        Parameters
        ----------
        lines: list of Line instances, with each line one of the lines to be drawn.
        ax: matplotlib Axes instance to plot in
        color: str, denoting colour of wedge. Default: black
        """
        for line in lines:
            self._plot_line(line, ax, color=color)

    def _plot_chiral_bond(self, orientation: str, chiral_center: Atom, line: Line, ax: Axes, midpoint: Vector) -> None:
        """
        Plot a chiral bond in matplotlib

        Parameters
        ----------
        orientation: str, 'front' for a forward-facing chiral bond, 'back' for a backward-facing one
        chiral_center: Atom instance, must be a chiral center
        line: Line instance, line denoting the position of a bond
        ax: matplotlib Axes instance
        midpoint: Vector instance, denoting the midpoint on the line. Used to split the line in two before drawing
            such that each bond half can be coloured separately
        """
        halflines = line.divide_in_two(midpoint)
        for halfline in halflines:
            truncated_line = halfline.get_truncated_line(self.options.short_bond_length)
            if orientation == 'front':
                bond_wedge = truncated_line.get_bond_wedge_front(self.options.chiral_bond_width, chiral_center)
                self._plot_chiral_bond_front(bond_wedge, ax, color=halfline.atom.draw.colour)
            else:
                bond_lines = halfline.get_bond_wedge_back(self.options.chiral_bond_width, chiral_center)
                self._plot_chiral_bond_back(bond_lines, ax, color=halfline.atom.draw.colour)

    @staticmethod
    def _draw_chiral_bond_front(polygon: List[Vector], color: str = 'black') -> str:
        """
        Create an SVG element for half of a forward-facing chiral bond, visualised as a solid-coloured wedge

        Parameters
        ----------
        polygon: list of Vector instances, with each vector one corner of the wedge.
            Must have length 3 for the first half of the wedge, and length 4 for the second half of the wedge
        color: str, denoting colour of wedge. Default: black
        """
        svg_string = '<polygon points="'
        polygon_strings = []
        for point in polygon:
            polygon_strings.append(f"{point.x},{point.y}")
        polygon_string = ' '.join(polygon_strings)
        svg_string += polygon_string + '"'
        if color != 'black':
            svg_string += f' fill="{color}"'
        svg_string += ' />'

        return svg_string

    @staticmethod
    def _draw_chiral_bond_back(lines: List[Line], color: str = 'black') -> str:
        """
        Create an SVG element for half of a backward-facing chiral bond, visualised as a series of thin lines

        Parameters
        ----------
        lines: list of Line instances, with each line one of the lines to be drawn.
        color: str, denoting colour of wedge. Default: black
        """
        svg_string = f'<g id="chiral-bond-back" stroke="{color}" stroke-width="{1}">\n  '
        line_strings = []
        for line in lines:
            line_string = f'<line x1="{line.point_1.x}" y1="{line.point_1.y}" x2="{line.point_2.x}" y2="{line.point_2.y}" />'
            line_strings.append(line_string)

        line_string = '\n  '.join(line_strings)
        svg_string += line_string
        svg_string += '\n</g>'

        return svg_string

    def set_structure_id(self, structure_id: str) -> None:
        """
        Set a structure identifier by which SVG elements will be grouped

        Parameters
        ----------
        structure_id: str, structure identifier

        """
        self.structure_id = structure_id

    def _add_svg_element(self, svg_element: str, atom: Atom) -> None:
        """
        Store an SVG element representing a half-bond or atom within an atom group

        Parameters
        ----------
        svg_element: str, pre-computed SVG element
        atom: Atom instance, associated with the SVG element
        """

        if self.annotation is not None and self.annotation in atom.annotations.annotations:
            annotation_value = atom.annotations.get_annotation(self.annotation)
        else:
            annotation_value = 'unlabeled'

        if self.structure_id:
            annotation_value = f"{annotation_value}_{self.structure_id}"

        if annotation_value not in self.svg_groups:
            self.svg_groups[annotation_value] = {}

        atom_label = f'atom_{atom}'

        if self.structure_id:
            atom_label = f'atom_{atom}_{self.structure_id}'

        if atom_label not in self.svg_groups[annotation_value]:
            self.svg_groups[annotation_value][atom_label] = []

        self.svg_groups[annotation_value][atom_label].append(svg_element)

    def _draw_chiral_bond(self, orientation: str, chiral_center: Atom, line: Line, midpoint: Vector) -> None:
        """
        Add an SVG element representing a chiral bond to the drawing

        Parameters
        ----------
        orientation: str, 'front' for a forward-facing chiral bond, 'back' for a backward-facing one
        chiral_center: Atom instance, must be a chiral center
        line: Line instance, line denoting the position of a bond
        midpoint: Vector instance, denoting the midpoint on the line. Used to split the line in two before drawing
            such that each bond half can be coloured separately
        """
        halflines = line.divide_in_two(midpoint)
        for halfline in halflines:
            truncated_line = halfline.get_truncated_line(self.options.short_bond_length)
            if orientation == 'front':
                bond_wedge = truncated_line.get_bond_wedge_front(self.options.chiral_bond_width, chiral_center)
                svg_polygon = self._draw_chiral_bond_front(bond_wedge, color=halfline.atom.draw.colour)
                self._add_svg_element(svg_polygon, halfline.atom)

            elif orientation == 'back':
                bond_lines = halfline.get_bond_wedge_back(self.options.chiral_bond_width, chiral_center)
                svg_lines = self._draw_chiral_bond_back(bond_lines, color=halfline.atom.draw.colour)
                self._add_svg_element(svg_lines, halfline.atom)
            else:
                raise ValueError(f"Unrecognised chiral bond orientation: {orientation}.")

    def _draw_halflines(self, line: Line, midpoint: Vector, aromatic=False) -> None:
        """
        Add an SVG element representing one line of a bond to the drawing

        Parameters
        ----------
        line: Line instance, line denoting the position of one line of a bond
        midpoint: Vector instance, denoting the midpoint on the line. Used to split the line in two before drawing
            such that each bond half can be coloured separately
        aromatic: bool, if True, draw a dashed line; if False, draw a solid line
        """

        halflines = line.divide_in_two(midpoint)
        for halfline in halflines:
            truncated_line = halfline.get_truncated_line(self.options.short_bond_length)
            if not aromatic:
                svg_line = self._draw_line(truncated_line, color=halfline.atom.draw.colour)
            else:
                svg_line = self._draw_dashed_line(truncated_line, color=halfline.atom.draw.colour)
            self._add_svg_element(svg_line, halfline.atom)

    # TODO: Can we integrate this function with the one above?
    def _draw_halflines_double(self, line: Line, midpoint: Vector, aromatic=False) -> None:
        """
        Add an SVG element representing one line of a pre-truncated double bond to the drawing

        Parameters
        ----------
        line: Line instance, line denoting the position of one line of a bond. This line has been pre-truncated
        midpoint: Vector instance, denoting the midpoint on the line. Used to split the line in two before drawing
            such that each bond half can be coloured separately
        aromatic: bool, if True, draw a dashed line; if False, draw a solid line
        """

        halflines = line.divide_in_two(midpoint)
        for halfline in halflines:
            if not aromatic:
                svg_line = self._draw_line(halfline, color=halfline.atom.draw.colour)
            else:
                svg_line = self._draw_dashed_line(halfline, color=halfline.atom.draw.colour)
            self._add_svg_element(svg_line, halfline.atom)

    @staticmethod
    def _draw_dashed_line(line: Line, color: str = 'black') -> str:
        """
        Create an SVG element representing the dashed line of half of an aromatic bond

        Parameters
        ----------
        line: Line instance, line denoting the position of half of the aromatic line of a bond
        color: str, colour of the line

        Returns
        ----------
        svg_line: str, SVG element depicting half of an aromatic bond
        """
        # Create this loop to limit svg size: default colour is black
        if color != 'black':
            svg_line = f'<line x1="{line.point_1.x}" y1="{line.point_1.y}" x2="{line.point_2.x}" y2="{line.point_2.y}" stroke="{color}" stroke-dasharray="3,2"/>'
        else:
            svg_line = f'<line x1="{line.point_1.x}" y1="{line.point_1.y}" x2="{line.point_2.x}" y2="{line.point_2.y}" stroke-dasharray="3,2"/>'
        return svg_line

    @staticmethod
    def _draw_line(line: Line, color: str = 'black') -> str:
        """
        Create an SVG element representing the line of half of a bond

        Parameters
        ----------
        line: Line instance, line denoting the position of half of the line of a bond
        color: str, colour of the line

        Returns
        ----------
        svg_line: str, SVG element depicting half of a bond
        """
        # Create this loop to limit svg size: default colour is black
        svg_line = f'<line x1="{line.point_1.x}" y1="{line.point_1.y}" x2="{line.point_2.x}" y2="{line.point_2.y}" stroke-width="1" stroke="{color}"/>'
        return svg_line

    def _plot_halflines(self, line: Line, ax: Axes, midpoint: Vector, aromatic: bool = False) -> None:
        """
        Plot one line of a bond in matplotlib

        Parameters
        ----------
        line: Line instance, line denoting the position of one line of a bond
        midpoint: Vector instance, denoting the midpoint on the line. Used to split the line in two before drawing
            such that each bond half can be coloured separately
        ax: matplotlib Axes instance, canvas to draw on
        aromatic: bool, if True, draw a dashed line; if False, draw a solid line
        """

        halflines = line.divide_in_two(midpoint)
        for halfline in halflines:
            truncated_line = halfline.get_truncated_line(self.options.short_bond_length)
            self._plot_line(truncated_line, ax, color=halfline.atom.draw.colour, aromatic=aromatic)

    def _plot_halflines_double(self, line: Line, ax: Axes, midpoint: Vector, aromatic: bool = False) -> None:
        """
        Plot one line of a pre-truncated double bond in matplotlib

        Parameters
        ----------
        line: Line instance, line denoting the position of one line of a bond. This line has been pre-truncated
        midpoint: Vector instance, denoting the midpoint on the line. Used to split the line in two before drawing
            such that each bond half can be coloured separately
        ax: matplotlib Axes instance
        aromatic: bool, if True, draw a dashed line; if False, draw a solid line
        """
        halflines = line.divide_in_two(midpoint)
        for halfline in halflines:
            self._plot_line(halfline, ax, color=halfline.atom.draw.colour, aromatic=aromatic)

    def _plot_line(self, line: Line, ax: Axes, color: str = 'black', aromatic: bool = False) -> None:
        """
        Plot the line of half of a bond in matplotlib

        Parameters
        ----------
        line: Line instance, line denoting the position of half of the line of a bond
        ax: matplotlib Axes instance, canvas to draw to
        color: str, colour of the line
        aromatic: bool, draw dashed line if True, solid line otherwise

        Returns
        ----------
        svg_line: str, SVG element depicting half of a bond
        """
        if not aromatic:
            ax.plot([line.point_1.x, line.point_2.x],
                    [line.point_1.y, line.point_2.y], color=color, linewidth=self.options.bond_thickness)
        else:
            ax.plot([line.point_1.x, line.point_2.x],
                    [line.point_1.y, line.point_2.y], color=color, linewidth=self.options.bond_thickness,
                    linestyle='dashed', dashes=(2, 1))

    @staticmethod
    def get_image_as_array() -> np.ndarray:
        """
        Returns numpy array representing a pixelated version of the drawing
        """
        # Return image as np.ndarray that represents RGB image
        canvas = plt.gca().figure.canvas
        canvas.draw()
        image = np.frombuffer(canvas.tostring_rgb(), dtype='uint8')
        image = image.reshape(canvas.get_width_height()[::-1] + (3,))
        plt.close('all')
        return image

    @staticmethod
    def save_svg_matplotlib(out_file: str) -> None:
        """
        Save matplotlib canvas to an SVG file

        Parameters
        ----------
        out_file: str, path to output SVG
        """
        if out_file.endswith('.svg'):
            pass
        else:
            out_file += '.svg'
        plt.savefig(out_file)
        plt.clf()
        plt.close(plt.gcf())
        plt.close('all')

    @staticmethod
    def get_svg_string_matplotlib() -> str:
        """
        Returns str, SVG string of matplotlib canvas
        """
        svg_string = StringIO()
        plt.savefig(svg_string, format='svg')
        svg = svg_string.getvalue()

        plt.clf()
        plt.close(plt.gcf())
        plt.close('all')
        return svg

    @staticmethod
    def save_png_matplotlib(out_file: str) -> None:
        """
        Save matplotlib canvas to a PNG file

        Parameters
        ----------
        out_file: str, path to output PNG
        """
        if out_file.endswith('.png'):
            pass
        else:
            out_file += '.png'
        plt.savefig(out_file)
        plt.clf()
        plt.close()

    @staticmethod
    def show_molecule() -> None:
        """
        Show canvas in matplotlib GUI
        """
        plt.show()
        plt.clf()
        plt.close()

    @staticmethod
    def _chiral_bond_drawn_correctly(bond: Bond) -> bool:
        """
        Returns True if a chiral double bond is represented correctly in the drawing, False otherwise

        Parameters
        ----------
        bond: Bond instance, must be a chiral double bond
        """
        assert bond.chiral

        must_be_fixed = False

        for neighbour_1 in bond.atom_1.neighbours:
            if neighbour_1 != bond.atom_2:
                for neighbour_2 in bond.atom_2.neighbours:
                    if neighbour_2 != bond.atom_1:
                        if neighbour_1.draw.is_drawn and neighbour_2.draw.is_drawn:

                            placement_1 = Vector.get_position_relative_to_line(bond.atom_1.draw.position,
                                                                               bond.atom_2.draw.position,
                                                                               neighbour_1.draw.position)
                            placement_2 = Vector.get_position_relative_to_line(bond.atom_1.draw.position,
                                                                               bond.atom_2.draw.position,
                                                                               neighbour_2.draw.position)

                            orientation = bond.chiral_dict[neighbour_1][neighbour_2]

                            if orientation == 'cis':
                                if placement_1 != placement_2:
                                    must_be_fixed = True
                            else:
                                if placement_1 == placement_2:
                                    must_be_fixed = True

        if must_be_fixed:
            return False

        else:
            return True

    def _fix_chiral_bond(self, double_bond: Bond) -> None:
        """
        Correct the position of an incorrectly represented chiral double bond

        Parameters
        ----------
        double_bond: Bond instance, chiral double bond

        """
        # If the bond is in a ring, mirror one of the atoms adjacent to the bond to the inside of the ring

        if len(double_bond.atom_1.draw.rings) and len(double_bond.atom_2.draw.rings) and \
                len(set(double_bond.atom_1.draw.rings).intersection(set(double_bond.atom_2.draw.rings))) >= 1:
            self._flip_stereobond_in_ring(double_bond)

        # Otherwise, rotate subtrees in the structure around adjacent bonds

        else:
            if len(double_bond.atom_1.draw.rings) > 0:
                parent_atom = double_bond.atom_2
                root_atom = double_bond.atom_1

            else:
                parent_atom = double_bond.atom_1
                root_atom = double_bond.atom_2

            neighbours = parent_atom.drawn_neighbours[:]
            neighbours.remove(root_atom)

            if len(neighbours) == 1:
                neighbour = neighbours[0]
                self._flip_subtree(neighbour, root_atom, parent_atom)

            # Only need to flip once if both neighbours are in the same ring,
            # as then the neighbours occur in the same subtree

            elif len(neighbours) == 2 and len(
                    set(neighbours[0].draw.rings).intersection(set(neighbours[1].draw.rings))) >= 1:

                self._flip_subtree(neighbours[0], root_atom, parent_atom)

            elif len(neighbours) == 2:
                neighbour_1 = neighbours[0]
                neighbour_2 = neighbours[1]

                self._flip_subtree(neighbour_1, root_atom, parent_atom)
                self._flip_subtree(neighbour_2, root_atom, parent_atom)

            self.fixed_chiral_bonds.add(double_bond)

    def _fix_chiral_bonds_in_rings(self) -> None:
        """
        Iterate over all sequences of alternating chiral double bonds in rings and correct them in order
        """
        double_bond_sequences = self.structure.find_double_bond_sequences()

        for double_bond_sequence in double_bond_sequences:
            for double_bond in double_bond_sequence:
                chirality_correct = self._chiral_bond_drawn_correctly(double_bond)
                if chirality_correct:
                    self.fixed_chiral_bonds.add(double_bond)
                else:
                    self._fix_chiral_bond(double_bond)

        for bond in self.structure.bonds.values():
            if bond.type == 'double' and bond.chiral and bond not in self.fixed_chiral_bonds:
                chirality_correct = self._chiral_bond_drawn_correctly(bond)
                if chirality_correct:
                    self.fixed_chiral_bonds.add(bond)
                else:
                    self._fix_chiral_bond(bond)

    # TODO: Refactor this such that a dictionary of all drawn components is created first
    def plot_structure(self) -> None:
        """
        Plot the positioned atoms and bonds to a matplotlib canvas
        """

        # Find the plotting dimensions of the molecule such that the canvas can be scaled to fit the molecule

        min_x = 100000000
        max_x = -100000000
        min_y = 100000000
        max_y = -100000000

        for atom in self.structure.graph:
            if atom.draw.positioned:
                if atom.draw.position.x < min_x:
                    min_x = atom.draw.position.x
                if atom.draw.position.y < min_y:
                    min_y = atom.draw.position.y
                if atom.draw.position.x > max_x:
                    max_x = atom.draw.position.x
                if atom.draw.position.y > max_y:
                    max_y = atom.draw.position.y

        height: float = max_y - min_y
        width: float = max_x - min_x

        fig, ax = plt.subplots(figsize=((width + 2 * self.options.padding) / 50.0,
                                        (height + 2 * self.options.padding) / 50.0), dpi=100)

        ax.set_aspect('equal', adjustable='box')
        ax.axis('off')

        ax.set_xlim([min_x - self.options.padding, max_x + self.options.padding])
        ax.set_ylim([min_y - self.options.padding, max_y + self.options.padding])
        plt.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)

        params = {'mathtext.default': 'regular'}
        plt.rcParams.update(params)

        for ring in self.rings:
            self.set_ring_center(ring)

        for bond_nr, bond in self.structure.bonds.items():
            if bond.atom_1.draw.positioned and bond.atom_2.draw.positioned:
                line = Line(bond.atom_1.draw.position, bond.atom_2.draw.position, bond.atom_1, bond.atom_2)
                midpoint = line.get_midpoint()

                if bond.type == 'single':
                    if bond in self.chiral_bonds:
                        orientation, chiral_center = self.chiral_bond_to_orientation[bond]
                        self._plot_chiral_bond(orientation, chiral_center, line, ax, midpoint)
                    else:
                        self._plot_halflines(line, ax, midpoint)
                elif bond.type in {'double', 'aromatic'}:
                    aromatic = False
                    if bond.type == 'aromatic':
                        aromatic = True
                    
                    if not self._is_terminal(bond.atom_1) and not self._is_terminal(bond.atom_2):
                        self._plot_halflines(line, ax, midpoint)

                        common_ring_numbers = self._get_common_rings(bond.atom_1, bond.atom_2)

                        if common_ring_numbers:
                            common_rings = []
                            for ring_nr in common_ring_numbers:
                                common_rings.append(self.get_ring(ring_nr))

                            common_rings.sort(key=lambda x: len(x.members))
                            common_ring = common_rings[0]
                            ring_centre = common_ring.center
                            second_line = line.double_line_towards_center(ring_centre, self.options.bond_spacing, self.options.double_bond_length)
                            second_line_midpoint = second_line.get_midpoint()
                            self._plot_halflines_double(second_line, ax, second_line_midpoint, aromatic)

                        else:
                            bond_neighbours = bond.atom_1.drawn_neighbours + bond.atom_2.drawn_neighbours
                            if bond_neighbours:
                                vectors = [atom.draw.position for atom in bond_neighbours]
                                gravitational_point = Vector.get_average(vectors)
                                second_line = line.double_line_towards_center(gravitational_point, self.options.bond_spacing, self.options.double_bond_length)
                                second_line_midpoint = second_line.get_midpoint()
                                self._plot_halflines_double(second_line, ax, second_line_midpoint, aromatic)
                            else:
                                print("Shouldn't happen!")
                    else:
                        if self._is_terminal(bond.atom_1) and self._is_terminal(bond.atom_2):
                            dummy_1 = Vector(bond.atom_1.draw.position.x + 1, bond.atom_1.draw.position.y + 1)
                            dummy_2 = Vector(bond.atom_1.draw.position.x - 1, bond.atom_1.draw.position.y - 1)
                            double_bond_line_1 = line.double_line_towards_center(dummy_1,
                                                                                 self.options.bond_spacing / 2.0,
                                                                                 self.options.double_bond_length)
                            double_bond_line_1_midpoint = double_bond_line_1.get_midpoint()
                            double_bond_line_2 = line.double_line_towards_center(dummy_2,
                                                                                 self.options.bond_spacing / 2.0,
                                                                                 self.options.double_bond_length)
                            double_bond_line_2_midpoint = double_bond_line_2.get_midpoint()

                            self._plot_halflines_double(double_bond_line_1, ax, double_bond_line_1_midpoint)
                            self._plot_halflines_double(double_bond_line_2, ax, double_bond_line_2_midpoint, aromatic)

                        else:

                            if self._is_terminal(bond.atom_1):
                                terminal_atom = bond.atom_1
                                branched_atom = bond.atom_2
                            else:
                                terminal_atom = bond.atom_2
                                branched_atom = bond.atom_1

                            if len(branched_atom.drawn_neighbours) >= 3:
                                closest_two = self.get_sorted_distances_from_list(terminal_atom, branched_atom.drawn_neighbours)
                                closest_atom_1 = closest_two[0][1]
                                closest_atom_2 = closest_two[1][1]

                                line = Line(terminal_atom.draw.position, branched_atom.draw.position, terminal_atom, branched_atom)

                                double_bond_line_1, double_bond_line_2 = line.get_perpendicular_lines(self.options.bond_spacing / 2.0)
                                terminal_atom_pos_1 = double_bond_line_1.get_atom_coords(terminal_atom)
                                terminal_atom_pos_2 = double_bond_line_2.get_atom_coords(terminal_atom)

                                closest_atom_to_pos_1 = terminal_atom_pos_1.get_closest_atom(closest_atom_1, closest_atom_2)
                                closest_atom_to_pos_2 = terminal_atom_pos_2.get_closest_atom(closest_atom_1, closest_atom_2)

                                bond_1_line = Line(branched_atom.draw.position, closest_atom_to_pos_1.draw.position, branched_atom, closest_atom_to_pos_1)
                                bond_2_line = Line(branched_atom.draw.position, closest_atom_to_pos_2.draw.position, branched_atom, closest_atom_to_pos_2)

                                double_bond_line_1_midpoint = double_bond_line_1.get_midpoint()
                                double_bond_line_2_midpoint = double_bond_line_2.get_midpoint()

                                intersection_1 = double_bond_line_1.find_intersection(bond_1_line)
                                intersection_2 = double_bond_line_2.find_intersection(bond_2_line)

                                if terminal_atom.draw.position.x > branched_atom.draw.position.x:
                                    # check for parallel lines
                                    if intersection_1 and intersection_1.x < 100000 and intersection_1.y < 100000:
                                        double_bond_line_1.point_1 = intersection_1
                                    if intersection_2 and intersection_2.x < 100000 and intersection_2.y < 100000:
                                        double_bond_line_2.point_1 = intersection_2

                                else:
                                    # check for parallel lines
                                    if intersection_1 and intersection_1.x < 100000 and intersection_1.y < 100000:
                                        double_bond_line_1.point_2 = intersection_1
                                    if intersection_2 and intersection_2.x < 100000 and intersection_2.y < 100000:
                                        double_bond_line_2.point_2 = intersection_2

                                self._plot_halflines(double_bond_line_1, ax, double_bond_line_1_midpoint)
                                self._plot_halflines(double_bond_line_2, ax, double_bond_line_2_midpoint, aromatic)

                            else:
                                self._plot_halflines(line, ax, midpoint)

                                bond_neighbours = bond.atom_1.drawn_neighbours + bond.atom_2.drawn_neighbours
                                if bond_neighbours:
                                    vectors = [atom.draw.position for atom in bond_neighbours]
                                    gravitational_point = Vector.get_average(vectors)
                                    second_line = line.get_parallel_line(gravitational_point,
                                                                         self.options.bond_spacing)
                                    second_line_midpoint = second_line.get_midpoint()
                                    self._plot_halflines(second_line, ax, second_line_midpoint, aromatic)
                                else:
                                    print("Shouldn't happen!")

                elif bond.type == 'triple':
                    self._plot_halflines(line, ax, midpoint)
                    line_1, line_2 = line.get_parallel_lines(self.options.bond_spacing)
                    line_1_midpoint = line_1.get_midpoint()
                    line_2_midpoint = line_2.get_midpoint()
                    self._plot_halflines(line_1, ax, line_1_midpoint)
                    self._plot_halflines(line_2, ax, line_2_midpoint)

        for atom in self.structure.graph:
            if atom.draw.positioned:
                text_h = ''
                text_h_pos = None
                if atom.type != 'C' or atom.draw.draw_explicit:
                    if atom.type == 'C':
                        text = '.'
                    else:
                        text = self.set_r_group_indices_subscript(atom.type)
                else:
                    text = ''

                horizontal_alignment = 'center'

                orientation = self._get_hydrogen_text_orientation(atom)
                if orientation == 'H_above_atom':
                    text_h_pos = Vector(atom.draw.position.x, atom.draw.position.y + 6)
                if orientation == 'H_below_atom':
                    text_h_pos = Vector(atom.draw.position.x, atom.draw.position.y - 6)

                atom_draw_position = Vector(atom.draw.position.x, atom.draw.position.y)
                if text == '.':
                    atom_draw_position.y += 2

                if not atom.charge and (atom.type != 'C' or atom.draw.draw_explicit):

                    if atom.draw.has_hydrogen:
                        hydrogen_count = 0
                        for neighbour in atom.neighbours:
                            if neighbour.type == 'H' and not neighbour.draw.is_drawn:
                                hydrogen_count += 1

                        if hydrogen_count and atom.type != 'C':

                            if hydrogen_count > 1:
                                if orientation == 'H_before_atom':
                                    text = r'$H_{hydrogens}{atom_type}$'.format(hydrogens=hydrogen_count,
                                                                                atom_type=atom.type)
                                    horizontal_alignment = 'right'
                                    atom_draw_position.x += 3
                                elif orientation == 'H_below_atom' or orientation == 'H_above_atom':
                                    text = atom.type
                                    text_h = r'$H_{hydrogens}$'.format(hydrogens=hydrogen_count)

                                else:
                                    text = r'${atom_type}H_{hydrogens}$'.format(hydrogens=hydrogen_count,
                                                                                atom_type=atom.type)
                                    horizontal_alignment = 'left'
                                    atom_draw_position.x -= 3
                            elif hydrogen_count == 1:
                                if orientation == 'H_before_atom':
                                    text = f'H{atom.type}'
                                    horizontal_alignment = 'right'
                                    atom_draw_position.x += 3
                                elif orientation == 'H_below_atom' or orientation == 'H_above_atom':
                                    text = atom.type
                                    text_h = 'H'
                                else:
                                    text = f'{atom.type}H'
                                    horizontal_alignment = 'left'
                                    atom_draw_position.x -= 3

                elif atom.charge:
                    if atom.charge > 0:
                        charge_symbol = '+'
                    else:
                        charge_symbol = '-'

                    hydrogen_count = 0
                    for neighbour in atom.neighbours:
                        if neighbour.type == 'H' and not neighbour.draw.is_drawn:
                            hydrogen_count += 1

                    if not hydrogen_count:

                        if abs(atom.charge) > 1:
                            charge_repr = f"{abs(atom.charge)}{charge_symbol}"
                            text = '$' + atom.type + '^{' + charge_repr + '}$'

                            # text = r'${atom_type}^{charge}{charge_symbol}$'.format(charge=atom.charge,
                            #                                                        atom_type=atom.type,
                            #                                                        charge_symbol=charge_symbol)
                        elif abs(atom.charge) == 1:
                            text = r'${atom_type}^{charge_symbol}$'.format(atom_type=atom.type,
                                                                           charge_symbol=charge_symbol)

                        horizontal_alignment = 'left'
                        atom_draw_position.x -= 3
                    else:
                    # elif atom.type != 'C' or atom.draw.draw_explicit:

                        if hydrogen_count > 1:
                            if orientation == 'H_before_atom':
                                if abs(atom.charge) > 1:
                                    charge_repr = f"{abs(atom.charge)}{charge_symbol}"
                                    text = '$H_' + str(hydrogen_count) + atom.type + '^{' + charge_repr + '}$'
                                    # text = r'$H_{hydrogens}{atom_type}^{charge}{charge_symbol}$'.format(hydrogens=hydrogen_count,
                                    #                                                                     atom_type=atom.type,
                                    #                                                                     charge=abs(atom.charge),
                                    #                                                                     charge_symbol=charge_symbol)
                                elif abs(atom.charge) == 1:
                                    text = r'$H_{hydrogens}{atom_type}^{charge_symbol}$'.format(hydrogens=hydrogen_count,
                                                                                                atom_type=atom.type,
                                                                                                charge_symbol=charge_symbol)

                                horizontal_alignment = 'right'
                                atom_draw_position.x += 3
                            elif orientation == 'H_above_atom' or orientation == 'H_below_atom':
                                text_h = r'$H_{hydrogens}$'.format(hydrogens=hydrogen_count)
                                if abs(atom.charge) > 1:
                                    charge_repr = f"{abs(atom.charge)}{charge_symbol}"
                                    text = '$' + atom.type + '^{' + charge_repr + '}$'
                                    # text = r'${atom_type}^{charge}{charge_symbol}$'.format(atom_type=atom.type,
                                    #                                                        charge=abs(atom.charge),
                                    #                                                        charge_symbol=charge_symbol)
                                elif abs(atom.charge) == 1:
                                    text = r'${atom_type}^{charge_symbol}$'.format(atom_type=atom.type,
                                                                                   charge_symbol=charge_symbol)
                            else:
                                if abs(atom.charge) > 1:
                                    charge_repr = f"{abs(atom.charge)}{charge_symbol}"
                                    text = '$' + atom.type + 'H_' + str(hydrogen_count) + '^{' + charge_repr + '}$'

                                    # text = r'${atom_type}H_{hydrogens}^{charge}{charge_symbol}$'.format(hydrogens=hydrogen_count,
                                    #                                                                     atom_type=atom.type,
                                    #                                                                     charge=abs(atom.charge),
                                    #                                                                     charge_symbol=charge_symbol)
                                elif abs(atom.charge) == 1:
                                    text = r'${atom_type}H_{hydrogens}^{charge_symbol}$'.format(hydrogens=hydrogen_count,
                                                                                                atom_type=atom.type,
                                                                                                charge_symbol=charge_symbol)

                                horizontal_alignment = 'left'
                                atom_draw_position.x -= 3
                        elif hydrogen_count == 1:
                            if orientation == 'H_before_atom':
                                if abs(atom.charge) > 1:
                                    charge_repr = f"{abs(atom.charge)}{charge_symbol}"
                                    text = '$H' + atom.type + '^{' + charge_repr + '}$'

                                elif abs(atom.charge) == 1:
                                    text = r'$H{atom_type}^{charge_symbol}$'.format(atom_type=atom.type,
                                                                                    charge_symbol=charge_symbol)
                                horizontal_alignment = 'right'
                                atom_draw_position.x += 3
                            elif orientation == 'H_above_atom' or orientation == 'H_below_atom':
                                text_h = 'H'
                                if abs(atom.charge) > 1:
                                    charge_repr = f"{abs(atom.charge)}{charge_symbol}"
                                    text = '$' + atom.type + '^{' + charge_repr + '}$'

                                    # text = r'${atom_type}^{charge}{charge_symbol}$'.format(atom_type=atom.type,
                                    #                                                        charge=abs(atom.charge),
                                    #                                                        charge_symbol=charge_symbol)
                                elif abs(atom.charge) == 1:
                                    text = r'${atom_type}^{charge_symbol}$'.format(atom_type=atom.type,
                                                                                   charge_symbol=charge_symbol)

                            else:
                                if abs(atom.charge) > 1:
                                    charge_repr = f"{abs(atom.charge)}{charge_symbol}"
                                    text = '$' + atom.type + 'H^{' + charge_repr + '}$'
                                    # text = r'${atom_type}H^{charge}{charge_symbol}$'.format(atom_type=atom.type,
                                    #                                                         charge=abs(atom.charge),
                                    #                                                         charge_symbol=charge_symbol)

                                elif abs(atom.charge) == 1:
                                    text = r'${atom_type}H^{charge_symbol}$'.format(atom_type=atom.type,
                                                                                    charge_symbol=charge_symbol)
                                horizontal_alignment = 'left'
                                atom_draw_position.x -= 3

                if text:
                    plt.text(atom_draw_position.x, atom_draw_position.y,
                             text,
                             horizontalalignment=horizontal_alignment,
                             verticalalignment='center',
                             color=atom.draw.colour)
                if text_h:
                    plt.text(text_h_pos.x, text_h_pos.y,
                             text_h,
                             horizontalalignment='center',
                             verticalalignment='center',
                             color=atom.draw.colour)

    # TODO: make this work with direct SVG writing
    @staticmethod
    def set_r_group_indices_subscript(atom_text: str) -> str:
        """
        Change numbers to subscript in rest group text

        Parameters
        ----------
        atom_text: str, text representing the rest group

        Returns
        -------
        atom_text: str, text representing the rest group with subscript for numbers
        """
        # Take str and return the same str with subscript digits
        # (pattern is necessary to not get confused with isotopes)
        sub_translation = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
        match = re.search('[RXZ]\d+', atom_text)
        if match:
            matched_pattern = match.group()
            adapted_pattern = matched_pattern.translate(sub_translation)
            atom_text = atom_text.replace(matched_pattern, adapted_pattern)
        return atom_text

    @staticmethod
    def _is_terminal(atom: Atom) -> bool:
        """
        Returns True if the atom only has one drawn neighbour and thus is terminal, False otherwise

        Parameters
        ----------
        atom: Atom instance, must be a drawn atom
        """
        if len(atom.drawn_neighbours) <= 1:
            return True

        return False

    def _process_structure(self) -> None:
        """
        Position all atoms and resolve overlaps
        """
        self._position()
        self.structure.refresh_structure()
        self.restore_ring_information()
        self.restore_ring_information()

        self._fix_chiral_bonds_in_rings()

        self.resolve_primary_overlaps()

        self.total_overlap_score, sorted_overlap_scores, atom_to_scores = self.get_overlap_score()

        for i in range(self.options.overlap_resolution_iterations):
            for bond in self.drawn_bonds:
                if self.can_rotate_around_bond(bond):

                    tree_depth_1: int = self.get_subgraph_size(bond.atom_1, {bond.atom_2})
                    tree_depth_2: int = self.get_subgraph_size(bond.atom_2, {bond.atom_1})

                    # Check neighbouring bonds to ensure neither are chiral; only then the bond is rotatable
                    # at the end of the indicated atom

                    atom_1_rotatable: bool = True
                    atom_2_rotatable: bool = True

                    for neighbouring_bond in bond.atom_1.bonds:
                        if neighbouring_bond.type == 'double' and neighbouring_bond.chiral:
                            atom_1_rotatable = False

                    for neighbouring_bond in bond.atom_2.bonds:
                        if neighbouring_bond.type == 'double' and neighbouring_bond.chiral:
                            atom_2_rotatable = False

                    # If neither are rotatable, continue

                    if not atom_1_rotatable and not atom_2_rotatable:
                        continue
                    elif atom_1_rotatable and not atom_2_rotatable:
                        atom_2 = bond.atom_1
                        atom_1 = bond.atom_2
                    elif atom_2_rotatable and not atom_1_rotatable:
                        atom_1 = bond.atom_1
                        atom_2 = bond.atom_2
                    else:
                        atom_1 = bond.atom_2
                        atom_2 = bond.atom_1

                        if tree_depth_1 > tree_depth_2:
                            atom_1 = bond.atom_1
                            atom_2 = bond.atom_2

                    subtree_overlap_score, _ = self.get_subtree_overlap_score(atom_2, atom_1, atom_to_scores)
                    if subtree_overlap_score > self.options.overlap_sensitivity:
                        neighbours_2 = atom_2.drawn_neighbours[:]
                        neighbours_2.remove(atom_1)

                        if len(neighbours_2) == 1:
                            neighbour = neighbours_2[0]
                            angle = neighbour.draw.position.get_rotation_away_from_vector(atom_1.draw.position,
                                                                                          atom_2.draw.position,
                                                                                          math.radians(120))

                            self.rotate_subtree(neighbour, atom_2, angle, atom_2.draw.position)

                            new_overlap_score, _, _ = self.get_overlap_score()
                            if new_overlap_score > self.total_overlap_score:
                                self.rotate_subtree(neighbour, atom_2, -angle, atom_2.draw.position)
                            else:
                                self.total_overlap_score = new_overlap_score

                        elif len(neighbours_2) == 2:
                            if atom_2.draw.rings and atom_1.draw.rings:
                                continue

                            neighbour_1: Atom = neighbours_2[0]
                            neighbour_2: Atom = neighbours_2[1]

                            if len(neighbour_1.draw.rings) == 1 and len(neighbour_2.draw.rings) == 1:

                                # If the neighbours are in different rings, or in rings at all, do nothing
                                if neighbour_1.draw.rings[0] != neighbour_2.draw.rings[0]:
                                    continue

                            elif neighbour_1.draw.rings or neighbour_2.draw.rings:
                                continue
                            else:
                                angle_1 = neighbour_1.draw.position.get_rotation_away_from_vector(atom_1.draw.position,
                                                                                                  atom_2.draw.position,
                                                                                                  math.radians(120))
                                angle_2 = neighbour_2.draw.position.get_rotation_away_from_vector(atom_1.draw.position,
                                                                                                  atom_2.draw.position,
                                                                                                  math.radians(120))

                                self.rotate_subtree(neighbour_1, atom_2, angle_1, atom_2.draw.position)
                                self.rotate_subtree(neighbour_2, atom_2, angle_2, atom_2.draw.position)

                                new_overlap_score, _, _ = self.get_overlap_score()

                                if new_overlap_score > self.total_overlap_score:
                                    self.rotate_subtree(neighbour_1, atom_2, -angle_1, atom_2.draw.position)
                                    self.rotate_subtree(neighbour_2, atom_2, -angle_2, atom_2.draw.position)
                                else:
                                    self.total_overlap_score = new_overlap_score

                        self.total_overlap_score, sorted_overlap_scores, atom_to_scores = self.get_overlap_score()

        for _ in range(self.options.overlap_resolution_iterations):
            if self.options.finetune:
                self._finetune_overlap_resolution()
                self.total_overlap_score, sorted_overlap_scores, atom_to_scores = self.get_overlap_score()

        for i in range(self.options.overlap_resolution_iterations):

            self.resolve_secondary_overlaps(sorted_overlap_scores)

    def _position(self) -> None:
        """
        Position all atoms
        """
        start_atom: Optional[Atom] = None

        for atom in self.structure.graph:
            if atom.draw.bridged_ring is not None:
                start_atom = atom
                break

        for ring in self.rings:
            if ring.bridged:
                start_atom = ring.members[0]

        if len(self.rings) > 0 and start_atom is None:
            start_atom = self.rings[0].members[0]

        if start_atom is None:
            for atom in self.drawn_atoms:
                if self._is_terminal(atom):
                    start_atom = atom
                    break

        if start_atom is None:
            start_atom = list(self.drawn_atoms)[0]

        # Iteratively position all atoms

        self.create_next_bond(start_atom, None, 0.0)

    # TODO: improve docstring
    def create_next_bond(self, atom: Atom, previous_atom: Optional[Atom] = None, angle: float = 0.0,
                         previous_branch_shortest: bool = False) -> None:
        """
        Iteratively position atoms one bond at a time

        Parameters
        ----------
        atom: Atom instance, current atom
        previous_atom: Atom instance, previous atom
        angle: float, angle of the previous atom with respect to the previous bond
        previous_branch_shortest: bool, True if the previous branch has the shortest tree depth, False otherwise

        """

        # TODO: check if this is necessary

        if atom.draw.positioned:
            return

        if previous_atom is None:
            # Create a 'dummy' previous position if the atom is the first
            # TODO: Check that the dummy is placed and rotated correctly
            dummy: Vector = Vector(self.options.bond_length, 0)
            dummy.rotate(math.radians(-60.0))

            atom.draw.previous_position = dummy
            atom.draw.previous_atom = None
            # Place the first atom one bond length right to the origin
            atom.draw.set_position(Vector(self.options.bond_length, 0))

            atom.draw.angle = math.radians(-60.0)

            # If the atom is not part of a bridged ring, it is now positioned

            if atom.draw.bridged_ring is None:
                atom.draw.positioned = True

        # If the previous atom was part of a ring

        elif len(previous_atom.draw.rings) > 0:
            # note that, as atoms within a ring all get placed all at the same time, this means that the current atom
            # is never in the same ring as the previous atom
            neighbours = previous_atom.drawn_neighbours
            joined_vertex = None
            # Initialise position to the origin
            position: Vector = Vector(0, 0)

            # If the previous atom was not part of a bridged ring and the previous atom was part of more than one ring
            # Example: the central two carbon atoms in a tryptophan ring system

            if previous_atom.draw.bridged_ring is None and len(previous_atom.draw.rings) > 1:
                # Find the vertex adjoining the current bridged ring that is also in both ring systems. This is the
                # joined vertex.
                for neighbour in neighbours:
                    if len(set(neighbour.draw.rings) & set(previous_atom.draw.rings)) == len(previous_atom.draw.rings):
                        joined_vertex = neighbour
                        break

            # If there is no joined vertex because the atom is part of a bridged ring OR it is only part of one ring

            if not joined_vertex:
                # TODO: Check that this also works for bridged rings
                # Position the atom one bond length away from the previous atom, with the bond perpendicular
                # to the imaginary line between both ring neighbours

                # First find the ring neighbours:

                for neighbour in neighbours:

                    if neighbour.draw.positioned and self.atoms_are_in_same_ring(neighbour, previous_atom):
                        # Calculate a vector that is perpendicular to the imaginary line between both neighbours
                        # and points towards the center of gravity of the neighbours
                        position.add(Vector.subtract_vectors(neighbour.draw.position, previous_atom.draw.position))

                # Invert this vector around the origin to obtain a vector in the desired direction
                position.invert()

                # Scale the vector such that is has the length of a bond
                position.normalise()
                position.multiply_by_scalar(self.options.bond_length)

                # Add the position of the previous atom to the vector to get the position of the current atom
                position.add(previous_atom.draw.position)

            else:
                # If there IS a joined vertex, place the new atom exactly 180 degrees away from this vertex
                # with the previous atom as rotating point, placing the new atom exactly opposite the vertex

                position = joined_vertex.draw.position.copy()
                position.rotate_around_vector(math.pi, previous_atom.draw.position)

            atom.draw.set_previous_position(previous_atom)
            atom.draw.set_position(position)
            atom.draw.positioned = True

        else:
            # If the previous atom was not part of a ring, simply create a vector of length bond length,
            # rotate it around the origin, and translate it to the position of the previous atom to obtain
            # the desired position
            position: Vector = Vector(self.options.bond_length, 0)
            position.rotate(angle)
            position.add(previous_atom.draw.position)

            atom.draw.set_position(position)
            atom.draw.set_previous_position(previous_atom)
            atom.draw.positioned = True

        # After we have placed the atom, we have to:
        #   1. Also place the atoms from any ring systems that the atom is a part of.
        #   2. Find the desired angle between this atom and the next atom based on adjacent branches

        # If the atom is part of a ring: position the entire bridged ring at once
        if len(atom.draw.rings) > 0:
            # Prioritise the placement of bridged rings
            if atom.draw.bridged_ring:
                next_ring: Ring = self.id_to_ring[atom.draw.bridged_ring]
            else:
                next_ring: Ring = self.id_to_ring[atom.draw.rings[0]]

            # Only position the ring if it has not yet been positioned

            if not next_ring.positioned:
                # Get a vector in the direction of the previous atom
                next_center = Vector.subtract_vectors(atom.draw.previous_position, atom.draw.position)
                # Invert that vector around the origin to get a vector into the direction of the ring center
                next_center.invert()

                # Scale the vector to the radius of the ring
                next_center.normalise()
                radius: float = Polygon.find_polygon_radius(self.options.bond_length, len(next_ring.members))
                next_center.multiply_by_scalar(radius)

                # Add the position of the previous atom to translate the vector into the correct position
                next_center.add(atom.draw.position)

                # Finally, create the ring, specifying the current atom as an anchor
                # Note that this function will automatically call create_next_bond once the ring has been placed
                self.create_ring(next_ring, next_center, atom)

        # If the atom is not part of a ring, position just the atom
        else:
            neighbours: List[Atom] = atom.drawn_neighbours[:]

            if previous_atom:
                if previous_atom in neighbours:
                    neighbours.remove(previous_atom)

            previous_angle: float = atom.draw.get_angle()

            if len(neighbours) == 1:

                # Note: this atom cannot have been positioned before; if it were, it would be part of a ring
                # TODO: Should we add an assert statement here to confirm the atom has not been positioned?
                next_atom: Atom = neighbours[0]

                current_bond: Bond = self.structure.bond_lookup[atom][next_atom]
                previous_bond: Optional[Bond] = None

                if previous_atom:
                    previous_bond = self.structure.bond_lookup[previous_atom][atom]

                # There are certain bond configurations which lead to the molecule being linear around those bonds

                if current_bond.type == 'triple' or (previous_bond and previous_bond.type == 'triple') or \
                        (current_bond.type == 'double' and previous_bond and previous_bond.type == 'double' and
                         previous_atom and len(previous_atom.draw.rings) == 0 and
                         len(atom.neighbours) == 2):

                    # If the atoms connected by two subsequent double bonds are carbons, it is difficult to recognise
                    # these bonds as two separate bonds. Therefore, we must make sure that such carbons are explicitly
                    # drawn to break up the continuous line in the drawing

                    if current_bond.type == 'double' and previous_bond.type == 'double':

                        atom.draw.draw_explicit = True

                    if current_bond.type == 'triple':
                        atom.draw.draw_explicit = True
                        next_atom.draw.draw_explicit = True

                    # TODO: Use this attribute to make sure linear double bonds align

                    if previous_atom:
                        previous_bond.draw.center = True

                    current_bond.draw.center = True

                    # TODO: Is there ever a situation where this is not true?

                    if current_bond.type == 'double' or current_bond.type == 'triple' or \
                            (previous_atom and previous_bond.type == 'triple'):
                        next_atom.draw.angle = 0.0

                    self.create_next_bond(next_atom, atom, previous_angle + next_atom.draw.angle)

                # If the previous atom was in a ring

                elif previous_atom and len(previous_atom.draw.rings) > 0:

                    # Two possible positions for a single neighbour:
                    # at an angle of 60 or -60 degrees
                    proposed_angle_1: float = math.radians(60.0)
                    proposed_angle_2: float = proposed_angle_1 * -1

                    # Create vectors around the origin that have the length of the bond

                    proposed_vector_1 = Vector(self.options.bond_length, 0)
                    proposed_vector_2 = Vector(self.options.bond_length, 0)

                    # Rotate them into their desired directions.

                    proposed_vector_1.rotate(proposed_angle_1 + atom.draw.get_angle())
                    proposed_vector_2.rotate(proposed_angle_2 + atom.draw.get_angle())

                    # Move the vector relative to the previous atom

                    proposed_vector_1.add(atom.draw.position)
                    proposed_vector_2.add(atom.draw.position)

                    # Choose the position that is the furthest away from the centre of mass of the structure

                    centre_of_mass: Vector = self.get_current_centre_of_mass()
                    distance_1: float = proposed_vector_1.get_squared_distance(centre_of_mass)
                    distance_2: float = proposed_vector_2.get_squared_distance(centre_of_mass)

                    if distance_1 < distance_2:
                        next_atom.draw.angle = proposed_angle_2
                    else:
                        next_atom.draw.angle = proposed_angle_1

                    self.create_next_bond(next_atom, atom, previous_angle + next_atom.draw.angle)

                # If the previous atom was not in a ring

                else:

                    # Set the new angle to the old angle to start with, aligning the two atoms

                    proposed_angle: float = atom.draw.angle

                    # If the previous atom had more than 3 neighbours

                    if previous_atom and len(previous_atom.drawn_neighbours) > 3:
                        if round(proposed_angle, 2) > 0.00:
                            proposed_angle = min([math.radians(60), proposed_angle])
                        elif round(proposed_angle, 2) < 0.00:
                            proposed_angle = max([-math.radians(60), proposed_angle])
                        else:
                            proposed_angle = math.radians(60)

                    elif proposed_angle is None:
                        last_angled_atom = self.get_last_atom_with_angle(atom)
                        proposed_angle = last_angled_atom.draw.angle

                        if proposed_angle is None:
                            proposed_angle = math.radians(60)

                    rotatable = True

                    if previous_atom:

                        bond = self.structure.bond_lookup[previous_atom][atom]
                        if bond.type == 'double' and bond.chiral:
                            rotatable = False

                            previous_previous_atom = previous_atom.draw.previous_atom

                            if previous_previous_atom:

                                configuration = bond.chiral_dict[previous_previous_atom][next_atom]
                                if configuration == 'cis':

                                    proposed_angle = -proposed_angle

                    if rotatable:

                        if previous_branch_shortest:
                            next_atom.draw.angle = proposed_angle
                        else:
                            next_atom.draw.angle = -proposed_angle
                    else:
                        next_atom.draw.angle = -proposed_angle
                    #
                    # if round(math.degrees(next_atom.draw.angle), 0) == 360 or \
                    #         round(math.degrees(next_atom.draw.angle), 0) == -360 or \
                    #         round(math.degrees(next_atom.draw.angle), 0) == 0:
                    #     print("Because we have an angle of 0", atom, next_atom)
                    #     atom.draw.draw_explicit = True

                    self.create_next_bond(next_atom, atom, previous_angle + next_atom.draw.angle)

            elif len(neighbours) == 2:
                proposed_angle = atom.draw.angle
                if not proposed_angle:
                    proposed_angle = math.radians(60)

                neighbour_1, neighbour_2 = neighbours

                subgraph_1_size = self.get_subgraph_size(neighbour_1, {atom})
                subgraph_2_size = self.get_subgraph_size(neighbour_2, {atom})

                if previous_atom:
                    subgraph_3_size = self.get_subgraph_size(previous_atom, {atom})

                else:
                    subgraph_3_size = 0

                cis_atom_index = 0
                trans_atom_index = 1

                if neighbour_2.type == 'C' and neighbour_1.type != 'C' and subgraph_2_size > 1 and subgraph_1_size < 5:
                    cis_atom_index = 1
                    trans_atom_index = 0

                elif neighbour_2.type != 'C' and neighbour_1.type == 'C' and subgraph_1_size > 1 and subgraph_2_size < 5:
                    cis_atom_index = 0
                    trans_atom_index = 1

                elif subgraph_2_size > subgraph_1_size:
                    cis_atom_index = 1
                    trans_atom_index = 0

                cis_atom = neighbours[cis_atom_index]
                trans_atom = neighbours[trans_atom_index]

                previous_branch_shortest = False

                if subgraph_3_size < subgraph_2_size and subgraph_3_size < subgraph_1_size:
                    previous_branch_shortest = True

                trans_atom.draw.angle = proposed_angle
                cis_atom.draw.angle = -proposed_angle

                cis_bond = self.structure.bond_lookup[atom][cis_atom]
                trans_bond = self.structure.bond_lookup[atom][trans_atom]

                # if the cis bond and trans bond are single bonds it can be adjacent to a
                # chiral double bond and may have to be placed in a specific orientation
                if cis_bond.type == 'single' and trans_bond.type == 'single':

                    if previous_atom:
                        previous_bond = self.structure.bond_lookup[atom][previous_atom]

                        # checks if the previous bond was a chiral double bond
                        # TODO: make sure chiral bonds aren't drawn first!
                        if previous_bond.type == 'double' and previous_bond.chiral:
                            if previous_atom.draw.previous_atom:
                                configuration_cis_atom = previous_bond.chiral_dict[previous_atom.draw.previous_atom][cis_atom]
                                if configuration_cis_atom == 'cis':
                                    trans_atom.draw.angle = -proposed_angle
                                    cis_atom.draw.angle = proposed_angle

                self.create_next_bond(trans_atom, atom, previous_angle + trans_atom.draw.angle, previous_branch_shortest)
                self.create_next_bond(cis_atom, atom, previous_angle + cis_atom.draw.angle, previous_branch_shortest)

            elif len(neighbours) == 3:
                # This means the atom has four outgoing bonds. We can draw this in two different ways:
                #   1. As a regular cross, with equal angles between all atoms
                #   2. As a butterfly, with one 
                subgraph_1_size = self.get_subgraph_size(neighbours[0], {atom})
                subgraph_2_size = self.get_subgraph_size(neighbours[1], {atom})
                subgraph_3_size = self.get_subgraph_size(neighbours[2], {atom})

                straight_atom = neighbours[0]
                left_atom = neighbours[1]
                right_atom = neighbours[2]

                if subgraph_2_size > subgraph_1_size and subgraph_2_size > subgraph_3_size:
                    straight_atom = neighbours[1]
                    left_atom = neighbours[0]
                    right_atom = neighbours[2]

                elif subgraph_3_size > subgraph_1_size and subgraph_3_size > subgraph_2_size:
                    straight_atom = neighbours[2]
                    left_atom = neighbours[0]
                    right_atom = neighbours[1]

                if previous_atom and len(previous_atom.draw.rings) < 1\
                        and len(straight_atom.draw.rings) < 1\
                        and len(left_atom.draw.rings) < 1\
                        and len(right_atom.draw.rings) < 1\
                        and self.get_subgraph_size(left_atom, {atom}) == 1\
                        and self.get_subgraph_size(right_atom, {atom}) == 1\
                        and self.get_subgraph_size(straight_atom, {atom}) > 1:
                    straight_atom.draw.angle = atom.draw.angle * -1
                    if atom.draw.angle >= 0:
                        left_atom.draw.angle = math.radians(30)
                        right_atom.draw.angle = math.radians(90)
                    else:
                        left_atom.draw.angle = math.radians(-30)
                        right_atom.draw.angle = math.radians(-90)

                else:
                    straight_atom.draw.angle = 0.0
                    left_atom.draw.angle = math.radians(90)
                    right_atom.draw.angle = math.radians(-90)

                self.create_next_bond(straight_atom, atom, previous_angle + straight_atom.draw.angle)
                self.create_next_bond(left_atom, atom, previous_angle + left_atom.draw.angle)
                self.create_next_bond(right_atom, atom, previous_angle + right_atom.draw.angle)

            elif len(neighbours) == 4:
                subgraph_1_size = self.get_subgraph_size(neighbours[0], {atom})
                subgraph_2_size = self.get_subgraph_size(neighbours[1], {atom})
                subgraph_3_size = self.get_subgraph_size(neighbours[2], {atom})
                subgraph_4_size = self.get_subgraph_size(neighbours[3], {atom})

                atom_1 = neighbours[0]
                atom_2 = neighbours[1]
                atom_3 = neighbours[2]
                atom_4 = neighbours[3]

                if subgraph_2_size > subgraph_1_size and subgraph_2_size > subgraph_3_size\
                        and subgraph_2_size > subgraph_4_size:
                    atom_1 = neighbours[1]
                    atom_2 = neighbours[0]

                elif subgraph_3_size > subgraph_1_size and subgraph_3_size > subgraph_2_size\
                        and subgraph_3_size > subgraph_4_size:
                    atom_1 = neighbours[2]
                    atom_2 = neighbours[0]
                    atom_3 = neighbours[1]

                elif subgraph_4_size > subgraph_1_size and subgraph_4_size > subgraph_2_size\
                        and subgraph_4_size > subgraph_3_size:
                    atom_1 = neighbours[3]
                    atom_2 = neighbours[0]
                    atom_3 = neighbours[1]
                    atom_4 = neighbours[2]

                atom_1.draw.angle = math.radians(-36)
                atom_2.draw.angle = math.radians(36)
                atom_3.draw.angle = math.radians(-108)
                atom_4.draw.angle = math.radians(108)

                self.create_next_bond(atom_1, atom, previous_angle + atom_1.draw.angle)
                self.create_next_bond(atom_2, atom, previous_angle + atom_2.draw.angle)
                self.create_next_bond(atom_3, atom, previous_angle + atom_3.draw.angle)
                self.create_next_bond(atom_4, atom, previous_angle + atom_4.draw.angle)
                
    def show_atom_numbers(self, atom_list: List["Atom"]) -> None:
        """
        Superpose atom numbers onto the drawing

        Parameters
        ----------
        atom_list: list of [atom, ->], with each atom an atom instance. List of atoms to visualise atom numbers for

        """
        for atom in atom_list:
            if atom in self.structure.graph:
                atom_drawing = self.structure.get_atom(atom)
                atom_string = str(atom)
                text_x = atom_drawing.draw.position.x + 3
                text_y = atom_drawing.draw.position.y + 8
                text = self._draw_text(atom_string, text_x, text_y)
                self._add_svg_element(text, atom)

    def restore_ring_information(self) -> None:
        """
        Restore the original ring information prior to fusing of bridged abd joined rings
        """
        bridged_rings: List[Ring] = self.get_bridged_rings()

        self.rings: List[Ring] = []
        self.ring_overlaps: List[RingOverlap] = []

        for ring in bridged_rings:
            for subring in ring.subrings:
                self.original_rings[subring.id].center = subring.center

        for ring in self.original_rings:
            for i, atom in enumerate(ring.members):
                positioned_atom = self.atom_nr_to_atom[atom.nr]
                ring.members[i] = positioned_atom
            self.set_ring_center(ring)
            self.rings.append(ring)

        for ring_overlap in self.original_ring_overlaps:
            self.ring_overlaps.append(ring_overlap)

        for atom in self.structure.graph:
            atom.draw.restore_rings()

    @staticmethod
    def bond_is_rotatable(bond: Bond) -> bool:
        """
        Returns True if the bond is rotatable in the drawing, False otherwise

        Parameters
        ----------
        bond: Bond instance

        Returns
        -------
        True if the bond is rotatable, i.e. it is not fixed in place due to stereochemical restraints; False otherwise
        """
        # If a bond is part of a ring, we can't rotate it
        if bond.atom_1.draw.rings and \
                bond.atom_2.draw.rings and \
                len(set(bond.atom_1.draw.rings).intersection(set(bond.atom_2.draw.rings))) > 0:
            return False

        if bond.type != 'single':
            # If a non-single bond is stereochemically restrained, it is not rotatable
            if bond.chiral:
                return False

            if len(bond.atom_1.drawn_neighbours) > 1 and len(bond.atom_2.drawn_neighbours) > 1:
                return False

        chiral = False

        # Also bonds adjacent to stereochemically restrained bonds are not rotatable
        for bond_1 in bond.atom_1.bonds:
            if bond_1.chiral:
                chiral = True
                break

        for bond_2 in bond.atom_2.bonds:
            if bond_2.chiral:
                chiral = True
                break

        if chiral:
            return False
        
        if bond.chiral_symbol:
            return False
        
        return True

    # TODO: Either remove this function or the one above
    @staticmethod
    def can_rotate_around_bond(bond: Bond) -> bool:
        """
        More lenient version of bond_is_rotatable. Return True if the bond is rotatable, False otherwise

        Parameters
        ----------
        bond: Bond instance

        Returns
        -------
        True if the bond is rotatable, i.e. it is not fixed in place due to stereochemical restraints; False otherwise

        """

        if bond.type != 'single':
            return False

        # If bond is terminal, don't bother rotating.

        if len(bond.atom_1.drawn_neighbours) == 1 or len(bond.atom_2.drawn_neighbours) == 1:
            return False

        # Added this, needs extensive checking

        if bond.atom_1.draw.rings and \
                bond.atom_2.draw.rings and \
                len(set(bond.atom_1.draw.rings).intersection(set(bond.atom_2.draw.rings))) > 0:
            return False

        return True

    def resolve_primary_overlaps(self) -> None:
        """
        Resolve overlaps when a ring has two outgoing edges
        """
        overlaps: list[dict[str, Union[Atom, list[int], list[Atom]]]] = []

        # Keep track of which atoms are resolved
        resolved_atoms: dict[Atom, bool] = {}

        for atom in self.structure.graph:
            if atom.draw.is_drawn:
                resolved_atoms[atom] = False

        for ring in self.rings:
            for atom in ring.members:
                if resolved_atoms[atom]:
                    continue

                resolved_atoms[atom] = True

                if not atom._adjacent_to_stereobond():

                    non_ring_neighbours = self.get_non_ring_neighbours(atom)

                    # If the ring has more than one outgoing edge of this ring, or if it has one outgoing edge but is
                    # part of two rings, an overlap needs to be resolved

                    if len(non_ring_neighbours) > 1 or (len(non_ring_neighbours) == 1 and len(atom.draw.rings) == 2):
                        overlaps.append({'common': atom,
                                         'rings': atom.draw.rings,
                                         'vertices': non_ring_neighbours})

        for overlap in overlaps:
            branches_to_adjust = overlap['vertices']
            rings = overlap['rings']
            root = overlap['common']

            # Split the two outgoing branches
            if len(branches_to_adjust) == 2:

                atom_1, atom_2 = branches_to_adjust

                if not atom_1.draw.is_drawn or not atom_2.draw.is_drawn:
                    continue

                angle = (2 * math.pi - self.id_to_ring[rings[0]].get_angle()) / 6.0

                self.rotate_subtree(atom_1, root, angle, root.draw.position)
                self.rotate_subtree(atom_2, root, -angle, root.draw.position)

                total, sorted_scores, atom_to_score = self.get_overlap_score()
                subtree_overlap_atom_1_1, _ = self.get_subtree_overlap_score(atom_1, root, atom_to_score)
                subtree_overlap_atom_2_1, _ = self.get_subtree_overlap_score(atom_2, root, atom_to_score)
                total_score = subtree_overlap_atom_1_1 + subtree_overlap_atom_2_1

                self.rotate_subtree(atom_1, root, -2.0 * angle, root.draw.position)
                self.rotate_subtree(atom_2, root, 2.0 * angle, root.draw.position)

                total, sorted_scores, atom_to_score = self.get_overlap_score()
                subtree_overlap_atom_1_2, _ = self.get_subtree_overlap_score(atom_1, root, atom_to_score)
                subtree_overlap_atom_2_2, _ = self.get_subtree_overlap_score(atom_2, root, atom_to_score)
                total_score_2 = subtree_overlap_atom_1_2 + subtree_overlap_atom_2_2

                if total_score_2 > total_score:
                    self.rotate_subtree(atom_1, root, 2.0 * angle, root.draw.position)
                    self.rotate_subtree(atom_2, root, -2.0 * angle, root.draw.position)

            elif len(branches_to_adjust) == 1:
                # This situation is already resolved in the initial positioning
                if len(rings) == 2:
                    pass

    def resolve_secondary_overlaps(self, sorted_scores: List[Tuple[float, Atom]]) -> None:
        """
        Resolve overlaps between non-adjacent atoms that are less than a bond distance away by rotating the bonds of
        terminal atoms

        Parameters
        ----------
        sorted_scores: reverse sorted list of [(score, atom), ->], with score a float (the lower the better) and
            atom an Atom instance. The score indicates how many atoms the atom overlaps with

        Returns
        -------

        """
        for score, atom in sorted_scores:
            if score > self.options.overlap_sensitivity:
                if len(atom.drawn_neighbours) <= 1:
                    if atom.drawn_neighbours and atom.drawn_neighbours[0]._adjacent_to_stereobond():
                        continue

                    closest_atom = self.get_closest_atom(atom)

                    drawn_neighbours = closest_atom.drawn_neighbours

                    if len(drawn_neighbours) <= 1:
                        if not closest_atom.draw.previous_position:
                            closest_position = drawn_neighbours[0].draw.position
                        else:
                            closest_position = closest_atom.draw.previous_position

                    else:
                        if not closest_atom.draw.previous_position:
                            closest_position = drawn_neighbours[0].draw.position
                        else:
                            closest_position = closest_atom.draw.position

                    if not atom.draw.previous_position:
                        atom_previous_position = atom.drawn_neighbours[0].draw.position
                    else:
                        atom_previous_position = atom.draw.previous_position

                    atom.draw.position.rotate_away_from_vector(closest_position, atom_previous_position,
                                                               math.radians(20))

    # TODO: Probably not a necessary function - see if it is removable
    def get_atom_nr_to_atom(self) -> None:
        """
        Returns a dictionary of {atom_index: atom, ->}, representing all atoms in the structure
        """

        self.atom_nr_to_atom = {}
        for atom in self.structure.graph:
            self.atom_nr_to_atom[atom.nr] = atom

    def get_subtree_overlap_score(self, root: Atom, root_parent: Atom,
                                  atom_to_score: Dict[Atom, float]) -> Tuple[float, Vector]:
        """
        Returns the overlap score of a subtree of the drawn structure, as well as the centre of that subtree

        Parameters
        ----------
        root: Atom instance, root of the subtree in the structure. Cannot be in a ring.
        root_parent: Atom instance, first atom not to be included in the subtree. Cannot be in a ring.
        atom_to_score: dict of {atom: score, ->}, with atom an Atom instance and score a float (the lower the better)
            representing the number of other atoms the atom overlaps with

        Returns
        -------
        score / count: float, a score normalised by the number of overlapping atoms in the subtree
        center: Vector, midpoint of the subtree

        """
        score = 0.0
        center = Vector(0, 0)

        count = 0

        for atom in self.traverse_substructure(root, {root_parent}):

            subscore = atom_to_score[atom]
            if subscore > self.options.overlap_sensitivity:
                score += subscore
                count += 1

            position = atom.draw.position.copy()
            position.multiply_by_scalar(subscore)
            center.add(position)

        if score:
            center.divide(score)

        if count == 0:
            count = 1

        return score / count, center

    def get_overlap_score(self) -> Tuple[float, List[Tuple[float, Atom]], Dict[Atom, float]]:
        """
        Returns the total overlap score for the structure, the reverse sorted overlap scores per atom, and a dictionary
        of atom to atom-specific overlap scores

        Returns
        -------
        total: float, score representing total overlaps in the structure. The closer the atoms are, the greater the
            penalty
        sorted_overlaps: reverse sorted list of [(score, atom), ->], with score a float (the lower the better) and
            atom an Atom instance. The score indicates how many atoms the atom overlaps with
        overlap_scores: dict of {atom: score, ->}, with atom an Atom instance and score a float (the lower the better)
            representing the number of other atoms the atom overlaps with
        """
        total = 0.0

        overlap_scores = {}
        for atom in self.drawn_atoms:
            overlap_scores[atom] = 0.0

        for i, atom_1 in enumerate(self.drawn_atoms):
            for j in range(i + 1, len(self.drawn_atoms)):
                atom_2 = self.drawn_atoms[j]
                distance = Vector.subtract_vectors(atom_1.draw.position, atom_2.draw.position).get_squared_length()
                if distance < self.options.bond_length_squared:
                    weight = (self.options.bond_length - math.sqrt(distance)) / self.options.bond_length
                    total += weight
                    overlap_scores[atom_1] += weight
                    overlap_scores[atom_2] += weight

        sorted_overlaps: List[Tuple[float, Atom]] = []

        for atom in self.drawn_atoms:
            sorted_overlaps.append((overlap_scores[atom], atom))

        sorted_overlaps.sort(key=lambda x: x[0], reverse=True)

        return total, sorted_overlaps, overlap_scores

    @staticmethod
    def get_non_ring_neighbours(atom: Atom) -> List[Atom]:
        """
        Returns a list of atom neighbours that are not in the same ring as the given atom

        Parameters
        ----------
        atom: Atom instance, atom to determine the non-ring neighbours for

        Returns
        -------
        non_ring_neighbours: list of [atom, ->], with each atom an Atom instance, an atom that is not in the same ring
            as atom

        """
        non_ring_neighbours: list[Atom] = []

        for neighbour in atom.drawn_neighbours:
            nr_overlapping_rings = len(set(atom.draw.rings).intersection(set(neighbour.draw.rings)))
            if nr_overlapping_rings == 0 and not neighbour.draw.is_bridge:
                non_ring_neighbours.append(neighbour)

        return non_ring_neighbours

    def rotate_subtree(self, root: Atom, root_parent: Atom, angle: float, center: Vector):
        """
        Rotate part of the structure around a vector, where root_parent is NOT rotated and root IS rotated

        Parameters
        ----------
        root: Atom instance, first atom to be rotated
        root_parent: Atom instance, atom adjacent to root which is NOT rotated. Note that root and root_parent should
            not share a ring
        angle: float, angle in radians by which to rotate the structure
        center: Vector, point in a plane around which the subtree is rotated

        """

        for atom in self.traverse_substructure(root, {root_parent}):
            atom.draw.position.rotate_around_vector(angle, center)
            for anchored_ring in atom.draw.anchored_rings:
                if anchored_ring.center:
                    anchored_ring.center.rotate_around_vector(angle, center)

    def get_average_position(self) -> Vector:
        """
        Returns the average position (float) of the structure
        """
        sum_x: float = 0.0
        sum_y: float = 0.0

        for atom in self.drawn_atoms:
            sum_x += atom.draw.position.x
            sum_y += atom.draw.position.y

        nr_atoms: int = len(self.drawn_atoms)

        return Vector(sum_x / nr_atoms, sum_y / nr_atoms)

    def rotate_structure(self, angle, midpoint: Optional[Vector] = None) -> None:
        """
        Rotate the entire drawing around a midpoint if given, around its average position otherwise

        Parameters
        ----------
        angle: float, angle in radians by which to rotate the structure
        midpoint: Vector or None. If given, rotate around Vector. If None, calculate the midpoint of the structure as
            drawn and rotate around that
        """

        if midpoint is None:
            midpoint: Vector = self.get_average_position()

        for atom in self.drawn_atoms:
            atom.draw.position.rotate_around_vector(angle, midpoint)
            for anchored_ring in atom.draw.anchored_rings:
                if anchored_ring.center:
                    anchored_ring.center.rotate_around_vector(angle, midpoint)

    def traverse_substructure(self, atom: Atom, visited: Set[Atom]) -> Generator[Atom, None, None]:
        """
        Yield atoms of a subtree, starting at atom and ignoring all bonds leading to atoms in visited

        Parameters
        ----------
        atom: Atom instance, first atom to be yielded
        visited: set of Atom instances, atoms that have already been yielded (or are ignored)

        Returns
        -------
        atom: Atom instance
        """
        yield atom
        visited.add(atom)
        for neighbour in atom.drawn_neighbours:
            if neighbour not in visited:
                yield from self.traverse_substructure(neighbour, visited)

    def get_subgraph_size(self, atom: Atom, masked_atoms: Set[Atom]) -> int:
        """
        Returns the size of a subtree starting at atom and ignoring all bonds adjacent to atoms in masked_atoms

        Parameters
        ----------

        atom: Atom instance, starting point of subtree
        masked_atoms: set of {atom, ->}, atoms to be ignored

        Returns
        -------
        int, number of atoms in the subtree

        """
        masked_atoms.add(atom)

        for neighbour in atom.drawn_neighbours:
            if neighbour not in masked_atoms:
                self.get_subgraph_size(neighbour, masked_atoms)

        return len(masked_atoms) - 1

    @staticmethod
    def get_last_atom_with_angle(atom: Atom) -> Optional[Atom]:
        """
        Returns the angle of the last visited atom

        Parameters
        ----------
        atom: Atom instance, atom for which to determine the angle

        Returns
        -------
        parent_atom: Atom instance, last atom to be assigned an angle

        """
        parent_atom: Optional[Atom] = atom.draw.previous_atom
        angle: float = parent_atom.draw.angle

        while parent_atom and not angle:
            parent_atom = parent_atom.draw.previous_atom
            angle = parent_atom.draw.angle

        return parent_atom

    def create_ring(self, ring: Ring, center: Optional[Vector] = None, start_atom: Optional[Atom] = None,
                    previous_atom: Optional[Atom] = None) -> None:
        """
        Position all atoms in a ring

        Parameters
        ----------
        ring: Ring instance
        center: Vector, center of the ring. if not given, initialise at 0,0
        start_atom: Atom instance, first atom of the ring
        previous_atom: Atom instance, last atom to be positioned (not in a ring)

        """
        # if the ring was already positioned, return
        if ring.positioned:
            return

        if center is None:
            center = Vector(0, 0)

        ordered_neighbour_ids: list[int] = ring.get_ordered_neighbours(self.ring_overlaps)
        starting_angle: float = 0

        if start_atom:
            starting_angle = Vector.subtract_vectors(start_atom.draw.position, center).angle()

        ring_size: int = len(ring.members)

        radius: float = Polygon.find_polygon_radius(self.options.bond_length, ring_size)
        angle: float = Polygon.get_central_angle(ring_size)

        ring.central_angle = angle

        # if the starting atom is not in the ring members, reassign the first atom from the ring to this variable
        if start_atom not in ring.members:
            if start_atom:
                start_atom.draw.positioned = False
            start_atom = ring.members[0]

        # if the ring is bridged, i.e. contains multiple composite rings, position it with the KK-Layout algorithm

        if ring.bridged:
            KKLayout(self.structure, ring.members, center, start_atom,
                     self.options.bond_length, self.options.kk_threshold, self.options.kk_inner_threshold,
                     self.options.kk_max_iteration, self.options.kk_max_inner_iteration,
                     self.options.kk_max_energy)
            ring.positioned = True

            self.set_ring_center(ring)
            center = ring.center

            for subring in ring.subrings:
                self.set_ring_center(subring)

        # Otherwise, use simple polygon mathematics to position the ring
        else:
            ring.set_member_positions(self.structure, start_atom, previous_atom, center, starting_angle, radius, angle)

        ring.positioned = True
        ring.center = center

        # Also position any rings directly adjacent to the current ring

        for neighbour_id in ordered_neighbour_ids:
            neighbour = self.id_to_ring[neighbour_id]
            if neighbour.positioned:
                continue

            atoms = list(RingOverlap.get_vertices(self.ring_overlaps, ring.id, neighbour.id))

            if len(atoms) == 2:

                # This ring is fused, i.e. it shares a bond with another ring
                ring.fused = True
                neighbour.fused = True

                atom_1 = atoms[0]
                atom_2 = atoms[1]

                midpoint = Vector.get_midpoint(atom_1.draw.position, atom_2.draw.position)
                normals = Vector.get_normals(atom_1.draw.position, atom_2.draw.position)

                normals[0].normalise()
                normals[1].normalise()

                apothem = Polygon.get_apothem_from_side_length(self.options.bond_length,
                                                               len(neighbour.members))

                normals[0].multiply_by_scalar(apothem)
                normals[1].multiply_by_scalar(apothem)

                normals[0].add(midpoint)
                normals[1].add(midpoint)

                next_center = normals[0]

                distance_to_center_1 = Vector.subtract_vectors(center, normals[0]).get_squared_length()
                distance_to_center_2 = Vector.subtract_vectors(center, normals[1]).get_squared_length()

                # Draw the ring on the side of the bond that lies opposite the current ring!
                if distance_to_center_2 > distance_to_center_1:
                    next_center = normals[1]

                position_1 = Vector.subtract_vectors(atom_1.draw.position, next_center)
                position_2 = Vector.subtract_vectors(atom_2.draw.position, next_center)

                # We have to make sure to pass the two atoms of the shared bond in the right order, such that the next
                # ring s positioned such that it is properly aligned with the current one
                if position_1.get_clockwise_orientation(position_2) == 'clockwise':
                    if not neighbour.positioned:
                        self.create_ring(neighbour, next_center, atom_1, atom_2)

                else:
                    if not neighbour.positioned:
                        self.create_ring(neighbour, next_center, atom_2, atom_1)

            elif len(atoms) == 1:
                # This ring is a spiro, i.e. it shares a single atom with another ring
                ring.spiro = True
                neighbour.spiro = True

                atom = atoms[0]
                next_center = Vector.subtract_vectors(center, atom.draw.position)

                # In the case of a spiro, we can simply mirror the previous ring center in the shared atom,
                # normalise the distance, and then multiply it by the required polygon radius

                next_center.invert()
                next_center.normalise()

                distance_to_center = Polygon.find_polygon_radius(self.options.bond_length, len(neighbour.members))
                next_center.multiply_by_scalar(distance_to_center)
                next_center.add(atom.draw.position)

                if not neighbour.positioned:
                    self.create_ring(neighbour, next_center, atom)

        for atom in ring.members:
            for neighbour in atom.drawn_neighbours:
                if neighbour.draw.positioned:
                    continue

                atom.draw.connected_to_ring = True
                # If there are bonds coming out of the ring, position it
                self.create_next_bond(neighbour, atom, 0.0)

    @staticmethod
    def set_ring_center(ring: Ring) -> None:
        """
        Calculate and set the center of the ring

        Parameters
        ----------
        ring: Ring instance

        """
        total: Vector = Vector(0, 0)
        for atom in ring.members:
            total.add(atom.draw.position)

        total.divide(len(ring.members))

        ring.center = total

    # TODO: Check if we can use this function instead of get_average_midpoint
    def get_current_centre_of_mass(self) -> Vector:
        """
        Returns a vector representing the current centre of mass of all positioned atoms

        Returns
        -------
        total: Vector instance
        """
        total = Vector(0, 0)
        count = 0

        for atom in self.structure.graph:
            if atom.draw.positioned:
                total.add(atom.draw.position)
                count += 1

        total.divide(count)

        return total

    @staticmethod
    def atoms_are_in_same_ring(atom_1: Atom, atom_2: Atom) -> bool:
        """
        Returns True if atom_1 and atom_2 are in at least one shared ring, False otherwise

        Parameters
        ----------
        atom_1: Atom instance
        atom_2: Atom instance

        Returns
        -------
        bool, True if atom_1 and atom_2 share a ring, False otherwise
        """
        for ring_id_1 in atom_1.draw.rings:
            for ring_id_2 in atom_2.draw.rings:
                if ring_id_1 == ring_id_2:
                    return True
        return False

    def define_rings(self) -> None:
        """
        Characterise rings and ring systems for drawing purposes
        """

        rings = self.structure.cycles.find_sssr()

        if not rings:
            return

        for ring_members in rings:
            ring = Ring(ring_members)
            self.add_ring(ring)

            for atom in ring_members:
                structure_atom = self.structure.get_atom(atom)
                structure_atom.draw.rings.append(ring.id)

        # Define ring overlaps

        for i, ring_1 in enumerate(self.rings[:-1]):
            for ring_2 in self.rings[i + 1:]:
                ring_overlap = RingOverlap(ring_1, ring_2)

                if len(ring_overlap.atoms) > 0:
                    self.add_ring_overlap(ring_overlap)

        # Identify neighbouring rings

        for ring in self.rings:
            neighbouring_rings = find_neighbouring_rings(self.ring_overlaps, ring.id)
            ring.neighbouring_rings = neighbouring_rings

        # Anchor rings to an atom

        for ring in self.rings:
            anchor = ring.members[0]
            if ring not in anchor.draw.anchored_rings:
                anchor.draw.anchored_rings.append(ring)

        # Prior to merging ring systems, back-up the ring information such that it can be retrieved for double bond
        # and aromatic bond positioning later

        self.backup_ring_info()

        while True:
            ring_id: int = -1

            # If there is a ring that is part of a bridged system, created the system

            for ring in self.rings:
                if self.is_part_of_bridged_ring(ring.id) and not ring.bridged:
                    ring_id: int = ring.id

            if ring_id == -1:
                break

            ring: Ring = self.id_to_ring[ring_id]
            involved_ring_ids: Union[list[int], set[int]] = []
            self.get_bridged_ring_subrings(ring.id, involved_ring_ids)
            involved_ring_ids = set(involved_ring_ids)

            self.has_bridged_ring = True
            self.create_bridged_ring(involved_ring_ids)

            # Move the ring ids that participate in the bridged rings to the bridged systems, and remove them as
            # individual rings to prevent double positioning

            for involved_ring_id in involved_ring_ids:
                involved_ring = self.id_to_ring[involved_ring_id]
                self.remove_ring(involved_ring)

        # TODO: The stuff below seems to do the same as the stuff above? See if it is removable

        # This is new - if stuff breaks, remove
        bridged_systems = find_bridged_systems(self.rings, self.ring_overlaps)
        if bridged_systems and not self.has_bridged_ring:
            self.has_bridged_ring = True
            for bridged_system in bridged_systems:
                involved_ring_ids = set(bridged_system)
                self.create_bridged_ring(involved_ring_ids)
                for involved_ring_id in involved_ring_ids:
                    involved_ring = self.id_to_ring[involved_ring_id]
                    self.remove_ring(involved_ring)

    # TODO: This function seems to do a lot that it doesn't need to, and doesn't do what it needs to! Revisit!
    def hide_hydrogens(self) -> None:
        """
        Mark explicit hydrogens as drawn and implicit hydrogens as not drawn. Store these hydrogens in the adjacent,
        drawn atoms.
        """
        hidden: list[Atom] = []
        exposed: list[Atom] = []

        self.structure.refresh_structure()

        for atom in self.structure.graph:
            if atom.type != 'H':
                continue

            elif atom.charge != 0:
                continue

            neighbour = atom.neighbours[0]
            neighbour.draw.has_hydrogen = True
            atom.draw.is_drawn = False
            hidden.append(atom)

            if len(neighbour.draw.rings) < 2 and neighbour.draw.bridged_ring is None and \
                    neighbour.draw.bridged_ring is not None and len(neighbour.draw.original_rings) < 2:

                atom.draw.is_drawn = False
                neighbour.draw.has_hydrogen = True
                hidden.append(atom)
            else:
                exposed.append(atom)

        for atom in self.structure.graph:
            atom.set_drawn_neighbours()
            if atom.type == 'O':
                pass

        self.drawn_bonds = []

        for bond_nr, bond in self.structure.bonds.items():
            if bond.atom_1.draw.is_drawn and bond.atom_2.draw.is_drawn:
                self.drawn_bonds.append(bond)

        self.drawn_atoms = self.structure.get_drawn_atoms()

    def get_bridged_rings(self) -> list[Ring]:
        """
        Returns a list of bridged ring systems

        Returns
        -------
        bridged_rings: list of [ring, ->], with each ring a Ring instance of a bridged ring system, meaning they
            contain at least two rings which share an overlap of 2 or more bonds

        """
        bridged_rings: list[Ring] = []
        for ring in self.rings:
            if ring.bridged:
                bridged_rings.append(ring)

        return bridged_rings

    def get_bridged_ring_subrings(self, ring_id: int, involved_ring_ids: list[int]) -> None:
        """
        Edit list of ring IDs until it contains all rings that are involved in a bridged ring system with the original
            ring

        Parameters
        ----------
        ring_id: int, ring identifier of original ring
        involved_ring_ids: list of [ring_id, ->], with ring_id int, edited until it contains all rings that are part of
            a bridged system with the ring of label ring_id
        """
        involved_ring_ids.append(ring_id)
        ring = self.id_to_ring[ring_id]

        for neighbour_id in ring.neighbouring_rings:
            if neighbour_id not in involved_ring_ids and neighbour_id != ring_id and \
                    rings_connected_by_bridge(self.ring_overlaps, ring_id, neighbour_id):
                self.get_bridged_ring_subrings(neighbour_id, involved_ring_ids)

    def create_bridged_ring(self, involved_ring_ids: set[int]) -> None:
        """
        Create and store a bridged ring system from a list of involved ring identifiers

        Parameters
        ----------
        involved_ring_ids: set of [ring_id, ->], with each ring_id int. Identifiers of rings involved in the bridged
            ring system
        """
        # Keeps track of all atoms involved in the ring system
        atoms: set[Atom] = set()

        # Keeps track of all rings that neighbour the bridged ring system
        neighbours: set[int] = set()

        for ring_id in involved_ring_ids:
            ring: Ring = self.id_to_ring[ring_id]
            ring.subring_of_bridged = True

            for atom in ring.members:
                atoms.add(atom)

            for neighbour_id in ring.neighbouring_rings:
                neighbours.add(neighbour_id)

        leftovers = set()
        ring_members: set[Atom] = set()

        for atom in atoms:
            intersect = involved_ring_ids & set(atom.draw.rings)

            if len(atom.draw.rings) == 1 or len(intersect) == 1:
                ring_members.add(atom)
            else:
                leftovers.add(atom)

        for atom in leftovers:
            is_on_ring = False

            for bond in atom.bonds:
                bond_associated_rings = min(len(bond.atom_1.draw.rings),
                                            len(bond.atom_2.draw.rings))
                if bond_associated_rings == 1:
                    is_on_ring = True

            if is_on_ring:
                atom.draw.is_bridge_atom = True
                ring_members.add(atom)
            else:
                atom.draw.is_bridge = True
                ring_members.add(atom)

        bridged_ring = Ring(list(ring_members))
        self.add_ring(bridged_ring)
        bridged_ring.bridged = True
        bridged_ring.neighbouring_rings = list(neighbours)

        for ring_id in involved_ring_ids:
            ring = self.id_to_ring[ring_id]
            bridged_ring.subrings.append(ring.copy())

        for atom in ring_members:
            atom.draw.bridged_ring = bridged_ring.id
            for ring_id in involved_ring_ids:
                if ring_id in atom.draw.rings:
                    atom.draw.rings.remove(ring_id)

            atom.draw.rings.append(bridged_ring.id)

        involved_ring_ids = list(involved_ring_ids)

        for i, ring_id_1 in enumerate(involved_ring_ids):
            for ring_id_2 in involved_ring_ids[i + 1:]:
                self.remove_ring_overlaps_between(ring_id_1, ring_id_2)

        for neighbour_id in neighbours:
            ring_overlaps = self.get_ring_overlaps(neighbour_id, involved_ring_ids)
            for ring_overlap in ring_overlaps:
                ring_overlap.update_other(bridged_ring.id, neighbour_id)

            neighbour = self.id_to_ring[neighbour_id]
            neighbour.neighbouring_rings.append(bridged_ring.id)

    def backup_ring_info(self):
        self.original_rings = []

        for ring in self.rings:
            self.original_rings.append(ring.copy())

        self.original_ring_overlaps = []
        for ring_overlap in self.ring_overlaps:
            self.original_ring_overlaps.append(ring_overlap.copy())

        for atom in self.structure.graph:

            atom.draw.original_rings = []
            for ring in atom.draw.rings:
                atom.draw.original_rings.append(ring)

    def get_ring_index(self, ring_id):
        for i, ring in enumerate(self.rings):
            if ring.id == ring_id:
                return i

    def get_ring(self, ring_id):
        for ring in self.rings:
            if ring.id == ring_id:
                return ring

    def add_ring(self, ring):
        ring.id = self.ring_id_tracker
        self.rings.append(ring)
        self.id_to_ring[ring.id] = ring

        self.ring_id_tracker += 1

    def remove_ring(self, ring):
        self.rings.remove(ring)
        overlaps_to_remove = []
        for ring_overlap in self.ring_overlaps:
            if ring_overlap.ring_id_1 == ring.id or ring_overlap.ring_id_2 == ring.id:
                overlaps_to_remove.append(ring_overlap)

        for ring_overlap in overlaps_to_remove:
            self.ring_overlaps.remove(ring_overlap)

        for neighbouring_ring in self.rings:
            if ring.id in neighbouring_ring.neighbouring_rings:
                neighbouring_ring.neighbouring_rings.remove(ring.id)

    def add_ring_overlap(self, ring_overlap):
        ring_overlap.id = self.ring_overlap_id_tracker
        self.ring_overlaps.append(ring_overlap)
        self.ring_overlap_id_tracker += 1

    def get_ring_overlaps(self, ring_id, ring_ids):
        ring_overlaps = []

        for ring_overlap in self.ring_overlaps:
            for ring_id_2 in ring_ids:
                if (ring_overlap.ring_id_1 == ring_id and ring_overlap.ring_id_2 == ring_id_2) or\
                        (ring_overlap.ring_id_2 == ring_id and ring_overlap.ring_id_1 == ring_id_2):
                    ring_overlaps.append(ring_overlap)

        return ring_overlaps

    def remove_ring_overlaps_between(self, ring_id_1, ring_id_2):
        to_remove = []

        for ring_overlap in self.ring_overlaps:
            if (ring_overlap.ring_id_1 == ring_id_1 and ring_overlap.ring_id_2 == ring_id_2) or\
                    (ring_overlap.ring_id_2 == ring_id_1 and ring_overlap.ring_id_1 == ring_id_2):
                to_remove.append(ring_overlap)

        for ring_overlap in to_remove:
            self.ring_overlaps.remove(ring_overlap)

    def is_part_of_bridged_ring(self, ring_id):
        for ring_overlap in self.ring_overlaps:
            if ring_overlap.involves_ring(ring_id) and ring_overlap.is_bridge():
                return True

        return False

    def get_closest_atom(self, atom):
        minimal_distance = 9999999
        closest_atom = None

        for atom_2 in self.drawn_atoms:
            if atom == atom_2:
                continue

            squared_distance = atom.draw.position.get_squared_distance(atom_2.draw.position)

            if squared_distance < minimal_distance:
                minimal_distance = squared_distance
                closest_atom = atom_2

        return closest_atom

    @staticmethod
    def get_sorted_distances_from_list(atom, atom_list):
        atom_distances = []

        for atom_2 in atom_list:
            if atom == atom_2:
                continue

            squared_distance = atom.draw.position.get_squared_distance(atom_2.draw.position)
            atom_distances.append((squared_distance, atom_2))

        atom_distances.sort(key=lambda x: x[0])

        return atom_distances


def draw_multiple(structure: Structure, coords_only: bool = False, options: Union[None, Options] = None, kekulise=True) -> Drawer:
    if not options:
        options = Options()
    options_main = Options()
    options_main.finetune = False

    drawer = Drawer(structure, options=options_main, coords_only=True, multiple=True)
    structures = structure.split_disconnected_structures()
    max_x = -100000000

    for i, substructure in enumerate(structures):
        subdrawer = Drawer(substructure, options=options, coords_only=True, kekulise=kekulise)
        bounding_box = subdrawer.structure.get_bounding_box()
        min_x = bounding_box[0]
        diff_x = max_x - min_x

        if max_x != -100000000:
            subdrawer.move_structure(x=diff_x + 20)

        bounding_box = subdrawer.structure.get_bounding_box()
        max_x = bounding_box[2]

        drawer.chiral_bonds += subdrawer.chiral_bonds
        drawer.chiral_bond_to_orientation.update(subdrawer.chiral_bond_to_orientation)

        for atom in subdrawer.structure.graph:
            if atom.draw.is_drawn:
                atom_2 = drawer.structure.atoms[atom.nr]
                atom_2.draw.position = atom.draw.position
                atom_2.draw.positioned = True

    drawer.structure.refresh_structure()

    if not coords_only:
        drawer.plot_structure()

    return drawer

