#!/usr/bin/env python
from typing import Union, List, Tuple, Dict, Set, Generator
import copy
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
from pikachu.chem.structure import Structure
from pikachu.errors import DrawingError


class Options:
    def __init__(self):
        self.width = 500
        self.height = 500
        self.bond_thickness = 2
        self.bond_length = 15
        self.chiral_bond_width = self.bond_length * 0.1
        self.bond_length_squared = self.bond_length ** 2
        self.short_bond_length = 0.50
        self.double_bond_length = 0.80
        self.bond_spacing = 0.18 * self.bond_length
        self.isomeric = True
        self.padding = 30
        self.font_size_large = 5
        self.font_size_small = 3
        self.kk_threshold = 0.1
        self.kk_inner_threshold = 0.1
        self.kk_max_iteration = 2000
        self.kk_max_inner_iteration = 50
        self.kk_max_energy = 1e9
        self.overlap_sensitivity = 0.10
        self.overlap_resolution_iterations = 5
        self.background_color = 'white'
        self.draw_hydrogens = False
        self.finetune = True
        self.strict_mode = False
        self.svg_font = "verdana"
        self.svg_font_size = 8
        self.svg_font_size_small = 6
        self.svg_letter_spacing = -2


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

        a = 0.0

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
    def __init__(self, structure: Structure, options: Union[Options, None] = None,
                 coords_only: bool = False, multiple: bool = False) -> None:

        if options is None:
            self.options = Options()
        else:
            self.options = options

        self.structure = structure.kekulise()
        self.rings = []
        self.ring_overlaps = []
        self.original_rings = []
        self.original_ring_overlaps = []
        self.id_to_ring = {}
        self.drawn_atoms = []
        self.drawn_bonds = []
        self.bridged_ring = False
        self.total_overlap_score = 0
        self.atom_nr_to_atom = {}
        self.chiral_bonds = []
        self.chiral_bond_to_orientation = {}
        self.fixed_chiral_bonds = set()
        self.multiple = multiple

        self.ring_id_tracker = 0
        self.ring_overlap_id_tracker = 0

        self.svg_groups = {}
        self.annotation = None
        self.structure_id = None

        self.svg_style = """<style> line {stroke: black; stroke_width: 1px;} polygon {fill: black;} </style>"""

        self.draw(coords_only=coords_only)

    def set_annotation_for_grouping(self, annotation: str) -> None:
        self.svg_groups = {}
        self.annotation = annotation

    def find_shortest_path(self, atom_1: Atom, atom_2: Atom) -> List[Bond]:
        distances = {}
        previous_hop = {}
        unvisited = set()

        for atom in self.structure.graph:
            distances[atom] = float('inf')
            previous_hop[atom] = None
            unvisited.add(atom)

        distances[atom_1] = 0

        while unvisited:

            current_atom = None
            minimum = float('inf')
            for atom in unvisited:
                dist = distances[atom]
                if dist < minimum:
                    current_atom = atom

            unvisited.remove(current_atom)
            if current_atom == atom_2:
                break
            for neighbour in self.structure.graph[current_atom]:
                if neighbour in unvisited:
                    alternative_distance = distances[current_atom] + 1.0

                    if alternative_distance < distances[neighbour]:
                        distances[neighbour] = alternative_distance
                        previous_hop[neighbour] = current_atom

        path_atoms = []
        current_atom = atom_2
        if previous_hop[current_atom] or current_atom == atom_1:
            while current_atom:
                path_atoms.insert(0, current_atom)
                current_atom = previous_hop[current_atom]

        path = []

        for i in range(1, len(path_atoms)):
            atom_1 = path_atoms[i - 1]
            atom_2 = path_atoms[i]
            bond = self.structure.bond_lookup[atom_1][atom_2]
            path.append(bond)

        return path

    def finetune_overlap_resolution(self) -> None:

        if self.total_overlap_score > self.options.overlap_sensitivity:
            clashing_atoms = self.find_clashing_atoms()

            best_bonds = []
            for atom_1, atom_2 in clashing_atoms:
                shortest_path = self.find_shortest_path(atom_1, atom_2)
                rotatable_bonds = []
                distances = []

                for i, bond in enumerate(shortest_path):
                    distance_1 = i
                    distance_2 = len(shortest_path) - i

                    average_distance = len(shortest_path) / 2

                    distance_metric = abs(average_distance - distance_1) + abs(average_distance - distance_2)

                    if self.bond_is_rotatable(bond):
                        rotatable_bonds.append(bond)
                        distances.append(distance_metric)

                best_bond = None
                optimal_distance = float('inf')
                for i, distance in enumerate(distances):
                    if distance < optimal_distance:
                        best_bond = rotatable_bonds[i]
                        optimal_distance = distance

                if best_bond:
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

                    scores = [overlap_score]

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
    def find_ring_neighbour(atom: Atom, bond: Bond) -> Atom:
        rings = set(bond.atom_1.draw.rings).intersection(set(bond.atom_2.draw.rings))

        cyclic_neighbour = None

        for neighbour in atom.neighbours:
            if len(set(neighbour.draw.rings).intersection(rings)) > 0 and neighbour.draw.is_drawn and neighbour != bond.get_connected_atom(atom):
                cyclic_neighbour = neighbour
                break

        assert cyclic_neighbour

        return cyclic_neighbour

    def find_ring_branch_to_flip(self, bond: Bond, neighbours_1: List[Atom],
                                 neighbours_2: List[Atom]) -> Tuple[Union[None, Atom],
                                                                   Union[None, List[Atom]]]:

        rings = set(bond.atom_1.draw.rings).intersection(set(bond.atom_2.draw.rings))

        resolvable = True

        if len(neighbours_1) == 1:
            central_atom = bond.atom_1
            flanking_atoms = (neighbours_1[0], bond.atom_2)

        elif len(neighbours_2) == 1:
            central_atom = bond.atom_2
            flanking_atoms = (neighbours_2[0], bond.atom_1)

        else:

            subtree_1_size = None
            neighbour_1_in_cycle = False
            neighbour_1 = None

            subtree_2_size = None
            neighbour_2_in_cycle = False
            neighbour_2 = None

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
            # branch gets flipped

            if not neighbour_1_in_cycle and not neighbour_2_in_cycle:
                if subtree_2_size > subtree_1_size:
                    central_atom = bond.atom_1
                    ring_atom = self.find_ring_neighbour(bond.atom_1, bond)
                    flanking_atoms = (bond.atom_2, ring_atom)

                else:
                    central_atom = bond.atom_2
                    ring_atom = self.find_ring_neighbour(bond.atom_2, bond)
                    flanking_atoms = (bond.atom_1, ring_atom)

            elif neighbour_1_in_cycle and not neighbour_2_in_cycle:
                central_atom = bond.atom_2
                ring_atom = self.find_ring_neighbour(bond.atom_2, bond)
                flanking_atoms = (bond.atom_1, ring_atom)

            elif neighbour_2_in_cycle and not neighbour_1_in_cycle:
                central_atom = bond.atom_1
                ring_atom = self.find_ring_neighbour(bond.atom_1, bond)
                flanking_atoms = (bond.atom_2, ring_atom)

            else:
                central_atom = None
                flanking_atoms = None
                resolvable = False

        if resolvable:
            return central_atom, flanking_atoms
        else:
            return None, None

    def flip_stereobond_in_ring(self, bond: Bond) -> None:
        neighbours_1 = bond.atom_1.get_drawn_neighbours()[:]
        neighbours_2 = bond.atom_2.get_drawn_neighbours()[:]

        neighbours_1.remove(bond.atom_2)
        neighbours_2.remove(bond.atom_1)

        # get the rings that this bond is part of

        # rings = set(bond.atom_1.draw.rings).intersection(set(bond.atom_2.draw.rings))

        central_atom = None

        # Atoms that are members of stereobonds have at least two neighbours: the atom at the other side
        # of the bond, and an atom that determines the cis/trans configuration.

        resolvable = True
        flanking_atoms = []

        # Check if the neighbouring atoms are adjacent to stereobonds, and if those
        # stereobonds have already been fixed

        neighbours_1_adjacent_to_stereobond = False
        neighbours_2_adjacent_to_stereobond = False

        neighbours_1_adjacent_to_fixed_stereobond = False
        neighbours_2_adjacent_to_fixed_stereobond = False

        for neighbour_1 in neighbours_1:
            if neighbour_1.adjacent_to_stereobond():
                neighbours_1_adjacent_to_stereobond = True
                for bond_1 in neighbour_1.bonds:
                    if bond_1.chiral and bond_1 in self.fixed_chiral_bonds:
                        neighbours_1_adjacent_to_fixed_stereobond = True

        for neighbour_2 in neighbours_2:
            if neighbour_2.adjacent_to_stereobond():
                neighbours_2_adjacent_to_stereobond = True
                for bond_2 in neighbour_2.bonds:
                    if bond_2.chiral and bond_2 in self.fixed_chiral_bonds:
                        neighbours_2_adjacent_to_fixed_stereobond = True

        if not neighbours_1_adjacent_to_stereobond and not neighbours_2_adjacent_to_stereobond:

            central_atom, flanking_atoms = self.find_ring_branch_to_flip(bond, neighbours_1, neighbours_2)
            if not central_atom:
                resolvable = False

        if neighbours_1_adjacent_to_stereobond and not neighbours_2_adjacent_to_stereobond:
            central_atom = bond.atom_2
            ring_neighbour = self.find_ring_neighbour(bond.atom_2, bond)
            flanking_atoms = (bond.atom_1, ring_neighbour)

        elif neighbours_2_adjacent_to_stereobond and not neighbours_1_adjacent_to_stereobond:
            central_atom = bond.atom_1
            ring_neighbour = self.find_ring_neighbour(bond.atom_1, bond)
            flanking_atoms = (bond.atom_2, ring_neighbour)

        elif neighbours_1_adjacent_to_stereobond and neighbours_2_adjacent_to_stereobond:
            if neighbours_1_adjacent_to_fixed_stereobond and not neighbours_2_adjacent_to_fixed_stereobond:
                central_atom = bond.atom_2
                ring_neighbour = self.find_ring_neighbour(bond.atom_2, bond)
                flanking_atoms = (bond.atom_1, ring_neighbour)

            elif neighbours_2_adjacent_to_fixed_stereobond and not neighbours_1_adjacent_to_fixed_stereobond:
                central_atom = bond.atom_1
                ring_neighbour = self.find_ring_neighbour(bond.atom_1, bond)
                flanking_atoms = (bond.atom_2, ring_neighbour)
            elif not neighbours_1_adjacent_to_fixed_stereobond and not neighbours_2_adjacent_to_fixed_stereobond:
                central_atom, flanking_atoms = self.find_ring_branch_to_flip(bond, neighbours_1, neighbours_2)
                if not central_atom:
                    resolvable = False
            else:
                resolvable = False

        if resolvable:
            self.flip_subtree(central_atom, flanking_atoms[0], flanking_atoms[1])

        else:
            if self.options.strict_mode:
                raise DrawingError('chiral bond ring')
            else:
                print("Warning! Cis/trans stereochemistry of cyclic system incorrectly drawn.")
            
    def flip_subtree(self, root: Atom, atom_1: Atom, atom_2: Atom) -> None:

        for atom in self.traverse_substructure(root, {atom_1, atom_2}):
            atom.draw.position.mirror_about_line(atom_1.draw.position, atom_2.draw.position)
            for anchored_ring in atom.draw.anchored_rings:
                if anchored_ring.center:
                    anchored_ring.center.mirror_about_line(atom_1.draw.position, atom_2.draw.position)

    def find_clashing_atoms(self) -> List[Tuple[Atom, Atom]]:
        clashing_atoms = []
        for i, atom_1 in enumerate(self.drawn_atoms):
            for j in range(i + 1, len(self.drawn_atoms)):
                atom_2 = self.drawn_atoms[j]
                if not self.structure.bond_exists(atom_1, atom_2):
                    distance = Vector.subtract_vectors(atom_1.draw.position, atom_2.draw.position).get_squared_length()
                    if distance < self.options.bond_length_squared:
                        clashing_atoms.append((atom_1, atom_2))

        return clashing_atoms

    def prioritise_chiral_bonds(self, chiral_center: Atom) -> List[Atom]:

        subtree_1_size = self.get_subgraph_size(chiral_center.neighbours[0], {chiral_center})
        subtree_2_size = self.get_subgraph_size(chiral_center.neighbours[1], {chiral_center})
        subtree_3_size = self.get_subgraph_size(chiral_center.neighbours[2], {chiral_center})

        sizes_and_atoms = [(subtree_1_size, chiral_center.neighbours[0]),
                           (subtree_2_size, chiral_center.neighbours[1]),
                           (subtree_3_size, chiral_center.neighbours[2])]

        if len(chiral_center.neighbours) == 4:
            subtree_4_size = self.get_subgraph_size(chiral_center.neighbours[3], {chiral_center})

            sizes_and_atoms = [(subtree_1_size, chiral_center.neighbours[0]),
                               (subtree_2_size, chiral_center.neighbours[1]),
                               (subtree_3_size, chiral_center.neighbours[2]),
                               (subtree_4_size, chiral_center.neighbours[3])]

        sizes_and_atoms.sort(key=lambda x: (x[0], ATOM_PROPERTIES.element_to_atomic_nr[x[1].type]))

        options_h = []
        options = []
        backup_options_rings = []
        backup_options_chiral_noring = []
        backup_options_chiral_ring = []
        backup_options_chiral_neighbour = []
        backup_options_chiral_neighbour_ring = []
        non_options = []

        for neighbour in [size_and_atom[1] for size_and_atom in sizes_and_atoms]:
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

                elif not other_chiral_centre and not in_ring and not neighbour.chiral:
                    options.append(neighbour)

                elif other_chiral_centre and not neighbour.chiral and not in_ring:
                    backup_options_chiral_noring.append(neighbour)

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
    def place_eclipsed_bond(hydrogen: Atom, angles_between_lines: List[float], atom_order: List[Atom],
                            wedge_atom: Atom) -> None:
        position = None

        for i, angle in enumerate(angles_between_lines):

            if round(angle, 3) >= round(math.pi, 3):
                position = i

        if position is not None:
            atom_order.insert(position, hydrogen)
        else:
            wedge_index = atom_order.index(wedge_atom)
            atom_order.insert(wedge_index + 1, hydrogen)

    @staticmethod
    def reorder(atom_order: List[Atom], wedge_atom: Atom) -> List[Atom]:
        first_index = atom_order.index(wedge_atom)
        new_order = atom_order[first_index:] + atom_order[:first_index]
        return new_order

    def move_structure(self, x: float = 0.0, y: float = 0.0) -> None:
        for atom in self.structure.graph:
            if atom.draw.is_drawn:
                atom.draw.position.x += x
                atom.draw.position.y += y
    
    def determine_chirality(self, chiral_center: Atom) -> Tuple[Bond, str]:

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

        priority = self.prioritise_chiral_bonds(chiral_center)
        
        if len(atom_order) == 3:

            if chiral_center.has_neighbour('H'):
                hydrogen = chiral_center.get_neighbour('H')
                assert not hydrogen.draw.is_drawn
                eclipsed_element = hydrogen
                wedge_atom = priority[1]

            elif chiral_center.lone_pairs:
                lone_pair = chiral_center.lone_pairs[0]
                eclipsed_element = lone_pair
                wedge_atom = priority[0]
            else:
                raise DrawingError("chiral center")

            self.place_eclipsed_bond(eclipsed_element, angles_between_lines, atom_order, wedge_atom)

        else:
            wedge_atom = priority[0]

        assert len(atom_order) == 4

        atom_order = self.reorder(atom_order, wedge_atom)

        original_order = chiral_center.neighbours[:]

        if len(original_order) == 3:
            if chiral_center.lone_pairs:
                original_order.append(chiral_center.lone_pairs[0])
            else:
                raise DrawingError("chiral center")

        order_matches_chirality = False

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

        wedge_bond = self.structure.bond_lookup[wedge_atom][chiral_center]

        return wedge_bond, wedge

    def set_chiral_bonds(self) -> None:
        for atom in self.structure.graph:
            if atom.chiral:
                bond, wedge = self.determine_chirality(atom)

                self.chiral_bonds.append(bond)
                self.chiral_bond_to_orientation[bond] = (wedge, atom)

    def flip_y_axis(self) -> None:
        for atom in self.structure.graph:
            if atom.draw.positioned:
                atom.draw.position.y = -atom.draw.position.y

    def convert_to_int(self) -> None:
        for atom in self.structure.graph:
            if atom.draw.positioned:
                atom.draw.position.x = int(atom.draw.position.x)
                atom.draw.position.y = int(atom.draw.position.y)

    def move_to_positive_coords(self) -> None:
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

    def draw_text(self, text: Union[str, int], x: float, y: float, font_size: Union[int, None] = None) -> str:
        if font_size is None:
            font_size = self.options.svg_font_size
        svg_text = f"""<text x="{x}" y="{y}" text-anchor="middle" font-family="{self.options.svg_font}" font-size="{font_size}"><tspan y="{y}" dy="0.35em">{text}</tspan></text>"""

        return svg_text

    def draw_svg(self, annotation: Union[None, str] = None) -> str:

        self.set_annotation_for_grouping(annotation)

        ring_centers_x = []
        ring_centers_y = []

        for ring in self.rings:
            self.set_ring_center(ring)

            ring_centers_x.append(ring.center.x)
            ring_centers_y.append(ring.center.y)

        for bond_nr, bond in self.structure.bonds.items():
            if bond.atom_1.draw.positioned and bond.atom_2.draw.positioned:
                line = Line(bond.atom_1.draw.position, bond.atom_2.draw.position, bond.atom_1, bond.atom_2)
                midpoint = line.get_midpoint()

                if bond.type == 'single':
                    if bond in self.chiral_bonds:

                        orientation, chiral_center = self.chiral_bond_to_orientation[bond]
                        self.draw_chiral_bond(orientation, chiral_center, line, midpoint)
                    else:
                        self.draw_halflines(line, midpoint)
                elif bond.type == 'double':
                    if not self.is_terminal(bond.atom_1) and not self.is_terminal(bond.atom_2):
                        self.draw_halflines(line, midpoint)

                        common_ring_numbers = self.get_common_rings(bond.atom_1, bond.atom_2)

                        if common_ring_numbers:
                            common_rings = []
                            for ring_nr in common_ring_numbers:
                                common_rings.append(self.get_ring(ring_nr))

                            common_rings.sort(key=lambda x: len(x.members))
                            common_ring = common_rings[0]
                            ring_centre = common_ring.center
                            second_line = line.double_line_towards_center(ring_centre, self.options.bond_spacing,
                                                                          self.options.double_bond_length)
                            second_line_midpoint = second_line.get_midpoint()
                            self.draw_halflines_double(second_line, second_line_midpoint)

                        else:
                            bond_neighbours = bond.atom_1.drawn_neighbours + bond.atom_2.drawn_neighbours
                            if bond_neighbours:
                                vectors = [atom.draw.position for atom in bond_neighbours]
                                gravitational_point = Vector.get_average(vectors)
                                second_line = line.double_line_towards_center(gravitational_point,
                                                                              self.options.bond_spacing,
                                                                              self.options.double_bond_length)
                                second_line_midpoint = second_line.get_midpoint()
                                self.draw_halflines_double(second_line, second_line_midpoint)
                            else:
                                print("Shouldn't happen!")
                    else:
                        if self.is_terminal(bond.atom_1) and self.is_terminal(bond.atom_2):
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

                            self.draw_halflines_double(double_bond_line_1, double_bond_line_1_midpoint)
                            self.draw_halflines_double(double_bond_line_2, double_bond_line_2_midpoint)

                        else:

                            if self.is_terminal(bond.atom_1):
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

                                self.draw_halflines(double_bond_line_1, double_bond_line_1_midpoint)
                                self.draw_halflines(double_bond_line_2, double_bond_line_2_midpoint)

                            else:
                                self.draw_halflines(line, midpoint)

                                bond_neighbours = bond.atom_1.drawn_neighbours + bond.atom_2.drawn_neighbours
                                if bond_neighbours:
                                    vectors = [atom.draw.position for atom in bond_neighbours]
                                    gravitational_point = Vector.get_average(vectors)
                                    second_line = line.get_parallel_line(gravitational_point,
                                                                         self.options.bond_spacing)
                                    second_line_midpoint = second_line.get_midpoint()
                                    self.draw_halflines(second_line, second_line_midpoint)
                                else:
                                    print("Shouldn't happen!")

                elif bond.type == 'triple':
                    self.draw_halflines(line, midpoint)
                    line_1, line_2 = line.get_parallel_lines(self.options.bond_spacing)
                    line_1_midpoint = line_1.get_midpoint()
                    line_2_midpoint = line_2.get_midpoint()
                    self.draw_halflines(line_1, line_1_midpoint)
                    self.draw_halflines(line_2, line_2_midpoint)

        for atom in self.structure.graph:
            if atom.draw.positioned:
                svg_text = ''
                svg_h_text = ''
                svg_charge_text = ''
                svg_h_count_text = ''

                if atom.type != 'C' or atom.draw.draw_explicit or atom.charge:
                    if atom.type == 'C' and not atom.charge:
                        svg_text = self.draw_text('.', atom.draw.position.x, atom.draw.position.y - 2)
                    else:
                        svg_text = self.draw_text(atom.type, atom.draw.position.x, atom.draw.position.y)

                        # TODO: Make this possible in svg writing
                        # text = self.set_r_group_indices_subscript(atom.type)

                orientation = self.get_hydrogen_text_orientation(atom)

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
                        svg_h_text = self.draw_text('H', h_pos.x, h_pos.y)
                        if hydrogen_count > 1:
                            svg_h_count_text = self.draw_text(hydrogen_count, h_subscript_pos.x, h_subscript_pos.y,
                                                              font_size=self.options.svg_font_size_small)
                    if atom.charge:
                        if atom.charge > 0:
                            charge_symbol = '+'
                        else:
                            charge_symbol = '-'

                        if abs(atom.charge) > 1:
                            charge_text = f"{abs(atom.charge)}{charge_symbol}"
                        else:
                            charge_text = charge_symbol

                        svg_charge_text = self.draw_text(charge_text, charge_pos.x, charge_pos.y,
                                                         font_size=self.options.svg_font_size_small)

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
                    self.add_svg_element(text_group, atom)

        svg = self.assemble_svg()
        return svg

    def write_svg(self, out_file: str, annotation: Union[str, None] = None) -> None:

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
        svg_string += self.draw_svg(annotation=annotation)
        svg_string += "</svg>"

        with open(out_file, 'w') as out:
            out.write(svg_string)

    def assemble_svg(self) -> str:

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

    def draw(self, coords_only: bool = False) -> None:

        if not self.options.draw_hydrogens:
            self.hide_hydrogens()

        self.get_atom_nr_to_atom()
        self.define_rings()

        if not self.multiple:
            self.process_structure()
            self.set_chiral_bonds()
            if not coords_only:
                self.draw_structure()
        else:
            self.restore_ring_information()

    @staticmethod
    def get_hydrogen_text_orientation(atom: Atom) -> str:
        four_positions = [Vector(atom.draw.position.x, atom.draw.position.y + 3),
                          Vector(atom.draw.position.x, atom.draw.position.y - 3),
                          Vector(atom.draw.position.x + 3, atom.draw.position.y),
                          Vector(atom.draw.position.x - 3, atom.draw.position.y)]

        positions_to_angles = [[], [], [], []]

        for neighbour in atom.drawn_neighbours:
            for i, position in enumerate(four_positions):
                angle = Vector.get_angle_between_vectors(position, neighbour.draw.position, atom.draw.position)
                positions_to_angles[i].append(angle)

        orientation = None

        if not positions_to_angles[0]:
            orientation = 'H_after_atom'
        elif min(positions_to_angles[2]) > 1.57078 or min(positions_to_angles[3]) > 1.57078:
            if min(positions_to_angles[2]) >= min(positions_to_angles[3]):
                orientation = 'H_after_atom'

            else:
                orientation = 'H_before_atom'
        else:

            smallest_angles = [min(angles) for angles in positions_to_angles]

            angle = 0
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

        return orientation

    @staticmethod
    def in_same_ring(atom_1: Atom, atom_2: Atom) -> bool:
        if atom_1.draw.rings and atom_2.draw.rings:
            joined_rings = list(set(atom_1.draw.rings).intersection(set(atom_2.draw.rings)))
            if joined_rings:
                return True

        return False

    @staticmethod
    def get_common_rings(atom_1: Atom, atom_2: Atom) -> Union[List[int], None]:
        if atom_1.draw.rings and atom_2.draw.rings:
            joined_rings = list(set(atom_1.draw.rings).intersection(set(atom_2.draw.rings)))
            return joined_rings

        return None

    @staticmethod
    def plot_chiral_bond_front(polygon: List[Vector], ax: Axes, color: str = 'black') -> None:
        x = []
        y = []
        for point in polygon:
            x.append(point.x)
            y.append(point.y)

        ax.fill(x, y, color=color)

    def plot_chiral_bond_back(self, lines: List[Line], ax: Axes, color: str = 'black') -> None:
        for line in lines:
            self.plot_line(line, ax, color=color)

    def plot_chiral_bond(self, orientation: str, chiral_center: Atom, line: Line, ax: Axes, midpoint: Vector) -> None:
        halflines = line.divide_in_two(midpoint)
        for halfline in halflines:
            truncated_line = halfline.get_truncated_line(self.options.short_bond_length)
            if orientation == 'front':
                bond_wedge = truncated_line.get_bond_wedge_front(self.options.chiral_bond_width, chiral_center)
                self.plot_chiral_bond_front(bond_wedge, ax, color=halfline.atom.draw.colour)
            else:
                bond_lines = halfline.get_bond_wedge_back(self.options.chiral_bond_width, chiral_center)
                self.plot_chiral_bond_back(bond_lines, ax, color=halfline.atom.draw.colour)

    def draw_chiral_bond_front(self, polygon: List[Vector], color: str = 'black') -> str:
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

    def draw_chiral_bond_back(self, lines: List[Line], color: str = 'black') -> str:
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
        self.structure_id = structure_id

    def add_svg_element(self, svg_element: str, atom: Atom) -> None:

        if self.annotation is not None and self.annotation in atom.annotations.annotations:
            annotation_value = atom.annotations.get_attribute(self.annotation)
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

    def draw_chiral_bond(self, orientation: str, chiral_center: Atom, line: Line, midpoint: Vector) -> None:
        halflines = line.divide_in_two(midpoint)
        for halfline in halflines:
            truncated_line = halfline.get_truncated_line(self.options.short_bond_length)
            if orientation == 'front':
                bond_wedge = truncated_line.get_bond_wedge_front(self.options.chiral_bond_width, chiral_center)
                svg_polygon = self.draw_chiral_bond_front(bond_wedge, color=halfline.atom.draw.colour)
                self.add_svg_element(svg_polygon, halfline.atom)

            elif orientation == 'back':
                bond_lines = halfline.get_bond_wedge_back(self.options.chiral_bond_width, chiral_center)
                svg_lines = self.draw_chiral_bond_back(bond_lines, color=halfline.atom.draw.colour)
                self.add_svg_element(svg_lines, halfline.atom)
            else:
                raise ValueError(f"Unrecognised chiral bond orientation: {orientation}.")

    def draw_halflines(self, line: Line, midpoint: Vector) -> None:

        halflines = line.divide_in_two(midpoint)
        for halfline in halflines:
            truncated_line = halfline.get_truncated_line(self.options.short_bond_length)
            svg_line = self.draw_line(truncated_line, color=halfline.atom.draw.colour)
            self.add_svg_element(svg_line, halfline.atom)

    def draw_halflines_double(self, line: Line, midpoint: Vector) -> None:

        halflines = line.divide_in_two(midpoint)
        for halfline in halflines:
            svg_line = self.draw_line(halfline, color=halfline.atom.draw.colour)
            self.add_svg_element(svg_line, halfline.atom)

    def draw_line(self, line: Line, color: str = 'black') -> str:
        if color != 'black':
            svg_line = f'<line x1="{line.point_1.x}" y1="{line.point_1.y}" x2="{line.point_2.x}" y2="{line.point_2.y}" stroke="{color}"/>'
        else:
            svg_line = f'<line x1="{line.point_1.x}" y1="{line.point_1.y}" x2="{line.point_2.x}" y2="{line.point_2.y}" />'
        return svg_line

    def plot_halflines(self, line: Line, ax: Axes, midpoint: Vector) -> None:
        halflines = line.divide_in_two(midpoint)
        for halfline in halflines:
            truncated_line = halfline.get_truncated_line(self.options.short_bond_length)
            self.plot_line(truncated_line, ax, color=halfline.atom.draw.colour)

    def plot_halflines_double(self, line: Line, ax: Axes, midpoint: Vector) -> None:
        halflines = line.divide_in_two(midpoint)
        for halfline in halflines:
            self.plot_line(halfline, ax, color=halfline.atom.draw.colour)

    def plot_line(self, line: Line, ax: Axes, color: str = 'black') -> None:
        ax.plot([line.point_1.x, line.point_2.x],
                [line.point_1.y, line.point_2.y], color=color, linewidth=self.options.bond_thickness)

    @staticmethod
    def get_image_as_array() -> np.ndarray:
        # Return image as np.ndarray that represents RGB image
        canvas = plt.gca().figure.canvas
        canvas.draw()
        image = np.frombuffer(canvas.tostring_rgb(), dtype='uint8')
        image = image.reshape(canvas.get_width_height()[::-1] + (3,))
        plt.close('all')
        return image

    @staticmethod
    def save_svg(out_file: str) -> None:
        if out_file.endswith('.svg'):
            pass
        else:
            out_file += '.svg'
        plt.savefig(out_file)
        plt.clf()
        plt.close(plt.gcf())
        plt.close('all')

    @staticmethod
    def save_svg_string() -> str:
        svg_string = StringIO()
        plt.savefig(svg_string, format='svg')
        svg = svg_string.getvalue()

        plt.clf()
        plt.close(plt.gcf())
        plt.close('all')
        return svg

    @staticmethod
    def save_png(out_file: str) -> None:
        if out_file.endswith('.png'):
            pass
        else:
            out_file += '.png'
        plt.savefig(out_file)
        plt.clf()
        plt.close()

    @staticmethod
    def show_molecule() -> None:
        plt.show()
        plt.clf()
        plt.close()

    @staticmethod
    def chirality_correct(bond: Bond) -> bool:
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

    def fix_chiral_bond(self, double_bond: Bond) -> None:
        if len(double_bond.atom_1.draw.rings) and len(double_bond.atom_2.draw.rings) and \
                len(set(double_bond.atom_1.draw.rings).intersection(set(double_bond.atom_2.draw.rings))) >= 1:
            self.flip_stereobond_in_ring(double_bond)

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
                self.flip_subtree(neighbour, root_atom, parent_atom)

            # Only need to flip once if both neighbours are in the same ring

            elif len(neighbours) == 2 and len(
                    set(neighbours[0].draw.rings).intersection(set(neighbours[1].draw.rings))) >= 1:

                self.flip_subtree(neighbours[0], root_atom, parent_atom)

            elif len(neighbours) == 2:
                neighbour_1 = neighbours[0]
                neighbour_2 = neighbours[1]

                self.flip_subtree(neighbour_1, root_atom, parent_atom)
                self.flip_subtree(neighbour_2, root_atom, parent_atom)

            self.fixed_chiral_bonds.add(double_bond)

    def fix_chiral_bonds_in_rings(self) -> None:
        double_bond_sequences = self.structure.find_double_bond_sequences()

        for double_bond_sequence in double_bond_sequences:
            for double_bond in double_bond_sequence:
                chirality_correct = self.chirality_correct(double_bond)
                if chirality_correct:
                    self.fixed_chiral_bonds.add(double_bond)
                else:
                    self.fix_chiral_bond(double_bond)

        for bond in self.structure.bonds.values():
            if bond.type == 'double' and bond.chiral and bond not in self.fixed_chiral_bonds:
                chirality_correct = self.chirality_correct(bond)
                if chirality_correct:
                    self.fixed_chiral_bonds.add(bond)
                else:
                    self.fix_chiral_bond(bond)

    def draw_structure(self) -> None:

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

        height = max_y - min_y
        width = max_x - min_x

        fig, ax = plt.subplots(figsize=((width + 2 * self.options.padding) / 50.0,
                                        (height + 2 * self.options.padding) / 50.0), dpi=100)

        ax.set_aspect('equal', adjustable='box')
        ax.axis('off')

        ax.set_xlim([min_x - self.options.padding, max_x + self.options.padding])
        ax.set_ylim([min_y - self.options.padding, max_y + self.options.padding])
        plt.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)

        params = {'mathtext.default': 'regular'}
        plt.rcParams.update(params)

        ring_centers_x = []
        ring_centers_y = []

        for ring in self.rings:
            self.set_ring_center(ring)

            ring_centers_x.append(ring.center.x)
            ring_centers_y.append(ring.center.y)

        for bond_nr, bond in self.structure.bonds.items():
            if bond.atom_1.draw.positioned and bond.atom_2.draw.positioned:
                line = Line(bond.atom_1.draw.position, bond.atom_2.draw.position, bond.atom_1, bond.atom_2)
                midpoint = line.get_midpoint()

                if bond.type == 'single':
                    if bond in self.chiral_bonds:
                        orientation, chiral_center = self.chiral_bond_to_orientation[bond]
                        self.plot_chiral_bond(orientation, chiral_center, line, ax, midpoint)
                    else:
                        self.plot_halflines(line, ax, midpoint)
                elif bond.type == 'double':
                    if not self.is_terminal(bond.atom_1) and not self.is_terminal(bond.atom_2):
                        self.plot_halflines(line, ax, midpoint)

                        common_ring_numbers = self.get_common_rings(bond.atom_1, bond.atom_2)

                        if common_ring_numbers:
                            common_rings = []
                            for ring_nr in common_ring_numbers:
                                common_rings.append(self.get_ring(ring_nr))

                            common_rings.sort(key=lambda x: len(x.members))
                            common_ring = common_rings[0]
                            ring_centre = common_ring.center
                            second_line = line.double_line_towards_center(ring_centre, self.options.bond_spacing, self.options.double_bond_length)
                            second_line_midpoint = second_line.get_midpoint()
                            self.plot_halflines_double(second_line, ax, second_line_midpoint)

                        else:
                            bond_neighbours = bond.atom_1.drawn_neighbours + bond.atom_2.drawn_neighbours
                            if bond_neighbours:
                                vectors = [atom.draw.position for atom in bond_neighbours]
                                gravitational_point = Vector.get_average(vectors)
                                second_line = line.double_line_towards_center(gravitational_point, self.options.bond_spacing, self.options.double_bond_length)
                                second_line_midpoint = second_line.get_midpoint()
                                self.plot_halflines_double(second_line, ax, second_line_midpoint)
                            else:
                                print("Shouldn't happen!")
                    else:
                        if self.is_terminal(bond.atom_1) and self.is_terminal(bond.atom_2):
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

                            self.plot_halflines_double(double_bond_line_1, ax, double_bond_line_1_midpoint)
                            self.plot_halflines_double(double_bond_line_2, ax, double_bond_line_2_midpoint)

                        else:

                            if self.is_terminal(bond.atom_1):
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

                                self.plot_halflines(double_bond_line_1, ax, double_bond_line_1_midpoint)
                                self.plot_halflines(double_bond_line_2, ax, double_bond_line_2_midpoint)

                            else:
                                self.plot_halflines(line, ax, midpoint)

                                bond_neighbours = bond.atom_1.drawn_neighbours + bond.atom_2.drawn_neighbours
                                if bond_neighbours:
                                    vectors = [atom.draw.position for atom in bond_neighbours]
                                    gravitational_point = Vector.get_average(vectors)
                                    second_line = line.get_parallel_line(gravitational_point,
                                                                         self.options.bond_spacing)
                                    second_line_midpoint = second_line.get_midpoint()
                                    self.plot_halflines(second_line, ax, second_line_midpoint)
                                else:
                                    print("Shouldn't happen!")

                elif bond.type == 'triple':
                    self.plot_halflines(line, ax, midpoint)
                    line_1, line_2 = line.get_parallel_lines(self.options.bond_spacing)
                    line_1_midpoint = line_1.get_midpoint()
                    line_2_midpoint = line_2.get_midpoint()
                    self.plot_halflines(line_1, ax, line_1_midpoint)
                    self.plot_halflines(line_2, ax, line_2_midpoint)

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

                orientation = self.get_hydrogen_text_orientation(atom)
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

    @staticmethod
    def set_r_group_indices_subscript(atom_text: str) -> str:
        # Take str and return the same str with subscript digits
        # (pattern is necessary to not to get confused with isotopes)
        sub_translation = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
        match = re.search('[RXZ]\d+', atom_text)
        if match:
            matched_pattern = match.group()
            adapted_pattern = matched_pattern.translate(sub_translation)
            atom_text = atom_text.replace(matched_pattern, adapted_pattern)
        return atom_text

    @staticmethod
    def is_terminal(atom: Atom) -> bool:
        if len(atom.drawn_neighbours) <= 1:
            return True

        return False

    def process_structure(self) -> None:
        self.position()
        self.structure.refresh_structure()
        self.restore_ring_information()

        self.fix_chiral_bonds_in_rings()

        self.resolve_primary_overlaps()

        self.total_overlap_score, sorted_overlap_scores, atom_to_scores = self.get_overlap_score()

        for i in range(self.options.overlap_resolution_iterations):
            for bond in self.drawn_bonds:
                if self.can_rotate_around_bond(bond):

                    tree_depth_1 = self.get_subgraph_size(bond.atom_1, {bond.atom_2})
                    tree_depth_2 = self.get_subgraph_size(bond.atom_2, {bond.atom_1})

                    # Check neighbouring bonds to ensure neither are chiral; only then the bond is rotatable at the end of the indicated atom

                    atom_1_rotatable = True
                    atom_2_rotatable = True

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
                            angle = neighbour.draw.position.get_rotation_away_from_vector(atom_1.draw.position, atom_2.draw.position, math.radians(120))

                            self.rotate_subtree(neighbour, atom_2, angle, atom_2.draw.position)

                            new_overlap_score, _, _ = self.get_overlap_score()
                            if new_overlap_score > self.total_overlap_score:
                                self.rotate_subtree(neighbour, atom_2, -angle, atom_2.draw.position)
                            else:
                                self.total_overlap_score = new_overlap_score

                        elif len(neighbours_2) == 2:
                            if atom_2.draw.rings and atom_1.draw.rings:
                                continue

                            neighbour_1 = neighbours_2[0]
                            neighbour_2 = neighbours_2[1]

                            if len(neighbour_1.draw.rings) == 1 and len(neighbour_2.draw.rings) == 1:
                                # If the neighbours are in different rings, or in rings at all, do nothing
                                if neighbour_1.draw.rings[0] != neighbour_2.draw.rings[0]:
                                    continue
                            elif neighbour_1.draw.rings or neighbour_2.draw.rings:
                                continue
                            else:
                                angle_1 = neighbour_1.draw.position.get_rotation_away_from_vector(atom_1.draw.position, atom_2.draw.position, math.radians(120))
                                angle_2 = neighbour_2.draw.position.get_rotation_away_from_vector(atom_1.draw.position, atom_2.draw.position, math.radians(120))

                                self.rotate_subtree(neighbour_1, atom_2, angle_1, atom_2.draw.position)
                                self.rotate_subtree(neighbour_2, atom_2, angle_2, atom_2.draw.position)

                                new_overlap_score, _, _ = self.get_overlap_score()

                                if new_overlap_score > self.total_overlap_score:
                                    self.rotate_subtree(neighbour_1, atom_2, -angle_1, atom_2.draw.position)
                                    self.rotate_subtree(neighbour_2, atom_2, -angle_2, atom_2.draw.position)
                                else:
                                    self.total_overlap_score = new_overlap_score

                        self.total_overlap_score, sorted_overlap_scores, atom_to_scores = self.get_overlap_score()

        if self.options.finetune:
            self.finetune_overlap_resolution()
            self.total_overlap_score, sorted_overlap_scores, atom_to_scores = self.get_overlap_score()

        for i in range(self.options.overlap_resolution_iterations):

            self.resolve_secondary_overlaps(sorted_overlap_scores)

    def position(self) -> None:
        start_atom = None

        for atom in self.structure.graph:
            if atom.draw.bridged_ring is not None:
                start_atom = atom
                break

        # is this necessary?

        for ring in self.rings:
            if ring.bridged:
                start_atom = ring.members[0]

        if len(self.rings) > 0 and start_atom is None:
            start_atom = self.rings[0].members[0]

        if start_atom is None:
            start_atom = list(self.drawn_atoms)[0]

        self.create_next_bond(start_atom, None, 0.0)

    def create_next_bond(self, atom, previous_atom=None, angle=0.0,
                         previous_branch_shortest=False, skip_positioning=False):

        if atom.draw.positioned and not skip_positioning:
            return

        if not skip_positioning:
            if not previous_atom:
                dummy = Vector(self.options.bond_length, 0)
                dummy.rotate(math.radians(-60.0))

                atom.draw.previous_position = dummy
                atom.draw.previous_atom = None
                atom.draw.set_position(Vector(self.options.bond_length, 0))
                atom.draw.angle = math.radians(-60.0)

                if atom.draw.bridged_ring is None:
                    atom.draw.positioned = True

            # If the previous atom was part of a ring

            elif len(previous_atom.draw.rings) > 0:
                neighbours = previous_atom.drawn_neighbours
                joined_vertex = None
                # Initialise position to the origin
                position = Vector(0, 0)

                # If the previous atom was not part of a bridged ring and the previous atom was part of more than one ring

                if previous_atom.draw.bridged_ring is None and len(previous_atom.draw.rings) > 1:
                    # Find the vertex adjoining the current bridged ring that is also in both ring systems. This is the
                    # joined vertex.
                    for neighbour in neighbours:
                        if len(set(neighbour.draw.rings) & set(previous_atom.draw.rings)) == len(previous_atom.draw.rings):
                            joined_vertex = neighbour
                            break

                # If there is no joined vertex

                if not joined_vertex:
                    # For each neighbour that is in the same ring:
                    #
                    for neighbour in neighbours:

                        if neighbour.draw.positioned and self.atoms_are_in_same_ring(neighbour, previous_atom):
                            position.add(Vector.subtract_vectors(neighbour.draw.position, previous_atom.draw.position))

                    position.invert()
                    position.normalise()
                    position.multiply_by_scalar(self.options.bond_length)
                    position.add(previous_atom.draw.position)

                else:

                    position = joined_vertex.draw.position.copy()
                    position.rotate_around_vector(math.pi, previous_atom.draw.position)

                atom.draw.set_previous_position(previous_atom)
                atom.draw.set_position(position)
                atom.draw.positioned = True

            else:
                position = Vector(self.options.bond_length, 0)
                position.rotate(angle)
                position.add(previous_atom.draw.position)

                atom.draw.set_position(position)
                atom.draw.set_previous_position(previous_atom)
                atom.draw.positioned = True

        # If the atom is part of a bridged ring: position the entire bridged ring
        if atom.draw.bridged_ring is not None:
            next_ring = self.id_to_ring[atom.draw.bridged_ring]

            if not next_ring.positioned:
                next_center = Vector.subtract_vectors(atom.draw.previous_position, atom.draw.position)
                next_center.invert()
                next_center.normalise()
                scalar = Polygon.find_polygon_radius(self.options.bond_length, len(next_ring.members))

                next_center.multiply_by_scalar(scalar)
                next_center.add(atom.draw.position)

                self.create_ring(next_ring, next_center, atom)

        # If the atom is part of a ring: position the entire ring
        elif len(atom.draw.rings) > 0:
            next_ring = self.id_to_ring[atom.draw.rings[0]]

            if not next_ring.positioned:
                next_center = Vector.subtract_vectors(atom.draw.previous_position, atom.draw.position)
                next_center.invert()
                next_center.normalise()

                radius = Polygon.find_polygon_radius(self.options.bond_length, len(next_ring.members))

                next_center.multiply_by_scalar(radius)
                next_center.add(atom.draw.position)

                self.create_ring(next_ring, next_center, atom)

        # If the atom is not part of a ring, position just the atom
        else:
            neighbours = atom.drawn_neighbours[:]

            if previous_atom:
                if previous_atom in neighbours:
                    neighbours.remove(previous_atom)

            previous_angle = atom.draw.get_angle()

            if len(neighbours) == 1:
                next_atom = neighbours[0]

                current_bond = self.structure.bond_lookup[atom][next_atom]
                previous_bond = None

                if previous_atom:
                    previous_bond = self.structure.bond_lookup[previous_atom][atom]

                if current_bond.type == 'triple' or (previous_bond and previous_bond.type == 'triple') or \
                        (current_bond.type == 'double' and previous_bond and previous_bond.type == 'double' and
                         previous_atom and len(previous_atom.draw.rings) == 0 and
                         len(atom.neighbours) == 2):

                    if current_bond.type == 'double' and previous_bond.type == 'double':

                        atom.draw.draw_explicit = True

                    if previous_atom:
                        previous_bond.draw.center = True

                    current_bond.draw.center = True

                    if current_bond.type == 'double' or current_bond.type == 'triple' or (previous_atom and previous_bond.type == 'triple'):
                        next_atom.draw.angle = 0.0

                    # next_atom.draw.draw_explicit = True

                    self.create_next_bond(next_atom, atom, previous_angle + next_atom.draw.angle)

                elif previous_atom and len(previous_atom.draw.rings) > 0:

                    proposed_angle_1 = math.radians(60.0)
                    proposed_angle_2 = proposed_angle_1 * -1

                    proposed_vector_1 = Vector(self.options.bond_length, 0)
                    proposed_vector_2 = Vector(self.options.bond_length, 0)

                    proposed_vector_1.rotate(proposed_angle_1)
                    proposed_vector_2.rotate(proposed_angle_2)

                    proposed_vector_1.add(atom.draw.position)
                    proposed_vector_2.add(atom.draw.position)

                    centre_of_mass = self.get_current_centre_of_mass()
                    distance_1 = proposed_vector_1.get_squared_distance(centre_of_mass)
                    distance_2 = proposed_vector_2.get_squared_distance(centre_of_mass)

                    if distance_1 < distance_2:
                        next_atom.draw.angle = proposed_angle_2
                    else:
                        next_atom.draw.angle = proposed_angle_1

                    self.create_next_bond(next_atom, atom, previous_angle + next_atom.draw.angle)

                else:

                    a = atom.draw.angle

                    if previous_atom and len(previous_atom.drawn_neighbours) > 3:
                        if round(a, 2) > 0.00:
                            a = min([math.radians(60), a])
                        elif round(a, 2) < 0:
                            a = max([-math.radians(60), a])
                        else:
                            a = math.radians(60)

                    elif not a:
                        last_angled_atom = self.get_last_atom_with_angle(atom)
                        a = last_angled_atom.draw.angle

                        if not a:
                            a = math.radians(60)

                    rotatable = True

                    if previous_atom:

                        bond = self.structure.bond_lookup[previous_atom][atom]
                        if bond.type == 'double' and bond.chiral:
                            rotatable = False

                            previous_previous_atom = previous_atom.draw.previous_atom

                            if previous_previous_atom:

                                configuration = bond.chiral_dict[previous_previous_atom][next_atom]
                                if configuration == 'cis':

                                    a = -a

                    if rotatable:

                        if previous_branch_shortest:
                            next_atom.draw.angle = a
                        else:
                            next_atom.draw.angle = -a
                    else:
                        next_atom.draw.angle = -a

                    if round(math.degrees(next_atom.draw.angle), 0) == 360 or \
                            round(math.degrees(next_atom.draw.angle), 0) == -360 or \
                            round(math.degrees(next_atom.draw.angle), 0) == 0:
                        atom.draw.draw_explicit = True

                    self.create_next_bond(next_atom, atom, previous_angle + next_atom.draw.angle)

            elif len(neighbours) == 2:
                a = atom.draw.angle
                if not a:
                    a = math.radians(60)

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

                trans_atom.draw.angle = a
                cis_atom.draw.angle = -a

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
                                    trans_atom.draw.angle = -a
                                    cis_atom.draw.angle = a

                self.create_next_bond(trans_atom, atom, previous_angle + trans_atom.draw.angle, previous_branch_shortest)
                self.create_next_bond(cis_atom, atom, previous_angle + cis_atom.draw.angle, previous_branch_shortest)

            elif len(neighbours) == 3:
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

    def restore_ring_information(self) -> None:
        bridged_rings = self.get_bridged_rings()

        self.rings = []
        self.ring_overlaps = []

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
        if bond.atom_1.draw.rings and \
                bond.atom_2.draw.rings and \
                len(set(bond.atom_1.draw.rings).intersection(set(bond.atom_2.draw.rings))) > 0:
            return False
        
        if bond.type != 'single':
            if bond.chiral:
                return False

            if len(bond.atom_1.drawn_neighbours) > 1 and len(bond.atom_2.drawn_neighbours) > 1:
                return False

        chiral = False
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

    @staticmethod
    def can_rotate_around_bond(bond: Bond) -> bool:

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
        overlaps = []
        resolved_atoms = {}
        for atom in self.structure.graph:
            if atom.draw.is_drawn:
                resolved_atoms[atom] = False

        for ring in self.rings:
            for atom in ring.members:
                if resolved_atoms[atom]:
                    continue

                resolved_atoms[atom] = True

                if not atom.adjacent_to_stereobond():

                    non_ring_neighbours = self.get_non_ring_neighbours(atom)

                    if len(non_ring_neighbours) > 1 or (len(non_ring_neighbours) == 1 and len(atom.draw.rings) == 2):
                        overlaps.append({'common': atom,
                                         'rings': atom.draw.rings,
                                         'vertices': non_ring_neighbours})

        for overlap in overlaps:
            branches_to_adjust = overlap['vertices']
            rings = overlap['rings']
            root = overlap['common']

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
                if len(rings) == 2:
                    pass

    def resolve_secondary_overlaps(self, sorted_scores: List[Tuple[float, Atom]]) -> None:
        for score, atom in sorted_scores:
            if score > self.options.overlap_sensitivity:
                if len(atom.drawn_neighbours) <= 1:
                    if atom.drawn_neighbours and atom.drawn_neighbours[0].adjacent_to_stereobond():
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

    def get_atom_nr_to_atom(self) -> None:
        self.atom_nr_to_atom = {}
        for atom in self.structure.graph:
            self.atom_nr_to_atom[atom.nr] = atom

    def get_subtree_overlap_score(self, root: Atom, root_parent: Atom,
                                  atom_to_score: Dict[Atom, float]) -> Tuple[float, Vector]:
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

        sorted_overlaps = []

        for atom in self.drawn_atoms:
            sorted_overlaps.append((overlap_scores[atom], atom))

        sorted_overlaps.sort(key=lambda x: x[0], reverse=True)

        return total, sorted_overlaps, overlap_scores

    @staticmethod
    def get_non_ring_neighbours(atom: Atom) -> List[Atom]:
        non_ring_neighbours = []

        for neighbour in atom.drawn_neighbours:
            nr_overlapping_rings = len(set(atom.draw.rings).intersection(set(neighbour.draw.rings)))
            if nr_overlapping_rings == 0 and not neighbour.draw.is_bridge:
                non_ring_neighbours.append(neighbour)

        return non_ring_neighbours

    def rotate_subtree(self, root: Atom, root_parent: Atom, angle: float, center: Vector):

        for atom in self.traverse_substructure(root, {root_parent}):
            atom.draw.position.rotate_around_vector(angle, center)
            for anchored_ring in atom.draw.anchored_rings:
                if anchored_ring.center:
                    anchored_ring.center.rotate_around_vector(angle, center)

    def rotate_subtree_independent(self, root: Atom, root_parent: Atom, masked_atoms: List[Atom],
                                   angle: float, center: Vector) -> None:

        masked_atoms.append(root_parent)
        masked_atoms = set(masked_atoms)

        for atom in self.traverse_substructure(root, masked_atoms):
            atom.draw.position.rotate_around_vector(angle, center)
            for anchored_ring in atom.draw.anchored_rings:
                if anchored_ring.center:
                    anchored_ring.center.rotate_around_vector(angle, center)

    def traverse_substructure(self, atom: Atom, visited: Set[Atom]) -> Generator[Atom, None, None]:
        yield atom
        visited.add(atom)
        for neighbour in atom.drawn_neighbours:
            if neighbour not in visited:
                yield from self.traverse_substructure(neighbour, visited)

    def get_subgraph_size(self, atom, masked_atoms):
        masked_atoms.add(atom)

        for neighbour in atom.drawn_neighbours:
            if neighbour not in masked_atoms:
                self.get_subgraph_size(neighbour, masked_atoms)

        return len(masked_atoms) - 1

    @staticmethod
    def get_last_atom_with_angle(atom):

        parent_atom = atom.draw.previous_atom
        angle = parent_atom.draw.angle

        while parent_atom and not angle:
            parent_atom = parent_atom.draw.previous_atom
            angle = parent_atom.draw.angle

        return parent_atom

    def create_ring(self, ring, center=None, start_atom=None, previous_atom=None):

        if ring.positioned:
            return

        if center is None:
            center = Vector(0, 0)

        ordered_neighbour_ids = ring.get_ordered_neighbours(self.ring_overlaps)
        starting_angle = 0

        if start_atom:
            starting_angle = Vector.subtract_vectors(start_atom.draw.position, center).angle()

        ring_size = len(ring.members)

        radius = Polygon.find_polygon_radius(self.options.bond_length, ring_size)
        angle = Polygon.get_central_angle(ring_size)

        ring.central_angle = angle

        if start_atom not in ring.members:
            if start_atom:
                start_atom.draw.positioned = False
            start_atom = ring.members[0]

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
        else:
            ring.set_member_positions(self.structure, start_atom, previous_atom, center, starting_angle, radius, angle)

        ring.positioned = True
        ring.center = center

        for neighbour_id in ordered_neighbour_ids:
            neighbour = self.id_to_ring[neighbour_id]
            if neighbour.positioned:
                continue

            atoms = list(RingOverlap.get_vertices(self.ring_overlaps, ring.id, neighbour.id))

            if len(atoms) == 2:

                # This ring is fused
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

                if distance_to_center_2 > distance_to_center_1:
                    next_center = normals[1]

                position_1 = Vector.subtract_vectors(atom_1.draw.position, next_center)
                position_2 = Vector.subtract_vectors(atom_2.draw.position, next_center)

                if position_1.get_clockwise_orientation(position_2) == 'clockwise':
                    if not neighbour.positioned:
                        self.create_ring(neighbour, next_center, atom_1, atom_2)

                else:
                    if not neighbour.positioned:
                        self.create_ring(neighbour, next_center, atom_2, atom_1)

            elif len(atoms) == 1:
                # This ring is a spiro
                ring.spiro = True
                neighbour.spiro = True

                atom = atoms[0]
                next_center = Vector.subtract_vectors(center, atom.draw.position)

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
                self.create_next_bond(neighbour, atom, 0.0)

    @staticmethod
    def set_ring_center(ring):
        total = Vector(0, 0)
        for atom in ring.members:
            total.add(atom.draw.position)

        total.divide(len(ring.members))

        ring.center = total

    def get_current_centre_of_mass(self):
        total = Vector(0, 0)
        count = 0

        for atom in self.structure.graph:
            if atom.draw.positioned:
                total.add(atom.draw.position)
                count += 1

        total.divide(count)

        return total

    @staticmethod
    def atoms_are_in_same_ring(atom_1, atom_2):
        for ring_id_1 in atom_1.draw.rings:
            for ring_id_2 in atom_2.draw.rings:
                if ring_id_1 == ring_id_2:
                    return True
        return False

    def define_rings(self):

        rings = self.structure.cycles.find_sssr()

        if not rings:
            return None

        for ring_members in rings:
            ring = Ring(ring_members)
            self.add_ring(ring)

            for atom in ring_members:
                atom.draw.rings.append(ring.id)

        for i, ring_1 in enumerate(self.rings[:-1]):
            for ring_2 in self.rings[i + 1:]:
                ring_overlap = RingOverlap(ring_1, ring_2)

                if len(ring_overlap.atoms) > 0:
                    self.add_ring_overlap(ring_overlap)

        for ring in self.rings:
            neighbouring_rings = find_neighbouring_rings(self.ring_overlaps, ring.id)
            ring.neighbouring_rings = neighbouring_rings

        for ring in self.rings:
            anchor = ring.members[0]
            if ring not in anchor.draw.anchored_rings:
                anchor.draw.anchored_rings.append(ring)

        self.backup_ring_info()

        while True:
            ring_id = -1

            for ring in self.rings:
                if self.is_part_of_bridged_ring(ring.id) and not ring.bridged:
                    ring_id = ring.id

            if ring_id == -1:
                break

            ring = self.id_to_ring[ring_id]
            involved_ring_ids = []
            self.get_bridged_ring_subrings(ring.id, involved_ring_ids)
            involved_ring_ids = set(involved_ring_ids)

            self.bridged_ring = True
            self.create_bridged_ring(involved_ring_ids)

            for involved_ring_id in involved_ring_ids:
                involved_ring = self.id_to_ring[involved_ring_id]
                self.remove_ring(involved_ring)

        # This is new - if stuff breaks, remove
        bridged_systems = find_bridged_systems(self.rings, self.ring_overlaps)
        if bridged_systems and not self.bridged_ring:
            self.bridged_ring = True
            for bridged_system in bridged_systems:
                involved_ring_ids = set(bridged_system)
                self.create_bridged_ring(involved_ring_ids)
                for involved_ring_id in involved_ring_ids:
                    involved_ring = self.id_to_ring[involved_ring_id]
                    self.remove_ring(involved_ring)

    def hide_hydrogens(self):
        hidden = []
        exposed = []

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

    def get_bridged_rings(self):
        bridged_rings = []
        for ring in self.rings:
            if ring.bridged:
                bridged_rings.append(ring)

        return bridged_rings

    def get_bridged_ring_subrings(self, ring_id, involved_ring_ids):
        involved_ring_ids.append(ring_id)
        ring = self.id_to_ring[ring_id]

        for neighbour_id in ring.neighbouring_rings:
            if neighbour_id not in involved_ring_ids and neighbour_id != ring_id and \
                    rings_connected_by_bridge(self.ring_overlaps, ring_id, neighbour_id):
                self.get_bridged_ring_subrings(neighbour_id, involved_ring_ids)

    def create_bridged_ring(self, involved_ring_ids):
        ring_members = set()
        atoms = set()
        neighbours = set()

        for ring_id in involved_ring_ids:
            ring = self.id_to_ring[ring_id]
            ring.subring_of_bridged = True

            for atom in ring.members:
                atoms.add(atom)

            for neighbour_id in ring.neighbouring_rings:
                neighbours.add(neighbour_id)

        leftovers = set()

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
            bridged_ring.subrings.append(copy.deepcopy(ring))

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
        self.original_rings = copy.deepcopy(self.rings)
        self.original_ring_overlaps = copy.deepcopy(self.ring_overlaps)

        for atom in self.structure.graph:
            atom.draw.original_rings = copy.deepcopy(atom.draw.rings)

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


def draw_multiple(structure: Structure, coords_only: bool = False, options: Union[None, Options] = None) -> Drawer:
    if not options:
        options = Options()
    options_main = Options()
    options_main.finetune = False

    drawer = Drawer(structure, options=options_main, coords_only=True, multiple=True)
    structures = structure.split_disconnected_structures()
    max_x = -100000000

    for i, substructure in enumerate(structures):
        subdrawer = Drawer(substructure, options=options, coords_only=True)
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
        drawer.draw_structure()

    return drawer

