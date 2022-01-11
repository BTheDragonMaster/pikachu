#!/usr/bin/env python
import time
import copy
import math
import matplotlib
# matplotlib.use('TkAgg')
from matplotlib import pyplot as plt

from pikachu.drawing.sssr import SSSR
from pikachu.drawing.rings import Ring, RingOverlap, find_neighbouring_rings, rings_connected_by_bridge
from pikachu.math_functions import Vector, Polygon, Line, Permutations
from pikachu.chem.chirality import find_chirality_from_nonh, get_chiral_permutations
from pikachu.chem.atom_properties import ATOM_PROPERTIES


class SpanningTree:
    def __init__(self, structure):
        self.spanning_tree = {}
        atom = list(structure.graph.keys())[0]
        self.visited = set()

        self.make_spanning_tree(atom)

    def make_spanning_tree(self, atom, parent=None):

        if parent:
            if parent not in self.spanning_tree:
                self.spanning_tree[parent] = []

            self.spanning_tree[parent].append(atom)

        self.visited.add(atom)
        for neighbour in atom.neighbours:
            if neighbour not in self.visited:
                self.make_spanning_tree(neighbour, atom)


class Vertex:
    def __init__(self, atom, parent):
        self.atom = atom
        self.parent = parent
        self.children = []
        if parent:
            self.parent.children.append(self)

    def __hash__(self):
        return self.atom.__hash__()

    def __repr__(self):
        return self.atom.__repr__()

    def __eq__(self, vertex):
        if self.atom == vertex.atom:
            return True
        else:
            return False


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

            if delta > max_energy and self.positioned[atom] == False:
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
            if not atom_1 in distance_matrix:
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
    def __init__(self, structure, options=None, save_png=None):
        if options == None:
            self.options = Options()

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

        self.ring_id_tracker = 0
        self.ring_overlap_id_tracker = 0

        self.draw()

    def prioritise_chiral_bonds(self, chiral_center):

        subtree_1_size = self.get_subgraph_size(chiral_center.neighbours[0], {chiral_center})
        subtree_2_size = self.get_subgraph_size(chiral_center.neighbours[1], {chiral_center})
        subtree_3_size = self.get_subgraph_size(chiral_center.neighbours[2], {chiral_center})
        subtree_4_size = self.get_subgraph_size(chiral_center.neighbours[3], {chiral_center})

        sizes_and_atoms = [(subtree_1_size, chiral_center.neighbours[0]),
                           (subtree_2_size, chiral_center.neighbours[1]),
                           (subtree_3_size, chiral_center.neighbours[2]),
                           (subtree_4_size, chiral_center.neighbours[3])]

        sizes_and_atoms.sort(key=lambda x: (x[0], ATOM_PROPERTIES.element_to_atomic_nr[x[1].type]))

        options_h = []
        options = []
        backup_options_rings = []
        backup_options_chiral = []

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

                if not other_chiral_centre and not in_ring and not neighbour.chiral:
                    options.append(neighbour)
                elif in_ring and not other_chiral_centre and not neighbour.chiral:
                    backup_options_rings.append(neighbour)
                else:
                    backup_options_chiral.append(neighbour)

        priority = options_h + options + backup_options_rings + backup_options_chiral
        assert len(priority) == 4

        return priority

    def determine_chirality(self, chiral_center):

        priority = self.prioritise_chiral_bonds(chiral_center)

        priority_matches_chirality = False

        if tuple(priority) in get_chiral_permutations(chiral_center.neighbours):
            priority_matches_chirality = True

        # First item in priority will always be a wedge in the opposite direction to the
        # second wedge, or won't be drawn.

        position_1 = priority[1].draw.position
        position_2 = priority[2].draw.position
        position_3 = priority[3].draw.position

        direction = Vector.get_directionality_triangle(position_1,
                                                       position_2,
                                                       position_3)

        chirality = chiral_center.chiral

        if chirality == 'clockwise':
            if direction == 'clockwise' and not priority_matches_chirality:
                wedge_1 = 'front'
            elif direction == 'clockwise' and priority_matches_chirality:
                wedge_1 = 'back'
            elif direction == 'counterclockwise' and not priority_matches_chirality:
                wedge_1 = 'back'
            else:
                wedge_1 = 'front'

        else:
            if direction == 'clockwise' and not priority_matches_chirality:
                wedge_1 = 'back'
            elif direction == 'clockwise' and priority_matches_chirality:
                wedge_1 = 'front'
            elif direction == 'counterclockwise' and not priority_matches_chirality:
                wedge_1 = 'front'
            else:
                wedge_1 = 'back'

        bonds_and_wedges = []

        bond_1 = self.structure.bond_lookup[chiral_center][priority[1]]
        bonds_and_wedges.append((bond_1, wedge_1))

        if priority[0].draw.is_drawn:

            bond_2 = self.structure.bond_lookup[chiral_center][priority[0]]

            if wedge_1 == 'front':
                wedge_2 = 'back'
            else:
                wedge_2 = 'front'

            bonds_and_wedges.append((bond_2, wedge_2))

        return bonds_and_wedges

    def set_chiral_bonds(self):
        for atom in self.structure.graph:
            if atom.chiral:
                bonds_and_wedges = self.determine_chirality(atom)
                for bond, wedge in bonds_and_wedges:
                    self.chiral_bonds.append(bond)
                    self.chiral_bond_to_orientation[bond] = (wedge, atom)

    def draw(self):
        start_time = time.time()

        if not self.options.draw_hydrogens:
            self.hide_hydrogens()
      #  print("Hiding hydrogens..")
        time_1 = time.time()
      #  print(time_1 - start_time)
        self.get_atom_nr_to_atom()
      #  print("Making atom dictionary..")
        time_2 = time.time()
      #  print(time_2 - time_1)
        self.define_rings()
     #   print("Defining rings..")
        time_3 = time.time()
      #  print(time_3 - time_2)
        self.process_structure()
     #   print("Processing structure..")
        time_4 = time.time()
     #   print(time_4 - time_3)
        self.set_chiral_bonds()
     #   print("Setting chiral bonds..")
        time_5 = time.time()
      #  print(time_5 - time_4)
        self.draw_structure()
      #  print("Drawing structure..")
        time_6 = time.time()
     #   print(time_6 - time_5)
        #self.draw_svg()
      #  self.draw_png()


    def get_hydrogen_text_orientation(self, atom):
        try:
            neighbour = atom.drawn_neighbours[0]
            if neighbour.draw.position.x > atom.draw.position.x + 3:
                orientation = 'H_before_atom'
            else:
                orientation = 'H_after_atom'

            return orientation
        except IndexError:
            return 'H_before_atom'

    @staticmethod
    def in_same_ring(atom_1, atom_2):
        if atom_1.draw.rings and atom_2.draw.rings:
            joined_rings = list(set(atom_1.draw.rings).intersection(set(atom_2.draw.rings)))
            if joined_rings:
                return True

        return False

    @staticmethod
    def get_common_rings(atom_1, atom_2):
        if atom_1.draw.rings and atom_2.draw.rings:
            joined_rings = list(set(atom_1.draw.rings).intersection(set(atom_2.draw.rings)))
            return joined_rings

        return None

    def plot_chiral_bond_front(self, polygon, ax, color='black'):
        x = []
        y = []
        for point in polygon:
            x.append(point.x)
            y.append(point.y)

        ax.fill(x, y, color=color)

    def plot_chiral_bond_back(self, lines, ax, color='black'):
        for line in lines:
            self.plot_line(line, ax, color=color)

    def plot_chiral_bond(self, orientation, chiral_center, line, ax, midpoint):
        halflines = line.divide_in_two(midpoint)
        for halfline in halflines:
            truncated_line = halfline.get_truncated_line(self.options.short_bond_length)
            if orientation == 'front':
                bond_wedge = truncated_line.get_bond_wedge_front(self.options.chiral_bond_width, chiral_center)
                self.plot_chiral_bond_front(bond_wedge, ax, color=halfline.atom.draw.colour)
            else:
                bond_lines = halfline.get_bond_wedge_back(self.options.chiral_bond_width, chiral_center)
                self.plot_chiral_bond_back(bond_lines, ax, color=halfline.atom.draw.colour)

    def plot_halflines(self, line, ax, midpoint):
        halflines = line.divide_in_two(midpoint)
        for halfline in halflines:
            truncated_line = halfline.get_truncated_line(self.options.short_bond_length)
            self.plot_line(truncated_line, ax, color=halfline.atom.draw.colour)

    def plot_halflines_double(self, line, ax, midpoint):
        halflines = line.divide_in_two(midpoint)
        for halfline in halflines:
            self.plot_line(halfline, ax, color=halfline.atom.draw.colour)

    def plot_line(self, line, ax, color='black'):
        ax.plot([line.point_1.x, line.point_2.x],
                 [line.point_1.y, line.point_2.y], color=color, linewidth=self.line_width)

    def save_svg(self, out_file):
        if out_file.endswith('.svg'):
            pass
        else:
            out_file += '.svg'
        plt.savefig(out_file)
        plt.clf()
        plt.close(plt.gcf())
        plt.close('all')

    def save_png(self, out_file):
        if out_file.endswith('.png'):
            pass
        else:
            out_file += '.png'
        plt.savefig(out_file)
        plt.clf()
        plt.close()

    def show_molecule(self):
        plt.show()
        plt.clf()
        plt.close()

    def draw_structure(self):

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

     #   print("Height:", height)
      #  print("Width:", width)

       # font_size = 3500 / height
        self.line_width = 2

        fig, ax = plt.subplots(figsize=((width + 2 * self.options.padding) / 50.0, (height + 2 * self.options.padding) / 50.0), dpi=100)
      #  fig, ax = plt.subplots()
        ax.set_aspect('equal', adjustable='box')
        ax.axis('off')

        ax.set_xlim([min_x - self.options.padding, max_x + self.options.padding])
        ax.set_ylim([min_y - self.options.padding, max_y + self.options.padding])
        plt.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)

      #  figure, ax = plt.subplots(figsize=(8, 8))
       # ax.patch.set_face_color(self.options.background_color)
      #  ax.set_aspect()('equal', adjustable='box')
      #  plt.gca().set_aspect('equal', adjustable='box')

        params = {'mathtext.default': 'regular',}
       #           'font.size': font_size}
        plt.rcParams.update(params)

        ring_centers_x = []
        ring_centers_y = []

        for ring in self.rings:
            self.set_ring_center(ring)

            ring_centers_x.append(ring.center.x)
            ring_centers_y.append(ring.center.y)

     #   ax.scatter(ring_centers_x, ring_centers_y, color='blue')

        for bond_nr, bond in self.structure.bonds.items():
            if bond.atom_1.draw.positioned and bond.atom_2.draw.positioned:
                line = Line(bond.atom_1.draw.position, bond.atom_2.draw.position, bond.atom_1, bond.atom_2)
                midpoint = line.get_midpoint()
                truncated_line = line.get_truncated_line(self.options.short_bond_length)
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
                if atom.type != 'C':
                    text = atom.type
                else:
                    text = ''

                horizontal_alignment = 'center'

                orientation = self.get_hydrogen_text_orientation(atom)

                if not atom.charge and atom.type != 'C':

                    if atom.draw.has_hydrogen:
                   # if len(atom.drawn_neighbours) == 1 and atom.draw.has_hydrogen:
                        hydrogen_count = 0
                        for neighbour in atom.neighbours:
                            if neighbour.type == 'H' and not neighbour.draw.is_drawn:
                                hydrogen_count += 1

                        if hydrogen_count:

                            if hydrogen_count > 1:
                                if orientation == 'H_before_atom':
                                    text = r'$H_{hydrogens}{atom_type}$'.format(hydrogens=hydrogen_count,
                                                                                atom_type=atom.type)
                                    horizontal_alignment = 'right'
                                    atom.draw.position.x += 3
                                else:
                                    text = r'${atom_type}H_{hydrogens}$'.format(hydrogens=hydrogen_count,
                                                                                atom_type=atom.type)
                                    horizontal_alignment = 'left'
                                    atom.draw.position.x -= 3
                            elif hydrogen_count == 1:
                                if orientation == 'H_before_atom':
                                    text = f'H{atom.type}'
                                    horizontal_alignment = 'right'
                                    atom.draw.position.x += 3
                                else:
                                    text = f'{atom.type}H'
                                    horizontal_alignment = 'left'
                                    atom.draw.position.x -= 3

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

                            text = r'${atom_type}^{charge}{charge_symbol}$'.format(charge=atom.charge,
                                                                                   atom_type=atom.type,
                                                                                   charge_symbol=charge_symbol)
                        elif abs(atom.charge) == 1:
                            text = r'${atom_type}^{charge_symbol}$'.format(atom_type=atom.type,
                                                                           charge_symbol=charge_symbol)

                        horizontal_alignment = 'left'
                        atom.draw.position.x -= 3
                    else:

                        if hydrogen_count > 1:
                            if orientation == 'H_before_atom':
                                if abs(atom.charge) > 1:
                                    text = r'$H_{hydrogens}{atom_type}^{charge}{charge_symbol}$'.format(hydrogens=hydrogen_count,
                                                                                                        atom_type=atom.type,
                                                                                                        charge=atom.charge,
                                                                                                        charge_symbol=charge_symbol)
                                elif abs(atom.charge) == 1:
                                    text = r'$H_{hydrogens}{atom_type}^{charge_symbol}$'.format(hydrogens=hydrogen_count,
                                                                                                atom_type=atom.type,
                                                                                                charge_symbol=charge_symbol)

                                horizontal_alignment = 'right'
                                atom.draw.position.x += 3

                            else:
                                if abs(atom.charge) > 1:
                                    text = r'${atom_type}H_{hydrogens}^{charge}{charge_symbol}$'.format(hydrogens=hydrogen_count,
                                                                                                        atom_type=atom.type,
                                                                                                        charge=atom.charge,
                                                                                                        charge_symbol=charge_symbol)
                                elif abs(atom.charge) == 1:
                                    text = r'${atom_type}H_{hydrogens}^{charge_symbol}$'.format(hydrogens=hydrogen_count,
                                                                                                atom_type=atom.type,
                                                                                                charge_symbol=charge_symbol)

                                horizontal_alignment = 'left'
                                atom.draw.position.x -= 3
                        elif hydrogen_count == 1:
                            if orientation == 'H_before_atom':
                                if abs(atom.charge) > 1:

                                    text = r'$H{atom_type}^{charge}{charge_symbol}$'.format(atom_type=atom.type,
                                                                                            charge=atom.charge,
                                                                                            charge_symbol=charge_symbol)
                                elif abs(atom.charge) == 1:
                                    text = r'$H{atom_type}^{charge_symbol}$'.format(atom_type=atom.type,
                                                                                    charge_symbol=charge_symbol)
                                horizontal_alignment = 'right'
                                atom.draw.position.x += 3
                            else:
                                if abs(atom.charge) > 1:
                                    text = r'${atom_type}H^{charge}{charge_symbol}$'.format(atom_type=atom.type,
                                                                                            charge=atom.charge,
                                                                                            charge_symbol=charge_symbol)

                                elif abs(atom.charge) == 1:
                                    text = r'${atom_type}H^{charge_symbol}$'.format(atom_type=atom.type,
                                                                                    charge_symbol=charge_symbol)
                                horizontal_alignment = 'left'
                                atom.draw.position.x -= 3

                if text:
                    plt.text(atom.draw.position.x, atom.draw.position.y,
                             text,
                             horizontalalignment=horizontal_alignment,
                             verticalalignment='center',
                             color=atom.draw.colour)

    def is_terminal(self, atom):
        if len(atom.drawn_neighbours) <= 1:
            return True

        return False

    def process_structure(self):
        a_time = time.time()
        self.position()
        pos_time_1 = time.time()
       # print("Positioning...")
       # print(pos_time_1 - a_time)
        self.structure.refresh_structure()
        pos_time_2 = time.time()
      #  print("Refreshing...")
      #  print(pos_time_2 - pos_time_1)
        self.restore_ring_information()

        self.resolve_primary_overlaps()

        self.total_overlap_score, sorted_overlap_scores, atom_to_scores = self.get_overlap_score()

        for i in range(self.options.overlap_resolution_iterations):
            for bond in self.drawn_bonds:
                if self.bond_is_rotatable(bond):

                    tree_depth_1 = self.get_subgraph_size(bond.atom_1, {bond.atom_2})
                    tree_depth_2 = self.get_subgraph_size(bond.atom_2, {bond.atom_1})

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
                                    angle_1 = neighbour_1.draw.position.get_rotation_away_from_vector(atom_1.position, atom_2.position, math.radians(120))
                                    angle_2 = neighbour_2.draw.position.get_rotation_away_from_vector(atom_1.position, atom_2.position, math.radians(120))

                                    self.rotate_subtree(neighbour_1, atom_2, angle_1, atom_2.position)
                                    self.rotate_subtree(neighbour_2, atom_2, angle_2, atom_2.position)

                                    new_overlap_score, _, _ = self.get_overlap_score()

                                    if new_overlap_score > self.total_overlap_score:
                                        self.rotate_subtree(neighbour_1, atom_2, -angle_1, atom_2.position)
                                        self.rotate_subtree(neighbour_2, atom_2, -angle_2, atom_2.position)
                                    else:
                                        self.total_overlap_score = new_overlap_score

                        self.total_overlap_score, sorted_overlap_scores, atom_to_scores = self.get_overlap_score()

        self.resolve_secondary_overlaps(sorted_overlap_scores)

        pos_time_3 = time.time()
      #  print("Overlap resolution...")
      #  print(pos_time_3 - pos_time_2)

    def position(self):
        start_atom = None

        for atom in self.structure.graph:
            if atom.draw.bridged_ring != None:
                start_atom = atom
                break

        #is this necessary?

        for ring in self.rings:
            if ring.bridged:
                start_atom = ring.members[0]

        if len(self.rings) > 0 and start_atom == None:
            start_atom = self.rings[0].members[0]

        if start_atom == None:
            start_atom = list(self.drawn_atoms)[0]

        self.create_next_bond(start_atom, None, 0.0)

    def create_next_bond(self, atom, previous_atom=None, angle=0.0,
                         previous_branch_shortest=False, skip_positioning=False):

     #   print(atom)

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

                if atom.draw.bridged_ring == None:
                    atom.draw.positioned = True

            elif len(previous_atom.draw.rings) > 0:
                neighbours = previous_atom.drawn_neighbours
                joined_vertex = None
                position = Vector(0, 0)

                if previous_atom.draw.bridged_ring == None and len(previous_atom.draw.rings) > 1:
                    for neighbour in neighbours:
                        if len(set(neighbour.draw.rings) & set(previous_atom.draw.rings)) == len(previous_atom.draw.rings):
                            joined_vertex = neighbour
                            break

                if not joined_vertex:
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
        if atom.draw.bridged_ring != None:
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

                    atom.draw.draw_explicit = False

                    if previous_atom:
                        previous_bond.draw.center = True

                    current_bond.draw.center = True

                    if current_bond.type == 'double' or current_bond.type == 'triple' or (previous_atom and previous_bond.type == 'triple'):
                        next_atom.draw.angle = 0.0

                    #TODO: if bond type is double, make sure that the carbon in the middle is drawn

                    next_atom.draw.draw_explicit = True

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
                        if a > 0:
                            a = min([math.radians(60), a])
                        elif a < 0:
                            a = max([-math.radians(60), a])
                        else:
                            a = math.radians(60)

                    elif not a:
                        last_angled_atom = self.get_last_atom_with_angle(atom)
                        a = last_angled_atom.draw.angle

                        if not a:
                            a = math.radians(60)

                    if previous_atom:

                        bond = self.structure.bond_lookup[previous_atom][atom]
                        if bond.type == 'double' and bond.chiral:

                            previous_previous_atom = previous_atom.draw.previous_atom
                            if previous_previous_atom:

                                configuration = bond.chiral_dict[previous_previous_atom][next_atom]
                                if configuration == 'cis':

                                    a = -a

                    if previous_branch_shortest:
                        next_atom.draw.angle = a
                    else:
                        next_atom.draw.angle = -a

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

    def restore_ring_information(self):
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

    def bond_is_rotatable(self, bond):
        if bond.type != 'single':
            return False

        #todo: If bond type IS single, make sure you can't rotate it if the adjacent bond is chiral

        # If bond is terminal, don't bother rotating.

        if len(bond.atom_1.drawn_neighbours) == 1 or len(bond.atom_2.drawn_neighbours) == 1:
            return False

        if bond.atom_1.draw.rings and \
                bond.atom_2.draw.rings and \
                len(set(bond.atom_1.draw.rings).intersection(set(bond.atom_2.draw.rings))) > 0:
            return False

        return True

    def resolve_primary_overlaps(self):
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

                non_ring_neighbours = self.get_non_ring_neighbours(atom)

                if len(non_ring_neighbours) > 1 or\
                    (len(non_ring_neighbours) == 1 and len(atom.draw.rings) == 2):
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

    def resolve_secondary_overlaps(self, sorted_scores):
        for score, atom in sorted_scores:
            if score > self.options.overlap_sensitivity:
                if len(atom.drawn_neighbours) <= 1:
                    closest_atom = self.get_closest_atom(atom)

                    if closest_atom:
                        closest_position = None

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

                    atom.draw.position.rotate_away_from_vector(closest_position, atom_previous_position, math.radians(20))


    def get_atom_nr_to_atom(self):
        self.atom_nr_to_atom = {}
        for atom in self.structure.graph:
            self.atom_nr_to_atom[atom.nr] = atom

    def get_subtree_overlap_score(self, root, root_parent, atom_to_score):
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

    def get_overlap_score(self):
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

    def get_non_ring_neighbours(self, atom):
        non_ring_neighbours = []

        for neighbour in atom.drawn_neighbours:
            nr_overlapping_rings = len(set(atom.draw.rings).intersection(set(neighbour.draw.rings)))
            if nr_overlapping_rings == 0 and not neighbour.draw.is_bridge:
                non_ring_neighbours.append(neighbour)

        return non_ring_neighbours

    def rotate_subtree(self, root, root_parent, angle, center):

        for atom in self.traverse_substructure(root, {root_parent}):
            atom.draw.position.rotate_around_vector(angle, center)
            for anchored_ring in atom.draw.anchored_rings:
                if anchored_ring.center:
                    anchored_ring.center.rotate_around_vector(angle, center)

    def rotate_subtree_independent(self, root, root_parent, masked_atoms, angle, center):

        masked_atoms.append(root_parent)
        masked_atoms = set(masked_atoms)

        for atom in self.traverse_substructure(root, {masked_atoms}):
            atom.draw.position.rotate_around_vector(angle, center)
            for anchored_ring in atom.draw.anchored_rings:
                if anchored_ring.center:
                    anchored_ring.center.rotate_around_vector(angle, center)

    def traverse_substructure(self, atom, visited):
        yield atom
        visited.add(atom)
        for neighbour in atom.drawn_neighbours:
           # if neighbour in self.structure.get_drawn_atoms():
            if neighbour not in visited:
                yield from self.traverse_substructure(neighbour, visited)

    def get_subgraph_size(self, atom, masked_atoms):
        masked_atoms.add(atom)

        for neighbour in atom.drawn_neighbours:
            if neighbour not in masked_atoms:
                self.get_subgraph_size(neighbour, masked_atoms)

        return len(masked_atoms) - 1

    def get_last_atom_with_angle(self, atom):

        parent_atom = atom.draw.previous_atom
        angle = parent_atom.draw.angle

        while parent_atom and not angle:
            parent_atom = parent_atom.draw.previous_atom
            angle = parent_atom.draw.angle

        return parent_atom

    def create_ring(self, ring, center = None, start_atom = None, previous_atom = None):

        if ring.positioned:
            return

        if center == None:
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
                #This ring is a spiro
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

    def set_ring_center(self, ring):
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

    def atoms_are_in_same_ring(self, atom_1, atom_2):
        for ring_id_1 in atom_1.draw.rings:
            for ring_id_2 in atom_2.draw.rings:
                if ring_id_1 == ring_id_2:
                    return True
        return False

    def define_rings(self):
        rings = SSSR(self.structure).get_rings()
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

    def hide_hydrogens(self):
        hidden = []
        exposed = []
        sometime = time.time()
        self.structure.refresh_structure()
      #  print('refreshing structure')
        newtime = time.time()
      #  print(newtime - sometime)
        for atom in self.structure.graph:
            if atom.type != 'H':
                continue

            neighbour = atom.neighbours[0]
            neighbour.draw.has_hydrogen = True
            atom.draw.is_drawn = False
            hidden.append(atom)

            if len(neighbour.draw.rings) < 2 and neighbour.draw.bridged_ring == None and \
                    neighbour.draw.bridged_ring != None and len(neighbour.draw.original_rings) < 2:

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
            for id in ring_ids:
                if (ring_overlap.ring_id_1 == ring_id and ring_overlap.ring_id_2 == id) or\
                        (ring_overlap.ring_id_2 == ring_id and ring_overlap.ring_id_1 == id):
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

    def get_sorted_distances_from_list(self, atom, atom_list):
        atom_distances = []

        for atom_2 in atom_list:
            if atom == atom_2:
                continue

            squared_distance = atom.draw.position.get_squared_distance(atom_2.draw.position)
            atom_distances.append((squared_distance, atom_2))

        atom_distances.sort(key=lambda x: x[0])

        return atom_distances




class Options:
    def __init__(self):
        self.width = 500
        self.height = 500
        self.bond_thickness = 0.6
        self.bond_length = 15
        self.chiral_bond_width = self.bond_length * 0.1
        self.bond_length_squared = self.bond_length ** 2
        self.short_bond_length = 0.50
        self.double_bond_length = 0.80
        self.bond_spacing = 0.18 * self.bond_length
        self.isomeric = True
        self.padding = 20
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
