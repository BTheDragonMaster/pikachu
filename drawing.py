#!/usr/bin/env python

import copy

from sssr import SSSR
import pikachu
from rings import Ring, RingOverlap, find_neighbouring_rings, rings_connected_by_bridge
from math_functions import Vector, Polygon
import math

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

        self.initialise_matrices()
        self.get_kk_layout()

    def initialise_matrices(self):

        self.distance_matrix = self.get_subgraph_distance_matrix(self.atoms)
        length = len(self.atoms)

        radius = Polygon.find_polygon_radius(500, length)
        angle = Polygon.get_central_angle(length)

        a = 0.0

        self.x_positions = {}
        self.y_positions = {}
        self.positioned = {}

        for atom in self.atoms:
            if not atom.draw.positioned:
                self.x_positions[atom] = center.x + math.cos(a) * radius
                self.y_positions[atom] = center.y + math.sin(a) * radius
            else:
                self.x_positions[atom] = atom.draw.position.x
                self.y_positions[atom] = atom.draw.position.y

            self.positioned[atom] = atom.draw.positioned
            a += angle

        self.length_matrix = {}
        self.spring_strengths = {}
        self.energy_matrix = {}

        self.energy_sums_x = {}
        self.energy_sums_y = {}

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

        while self.max_energy > self.threshold and self.max_iteration > iteration:
            iteration += 1
            max_energy_atom, max_energy, d_ex, d_ey = self.highest_energy()
            delta = max_energy
            inner_iteration = 0

            while delta > self.inner_threshold and self.max_inner_iteration > inner_iteration:
                inner_iteration += 1
                self.update(max_energy_atom, d_ex, d_ey)
                delta, d_ex, d_ey = energy(max_energy_atom)

        for atom in atoms:
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

            denom = 1.0 / (squared_xdiff * squared_ydiff) ** 1.5

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

            self.energy_sums_x += dx - previous_ex
            self.energy_sums_y += dy - previous_ey

        self.energy_sums_x = d_ex
        self.energy_sums_y = d_ey


    def get_subgraph_distance_matrix(self, atoms):
        adjacency_matrix = self.get_subgraph_adjacency_matrix(atoms)
        distance_matrix = {}

        for atom_1 in atoms:
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
    def __init__(self, structure, options=None):
        if options == None:
            self.options = Options()

        self.structure = structure
        self.rings = []
        self.ring_overlaps = []
        self.id_to_ring = {}
        self.bridged_ring = False

        self.ring_id_tracker = 0
        self.ring_overlap_id_tracker = 0

        self.draw()

    def draw(self):
        self.define_rings()
        self.hide_hydrogens()

    def process_structure(self):
        pass

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
            start_atom = list(self.structure.graph.keys())[0]

        self.create_next_bond(start_atom, None, 0.0)

    def create_next_bond(self, atom, previous_atom = None, angle = 0.0,
                         previous_branch_shortest = False, skip_positioning = False):
        if atom.draw.positioned and not skip_positioning:
            return

        double_bond_configuration_set = False

        if not skip_positioning:
            if not previous_atom:
                dummy = Vector(self.options.bond_length, 0)
                dummy.rotate(-60)

                atom.draw.previous_position = dummy
                atom.draw.previous_atom = None
                atom.draw.set_position(Vector(self.options.bond_length, 0))
                atom.draw.angle = -60

                if atom.draw.bridged_ring == None:
                    atom.draw.positioned = True

            elif len(previous_atom.draw.rings) > 0:
                neighbours = previous_atom.neighbours
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
                    position.add(previous_atom.position)

                else:

                    position = joined_vertex.copy()
                    position.rotate_around_vector(180, previous_atom.draw.position)

                atom.draw.set_previous_position(previous_atom)
                atom.draw.set_position(position)
                atom.draw.positioned = True

            else:
                position = Vector(self.options.bond_length, 0)
                position.rotate(angle)
                position.add(previous_atom.position)

                atom.draw.set_position(position)
                atom.draw.set_previous_position(previous_atom)
                atom.draw.positioned = True

        if atom.draw.bridged_ring != None:
            next_ring = self.id_to_ring(atom.draw.bridged_ring)

            if not next_ring.positioned:
                next_center = Vector.subtract_vectors(atom.draw.previous_position, atom.draw.position)
                next_center.invert()
                next_center.normalise()
                scalar = Polygon.find_polygon_radius(self.options.bond_length, len(next_ring.members))

                next_center.multiply_by_scalar(scalar)
                next_center.add(atom.draw.position)

                self.create_ring(next_ring, next_center, atom)

        elif len(atom.draw.rings) > 0:
            next_ring = self.id_to_ring(atom.draw.rings[0])

            if not next_ring.positioned:
                next_center = Vector.subtract_vectors(atom.draw.previous_position, atom.draw.position)
                next_center.invert()
                next_center.normalise()

                radius = Polygon.find_polygon_radius(self.options.bond_length, len(next_ring.members))

                next_center.multiply_by_scalar(radius)
                next_center.add(atom.position)

                self.create_ring(next_ring, next_center, vertex)

        else:
            neighbours = atom.get_drawn_neighbours()

            if previous_atom:
                if previous_atom in neighbours:
                    neighbours.remove(previous_atom)

            previous_angle = atom.get_angle()

            if len(neighbours) == 1:
                next_atom = neighbours[0]

                current_bond = self.structure.bond_lookup[atom][next_atom]
                previous_bond = None

                if previous_atom:
                    previous_bond = self.structure.bond_lookup[previous_atom][atom]


                if current_bond.type == 'triple' or previous_bond.type == 'triple' or \
                        (current_bond.type == 'double' and previous_bond.type == 'double' and
                         previous_atom and len(previous_atom.draw.rings) == 0):
                    #shouldn't this be true?
                    atom.draw.draw_explicit = False

                    if previous_atom:
                        previous_bond.draw.center = True

                    current_bond.draw.center = True

                    if current_bond.type == 'triple' or (previous_atom and previous_bond.type == 'triple'):
                        next_atom.draw.angle = 0.0

                    next_atom.draw.draw_explicit = True

                    self.create_next_bond(next_atom, atom, previous_angle + next_atom.draw.angle)

                elif previous_atom and len(previous_atom.draw.rings) > 0:
                    proposed_angle_1 = 60
                    proposed_angle_2 = -60

                    proposed_vector_1 = Vector(self.options.bond_length, 0)
                    proposed_vector_2 = Vector(self.options.bond_length, 0)

                    proposed_vector_1.rotate(proposed_angle_1).add(atom.draw.position)
                    proposed_vector_2.rotate(proposed_angle_2).add(atom.draw.position)

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

                    if previous_atom and len(previous_atom.neighbours) > 3:
                        if a > 0:
                            a = min([60, a])
                        elif a < 0:
                            a = max([-60, a])
                        else:
                            a = 60

                    elif a == None:
                        last_angled_atom = self.get_last_atom_with_angle(atom)
                        a = last_angled_atom.angle

                        if a == None:
                            #NOT DONE HERE
                            a = 60

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

                    self.create_next_bond(next_atom, atom, previous_angle + next_atom.angle)

            elif len(neighbours) == 2:
                a = atom.draw.angle
                if a == None:
                    a = 60

                neighbour_1, neighbour_2 = neighbours



                subgraph_1_size = self.get_subgraph_size(neighbour_1, atom)
                subgraph_2_size = self.get_subgraph_size(neighbour_2, atom)

                if previous_atom:

                    subgraph_3_size = self.get_subgraph_size(previous_atom, atom)

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

                cis_edge = self.structure.bond_lookup[cis_atom][atom]
                trans_edge = self.structure.bond_lookup[trans_atom][atom]

                previous_branch_shortest = False

                if subgraph_3_size < subgraph_2_size and subgraph_3_size < subgraph_1_size:
                    previous_branch_shortest = True

                trans_atom.draw.angle = a
                cis_atom.draw.angle = -a

                cis_bond = self.structure.bond_lookup[atom][cis_atom]
                trans_bond = self.structure.bond_lookup[atom][trans_atom]

                if cis_bond.type == 'double':
                    if cis_bond.chiral and previous_atom:























    def get_subgraph_size(self, atom, masked_atom):
        seen_atoms = {masked_atom, atom}

        for neighbour in atom.get_drawn_neighbours():
            if neighbour not in seen_atoms:
                self.get_subgraph_size(neighbour, seen_atoms)

        return len(seen_atoms) - 1


    def get_last_atom_with_angle(self, atom):
        angle = None

        parent_atom = atom.draw.previous_atom

        while parent_atom and angle == None:
            parent_atom = atom.draw.previous_atom
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
            starting_angle = Vector.subtract_vectors(start_atom.position, center).angle()

        ring_size = len(ring.members)

        radius = Polygon.find_polygon_radius(self.options.bond_length, ring_size)
        angle = Polygon.get_central_angle(ring_size)

        ring.central_angle = angle

        if start_atom not in ring.members:
            if start_atom:
                start_atom.positioned = False
            start_atom = ring.members[0]

        if ring.bridged:
            KKLayout(self.structure, ring.members, center, start_atom,
                     self.options.bond_length, self.options.kk_threshold, self.options.kk_inner_threshold,
                     self.options.kk_max_inner_iteration, self.options.kk_max_inner_iteration,
                     self.kk_max_energy)
            ring.positioned = True

            self.set_ring_center(ring)
            center = ring.center

            for subring in ring.subrings:
                self.set_ring_center(subring)
        else:
            ring.set_member_positions(self.structure, start_atom, previous_atom)

        ring.positioned = True
        ring.center = center

        for neighbour_id in ordered_neighbour_ids:
            neighbour = self.id_to_ring(neighbour_id)
            if neighbour.positioned:
                continue

            atoms = RingOverlap.get_vertices(self.ring_overlaps, ring.id, neighbour.id)

            if len(atoms) == 2:
                #This ring is fused
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

                normals[0].multiply_by_scalar(apothem).add(midpoint)
                normals[1].multiply_by_scalar(apothem).add(midpoint)

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
            for neighbour in atom.neighbours:
                if neighbour.positioned:
                    continue

                atom.draw.connected_to_ring = True
                self.create_next_bond(neighbour, atom, 0.0)


    def set_ring_center(self, ring):
        total = Vector(0, 0)
        for atom in ring.members:
            total.add(atom.draw.position)

        ring.center = total.divide(len(ring.members))

    def get_current_centre_of_mass(self):
        total = Vector(0, 0)
        count = 0

        for atom in self.structure.graph:
            if atom.draw.positioned:
                total.add(atom.draw.position)
                count += 1

        return total.divide(count)











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
            anchor.draw.ring_anchors.add(ring.id)

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
        for atom in self.structure.graph:
            if atom.type != 'H':
                continue

            neighbour = atom.neighbours[0]
            neighbour.draw.has_hydrogen = True
            if not neighbour.chiral or (len(neighbour.draw.rings) < 2 and neighbour.draw.bridged_ring == None) or \
                    (neighbour.draw.bridged_ring != None and len(neighbour.draw.original_rings) < 2):
                atom.draw.is_drawn = False
                hidden.append(atom)

            else:
                exposed.append(atom)

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
            intersect = involved_ring_ids & atom.draw.rings

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
            atom.draw.original_rings = copy.copy(atom.draw.rings)

    def get_ring_index(self, ring_id):
        for i, ring in enumerate(self.rings):
            if ring.id == ring_id:
                return i

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

class Options:
    def __init__(self):
        self.width = 500
        self.height = 500
        self.bond_thickness = 0.6
        self.bond_length = 15
        self.short_bond_length = 0.85
        self.bond_spacing = 0.18 * bond_length
        self.isomeric = True
        self.padding = 20
        self.font_size_large = 5
        self.font_size_small = 3
        self.kk_threshold = 0.1
        self.kk_inner_threshold = 0.1
        self.kk_max_iteration = 2000
        self.kk_max_inner_iteration = 50
        self.kk_max_energy = 1e9


if __name__ == "__main__":
    daptomycin_smiles = pikachu.Smiles('CCCCCCCCCC(=O)N[C@@H](CC1=CNC2=CC=CC=C21)C(=O)N[C@@H](CC(=O)N)C(=O)N[C@@H](CC(=O)O)C(=O)N[C@H]3[C@H](OC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@H](NC(=O)CNC(=O)[C@@H](NC(=O)[C@H](NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)CNC3=O)CCCN)CC(=O)O)C)CC(=O)O)CO)[C@H](C)CC(=O)O)CC(=O)C4=CC=CC=C4N)C')
    daptomycin_structure = daptomycin_smiles.smiles_to_structure()
    teicoplanin_smiles = pikachu.Smiles('CCCCCCCCCC(=O)N[C@@H]1[C@H]([C@@H]([C@H](O[C@H]1OC2=C3C=C4C=C2OC5=C(C=C(C=C5)[C@H]([C@H]6C(=O)N[C@@H](C7=C(C(=CC(=C7)O)OC8[C@H]([C@H]([C@@H]([C@H](O8)CO)O)O)O)C9=C(C=CC(=C9)[C@H](C(=O)N6)NC(=O)[C@@H]4NC(=O)[C@@H]1C2=CC(=CC(=C2)OC2=C(C=CC(=C2)[C@H](C(=O)N[C@H](CC2=CC(=C(O3)C=C2)Cl)C(=O)N1)N)O)O)O)C(=O)O)O[C@H]1[C@@H]([C@H]([C@@H]([C@H](O1)CO)O)O)NC(=O)C)Cl)CO)O)O')
    teicoplanin_structure = teicoplanin_smiles.smiles_to_structure()
    test_smiles = pikachu.Smiles('C1CCCC2CC1CCCC2')
    test_structure = test_smiles.smiles_to_structure()
    graph = SSSR(teicoplanin_structure)
    rings = graph.get_rings()

    drawer = Drawer(daptomycin_structure)
    print(drawer.rings)