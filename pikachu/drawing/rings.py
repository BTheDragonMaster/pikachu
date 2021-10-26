#!/usr/bin/env python

from pikachu.math_functions import Vector
import math

class Ring:
    def __init__(self, members):
        self.id = None
        self.members = members
        self.edges = []
        self.inside_vertices = []
        self.neighbouring_rings = []
        self.positioned = False
        self.center = Vector(0, 0)
        self.subrings = []
        self.bridged = False
        self.subring_of_bridged = False
        self.spiro = False
        self.fused = False
        self.central_angle = 0.0
        self.flippable = True

    def __hash__(self):
        return self.id

    def __eq__(self, other):
        return self.id == other.id

    def __repr__(self):
        return str(self.id) + ' ' + '-'.join([atom.__repr__() for atom in self.members])

    def get_angle(self):
        return math.pi - self.central_angle

    def get_ordered_neighbours(self, ring_overlaps):
        ordered_neighbours_and_atom_nrs = []

        for neighbour_id in self.neighbouring_rings:
            atoms = RingOverlap.get_vertices(ring_overlaps, self.id, neighbour_id)
            ordered_neighbours_and_atom_nrs.append((len(atoms), neighbour_id))

        ordered_neighbours_and_atom_nrs = sorted(ordered_neighbours_and_atom_nrs, key=lambda x: x[0], reverse=True)
        ordered_neighbour_ids = [x[1] for x in ordered_neighbours_and_atom_nrs]

        return ordered_neighbour_ids

    def set_member_positions(self, structure, start_atom, previous_atom, center, a, radius, angle):
        current_atom = start_atom
        iteration = 0

        while current_atom != None and iteration < 100:
            previous = current_atom

            if not previous.draw.positioned:
                x = center.x + math.cos(a) * radius
                y = center.y + math.sin(a) * radius
                previous.draw.set_position(Vector(x, y))

            a += angle

            if not self.bridged or len(self.subrings) < 3:
                previous.draw.angle = a
                previous.draw.positioned = True

            current_atom = structure.get_next_in_ring(self, current_atom, previous_atom)
            previous_atom = previous

            if current_atom == start_atom:
                current_atom = None

            iteration += 1


class RingOverlap:
    def __init__(self, ring_1, ring_2):
        self.id = None
        self.ring_id_1 = ring_1.id
        self.ring_id_2 = ring_2.id
        self.atoms = set()

        for atom_1 in ring_1.members:
            for atom_2 in ring_2.members:
                if atom_1 == atom_2:
                    self.atoms.add(atom_1)

    def __hash__(self):
        return self.id

    def __eq__(self, other):
        return self.id == other.id

    def involves_ring(self, ring_id):
        if self.ring_id_1 == ring_id or self.ring_id_2 == ring_id:
            return True
        else:
            return False

    def update_other(self, ring_id, other_ring_id):
        if self.ring_id_1 == other_ring_id:
            self.ring_id_2 = ring_id
        else:
            self.ring_id_1 = ring_id

    def is_bridge(self):
        if len(self.atoms) > 2:
            return True

        for atom in self.atoms:
            if len(atom.draw.rings) > 2:
                return True

        return False


    @staticmethod
    def get_vertices(ring_overlaps, ring_id_1, ring_id_2):
        for ring_overlap in ring_overlaps:
            if (ring_overlap.ring_id_1 == ring_id_1 and ring_overlap.ring_id_2 == ring_id_2) or\
                    (ring_overlap.ring_id_1 == ring_id_2 and ring_overlap.ring_id_2 == ring_id_1):
                return ring_overlap.atoms




def find_neighbouring_rings(ring_overlaps, ring_id):
    neighbouring_rings = []

    for ring_overlap in ring_overlaps:
        if ring_overlap.ring_id_1 == ring_id:
            neighbouring_rings.append(ring_overlap.ring_id_2)
        elif ring_overlap.ring_id_2 == ring_id:
            neighbouring_rings.append(ring_overlap.ring_id_1)

    return neighbouring_rings

def rings_connected_by_bridge(ring_overlaps, ring_id_1, ring_id_2):
    for ring_overlap in ring_overlaps:
        if ring_id_1 == ring_overlap.ring_id_1 and ring_id_2 == ring_overlap.ring_id_2:
            return ring_overlap.is_bridge()
        if ring_id_2 == ring_overlap.ring_id_1 and ring_id_1 == ring_overlap.ring_id_2:
            return ring_overlap.is_bridge()

    return False


