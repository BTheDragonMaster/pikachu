#!/usr/bin/env python

from pikachu.math_functions import Vector
import math

def ring_groups_have_overlap(group_1, group_2, ring_overlaps):
    for ring_1 in group_1:
        for ring_2 in group_2:
            if ring_1 in find_neighbouring_rings(ring_overlaps, ring_2):
                return True
    return False


def get_ring_groups(rings, ring_overlaps):
    ring_groups = []
    for ring in rings:
        ring_groups.append([ring.id])

    current_ring_nr = 0
    previous_ring_nr = -1

    while current_ring_nr != previous_ring_nr:
        previous_ring_nr = current_ring_nr
        indices = None
        new_group = None

        for i, ring_group_1 in enumerate(ring_groups):
            ring_group_1_found = False
            for j, ring_group_2 in enumerate(ring_groups):
                if i != j:
                    if ring_groups_have_overlap(ring_group_1, ring_group_2, ring_overlaps):
                        indices = [i, j]
                        new_group = list(set(ring_group_1 + ring_group_2))
                        ring_group_1_found = True
                        break
            if ring_group_1_found:
                break

        if new_group:
            indices.sort(reverse=True)
            for index in indices:
                ring_groups.pop(index)
            ring_groups.append(new_group)

        current_ring_nr = len(ring_groups)

    return ring_groups


def get_group_overlap_nr(ring_group, ring_overlaps):
    overlaps = 0
    ring_group = set(ring_group)
    for ring_overlap in ring_overlaps:
        if ring_overlap.ring_id_1 in ring_group and ring_overlap.ring_id_2 in ring_group:
            overlaps += 1

    return overlaps


def find_bridged_systems(rings, ring_overlaps):
    bridged_systems = []
    ring_groups = get_ring_groups(rings, ring_overlaps)
    for ring_group in ring_groups:
        ring_nr = len(ring_group)
        overlap_nr = get_group_overlap_nr(ring_group, ring_overlaps)
        if overlap_nr >= ring_nr:
            bridged_systems.append(ring_group)

    return bridged_systems


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


