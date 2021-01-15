#!/usr/bin/env python

import pikachu
from sssr import SSSR

class Ring:
    def __init__(self, members):
        self.id = None
        self.members = members
        self.edges = []
        self.inside_vertices = []
        self.neighbouring_rings = []
        self.positioned = False
        self.center = (0, 0)
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


