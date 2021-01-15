#!/usr/bin/env python

import copy

from sssr import SSSR
import pikachu
from rings import Ring, RingOverlap, find_neighbouring_rings, rings_connected_by_bridge
from math_functions import Vector, add_vectors, subtract_vectors

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
                         origin_is_shortest = False, skip_positioning = False):
        if atom.draw.positioned and not skip_positioning:
            return

        double_bond_configuration_set = False

        if not skip_positioning:
            if not previous_atom:
                dummy = Vector(self.options.bond_length, 0)
                dummy.rotate(-60)

                atom.draw.previous_position = dummy
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

                atom.draw.previous_position = previous_atom.position
                atom.draw.set_position(position)
                atom.draw.positioned = True

            else:
                position = Vector(self.options.bond_length, 0)
                position.rotate(angle)
                position.add(previous_atom.position)

                atom.set_position(position)
                atom.previous_position = previous_atom.position
                atom.positioned = True

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

    def create_ring(self, ring, center = None, start_atom = None, previous_atom = None):
        if ring.positioned:
            return

        if center == None:
            center = Vector(0, 0)

        ordered_neighbours = ring.get_ordered_neighbours(self.ring_overlaps)




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
            if not neighbour.chiral or (len(neighbour.draw.rings) < 2 and neighbour.draw.bridged_ring == None) or\
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