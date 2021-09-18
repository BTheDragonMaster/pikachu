from pprint import pprint
import copy

from pikachu.fingerprinting.daylight import Daylight
from pikachu.fingerprinting.hashing import hash_32_bit_integer
from pikachu.chem.bond_properties import BOND_PROPERTIES
from pikachu.chem.chirality import find_chirality_from_nonh


class ECFP:
    def __init__(self, structure, iterations=2):
        self.structure = structure
        self.iterations = iterations

        self.identifiers = {}
        self.bonds = {}
        self.fingerprint = set()
        self.seen_atoms = {}
        self.features = {}
        self.set_initial_identifiers()

        self.ecfp()

    def set_initial_identifiers(self):
        for atom in self.structure.graph:
            if atom.type != 'H' and atom.type != '*':
                self.identifiers[atom] = {}
                self.seen_atoms[atom] = {}
                self.seen_atoms[atom][0] = {atom}

                daylight_properties = Daylight(atom, self.structure)
                initial_identifier = hash_32_bit_integer(daylight_properties.daylight)

                self.identifiers[atom][0] = initial_identifier
                self.fingerprint.add(initial_identifier)
                bonds = set(atom.get_non_hydrogen_bonds())
                self.bonds[initial_identifier] = bonds

                feature = tuple(sorted(list(bonds) + [atom], key=lambda x: (x.nr, x.type)))
                self.features[feature] = (initial_identifier, 0, atom)

    def ecfp(self):

        for i in range(self.iterations):
            print("Iteration: ", i + 1)
            pprint(self.features)
            new_features = []
            for atom in self.identifiers:
                identifier = self.identifiers[atom][i]
                array = [i + 1, identifier]

                array_to_add = []

                neighbouring_bonds = set()

                self.seen_atoms[atom][i + 1] = copy.deepcopy(self.seen_atoms[atom][i])

                for neighbour in atom.get_non_hydrogen_neighbours():
                    bond = self.structure.bond_lookup[atom][neighbour]
                    bond_order = BOND_PROPERTIES.bond_type_to_order[bond.type]
                    neighbour_identifier = self.identifiers[neighbour][i]
                    neighbouring_bonds = neighbouring_bonds.union(self.bonds[neighbour_identifier])
                    array_to_add.append((bond_order, neighbour_identifier, neighbour))

                    for seen_atom in self.seen_atoms[neighbour][i]:
                        self.seen_atoms[atom][i + 1].add(seen_atom)

                array_to_add.sort(key=lambda x: (x[0], x[1]))

                attachment_order = []

                for bond_order, atom_id, neighbour in array_to_add:
                    array.append(bond_order)
                    array.append(atom_id)
                    attachment_order.append(neighbour)

                if atom.chiral:
                    chirality = find_chirality_from_nonh(atom.neighbours, attachment_order, atom.chiral)
                    if chirality == 'clockwise':
                        array.append(1)
                    else:
                        array.append(0)

                new_identifier = hash_32_bit_integer(array)

                self.identifiers[atom][i + 1] = new_identifier
                bonds_core_previous = self.bonds[identifier]
                bonds_attachment = atom.get_non_hydrogen_bonds()

                bond_set = bonds_core_previous.union(bonds_attachment)
                bond_set = bond_set.union(neighbouring_bonds)
                self.bonds[new_identifier] = bond_set

                feature = tuple(sorted(list(bond_set) + list(self.seen_atoms[atom][i + 1]), key=lambda x: (x.nr, x.type)))

                if feature not in self.features:
                    new_features.append((feature, new_identifier, atom))

            new_features.sort(key=lambda x: tuple([y.nr for y in x[0]] + [x[1]]))
            previous_feature = None
            previous_atom = None

            for new_feature, identifier, atom in new_features:
                if new_feature == previous_feature:
                    print('discarded: ', atom, 'in favour of', previous_atom)
                    print(new_feature, identifier)

                    continue
                else:
                    self.features[new_feature] = (identifier, i + 1, atom)
                    self.fingerprint.add(identifier)
                    previous_feature = new_feature
                    previous_atom = atom

        pprint(self.identifiers)
