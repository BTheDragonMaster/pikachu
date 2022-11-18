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
        self.disambiguated_chiral = {}
        self.seen_atoms = {}
        self.features = {}
        self.hash_to_feature = {}
        self.set_initial_identifiers()

        self.ecfp()

    def set_initial_identifiers(self):
        for atom in self.structure.graph:
            if atom.type != 'H' and atom.type != '*':
                if atom.chiral:
                    self.disambiguated_chiral[atom] = False

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
                self.hash_to_feature[initial_identifier] = feature

    def ecfp(self):

        for i in range(self.iterations):
            new_features = []
            for atom in self.identifiers:

                identifier = self.identifiers[atom][i]
                array = [i + 1, identifier]

                array_to_add = []

                neighbouring_bonds = []

                seen_atoms = list(self.seen_atoms[atom][i])

                for neighbour in atom.get_non_hydrogen_neighbours():
                    bond = self.structure.bond_lookup[atom][neighbour]
                    bond_order = BOND_PROPERTIES.bond_type_to_order[bond.type]
                    neighbour_identifier = self.identifiers[neighbour][i]
                    for neighbouring_bond in self.bonds[neighbour_identifier]:
                        neighbouring_bonds.append(neighbouring_bond)
                    array_to_add.append((bond_order, neighbour_identifier, neighbour))

                    for seen_atom in self.seen_atoms[neighbour][i]:
                        seen_atoms.append(seen_atom)

                self.seen_atoms[atom][i + 1] = set(seen_atoms)

                neighbouring_bonds = set(neighbouring_bonds)

                array_to_add.sort(key=lambda x: (x[0], x[1]))

                attachment_order = []

                for bond_order, atom_id, neighbour in array_to_add:
                    array.append(bond_order)
                    array.append(atom_id)
                    attachment_order.append(neighbour)

                if atom.chiral and not self.disambiguated_chiral[atom]:
                    neighbour_identifiers = []
                    neighbour_identifiers_sorted = []
                    for neighbour in attachment_order:
                        neighbour_identifier = self.identifiers[neighbour][i]
                        neighbour_identifiers_sorted.append(neighbour_identifier)

                    if len(neighbour_identifiers_sorted) == len(set(neighbour_identifiers_sorted)):

                        for neighbour in atom.neighbours:
                            if neighbour.type == 'H':
                                neighbour_identifier = 'dummy'
                            else:
                                neighbour_identifier = self.identifiers[neighbour][i]
                            neighbour_identifiers.append(neighbour_identifier)

                        chirality = find_chirality_from_nonh(neighbour_identifiers, neighbour_identifiers_sorted,
                                                             atom.chiral)
                        if chirality == 'clockwise':
                            array.append(1)
                        else:
                            array.append(0)

                        self.disambiguated_chiral[atom] = True

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
            #new_features.sort()
            previous_feature = None
            previous_atom = None

            for new_feature, identifier, atom in new_features:
                # TODO: Make a better feature representation, perhaps as SMILES string
                self.fingerprint.add(identifier)
                if new_feature == previous_feature:
                    continue
                else:
                    self.features[new_feature] = (identifier, i + 1, atom)
                    self.hash_to_feature[identifier] = new_feature

                    previous_feature = new_feature
                    previous_atom = atom


def build_ecfp_bitvector(structures, depth=2, bits=1024):
    fingerprints = []
    identifier_to_feature = {}

    for structure in structures:
        ecfp = ECFP(structure, iterations=depth)
        fingerprints.append(ecfp.fingerprint)
        identifier_to_feature.update(ecfp.hash_to_feature)

    substructure_to_count = {}

    for fingerprint in fingerprints:
        for identifier in fingerprint:
            if identifier not in substructure_to_count:
                substructure_to_count[identifier] = 0
            substructure_to_count[identifier] += 1

    substructures = sorted(list(substructure_to_count.items()), key=lambda x: x[1], reverse=True)
    bitvector_substructures = [x[0] for x in substructures[:bits]]
    bitvector_mapping = {}
    for substructure in bitvector_substructures:
        bitvector_mapping[substructure] = identifier_to_feature[substructure]

    return bitvector_substructures, bitvector_mapping




