from pikachu.fingerprinting.daylight import Daylight
from pikachu.fingerprinting.hashing import hash_32_bit_integer
from pikachu.chem.bond_properties import BOND_PROPERTIES
from pikachu.chem.chirality import find_chirality_from_nonh
from pikachu.chem.substructure import Substructure
from pikachu.chem.atom import Atom
from pikachu.chem.bond import Bond
from pikachu.drawing.drawing import Drawer


class Feature:
    def __init__(self, identifier, features):
        self.identifier = identifier
        self.features = features
        self.substructure = self.to_structure()

    def draw(self, svg_out):
        drawer = Drawer(self.substructure, coords_only=True, kekulise=False)
        drawer.write_svg(svg_out)

    def to_structure(self):

        bonds = []
        atoms = []
        bond_nr_to_bond = {}
        atom_nr_to_atom = {}

        for component in self.features:
            if type(component) == Bond:
                bonds.append(component)
            elif type(component) == Atom:
                atoms.append(component)
            else:
                raise ValueError(f"Expected bond or atom type, got {type(component)}")
        structure_graph = {}

        for atom in atoms:
            atom_copy = Atom(atom.type, atom.nr, None, atom.charge, atom.aromatic)
            atom_nr_to_atom[atom_copy.nr] = atom_copy

        for bond in bonds:
            atom_1, atom_2 = None, None
            if bond.atom_1.nr in atom_nr_to_atom:
                atom_1 = atom_nr_to_atom[bond.atom_1.nr]
            if bond.atom_2.nr in atom_nr_to_atom:
                atom_2 = atom_nr_to_atom[bond.atom_2.nr]

            if atom_1 is None:
                atom_1 = Atom('*', bond.atom_1.nr, chiral=None, charge=0, aromatic=bond.atom_1.aromatic)
            if atom_2 is None:
                atom_2 = Atom('*', bond.atom_2.nr, chiral=None, charge=0, aromatic=bond.atom_2.aromatic)

            if atom_1 not in structure_graph:
                structure_graph[atom_1] = []
            if atom_2 not in structure_graph:
                structure_graph[atom_2] = []

            structure_graph[atom_1].append(atom_2)
            structure_graph[atom_2].append(atom_1)

            bond_copy = Bond(atom_1, atom_2, bond.type, bond.nr)
            bond_nr_to_bond[bond.nr] = bond_copy

        substructure = Substructure(structure_graph, bond_nr_to_bond)
        substructure.make_bond_lookup()
            
        return substructure


class ECFP:
    def __init__(self, structure, iterations=2):
        self.structure = structure
        self.iterations = iterations

        self.identifiers = {}
        self.bonds = {}
        self.fingerprint = set()
        self.disambiguated_chiral = {}
        self.seen_atoms = {}
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
                self.bonds[atom] = {}
                self.bonds[atom][0] = bonds

                feature = sorted(list(bonds) + [atom], key=lambda x: (x.nr, x.type))

                self.hash_to_feature[initial_identifier] = Feature(initial_identifier, feature)

    def ecfp(self):
        # Nr of iterations determines the radius of the fingerprinting
        for i in range(self.iterations):
            new_features = []
            for atom in self.identifiers:

                # Start with the hash assigned to the atom in the previous iteration

                identifier = self.identifiers[atom][i]

                # Initialise array that will be hashed with the previous hash, as well as the index of the current hash
                array = [i + 1, identifier]

                array_to_add = []

                neighbouring_bonds = []

                seen_atoms = list(self.seen_atoms[atom][i])

                # Iterate over the atom's direct non-hydrogen neighbours

                for neighbour in atom.get_non_hydrogen_neighbours():

                    # Get the bond between the atom and its neighbour
                    bond = self.structure.bond_lookup[atom][neighbour]

                    # Turn bond type into an integer for hashing
                    bond_order = BOND_PROPERTIES.bond_type_to_order[bond.type]

                    # Obtain the identifier assigned to the atom's neighbour in the previous iteration.
                    neighbour_identifier = self.identifiers[neighbour][i]

                    # Store bonds that are attached to the neighbour.
                    for neighbouring_bond in self.bonds[neighbour][i]:
                        neighbouring_bonds.append(neighbouring_bond)

                    # Store atom and bond information
                    array_to_add.append((bond_order, neighbour_identifier, neighbour))

                    for seen_atom in self.seen_atoms[neighbour][i]:
                        seen_atoms.append(seen_atom)

                self.seen_atoms[atom][i + 1] = set(seen_atoms)

                neighbouring_bonds = set(neighbouring_bonds)

                # Sort the bonds by bond order first, and neighbour identifier next
                array_to_add.sort(key=lambda x: (x[0], x[1]))

                # Store the order in which atoms are added to the hash for chirality disambiguation
                attachment_order = []

                # Extend the initial to-be-hashed array with information on neighbouring atoms and bonds
                for bond_order, atom_id, neighbour in array_to_add:

                    # Store integer representing bond type
                    array.append(bond_order)

                    # Store hashed identifier of neighbouring atom
                    array.append(atom_id)

                    # Store atom for chirality disambiguation
                    attachment_order.append(neighbour)

                # If the atom is chiral and not yet disambiguated, we need to check if the chirality can be
                # disambiguated at this level
                if atom.chiral and not self.disambiguated_chiral[atom]:
                    neighbour_identifiers = []
                    neighbour_identifiers_sorted = []

                    for neighbour in attachment_order:
                        neighbour_identifier = self.identifiers[neighbour][i]
                        neighbour_identifiers_sorted.append(neighbour_identifier)

                    # Chirality is resolved when no duplicate neighbour identifiers exist
                    if len(neighbour_identifiers_sorted) == len(set(neighbour_identifiers_sorted)):

                        # Make a list of neighbour identifiers for chirality determination
                        for neighbour in atom.neighbours:
                            # Store hydrogens as dummies
                            if neighbour.type == 'H':
                                neighbour_identifier = 'dummy'
                            else:
                                neighbour_identifier = self.identifiers[neighbour][i]

                            neighbour_identifiers.append(neighbour_identifier)

                        # Determine chirality
                        chirality = find_chirality_from_nonh(neighbour_identifiers, neighbour_identifiers_sorted,
                                                             atom.chiral)
                        if chirality == 'clockwise':
                            array.append(1)
                        else:
                            array.append(0)

                        # Make sure atom only gets disambiguated once
                        self.disambiguated_chiral[atom] = True

                # New hash is made from all hashes from all previous states
                new_identifier = hash_32_bit_integer(array)

                # Store the new identifier for subsequent rounds
                self.identifiers[atom][i + 1] = new_identifier

                # List of bonds associated with previous atom
                bonds_core_previous = self.bonds[atom][i]
                bonds_attachment = atom.get_non_hydrogen_bonds()

                bond_set = bonds_core_previous.union(bonds_attachment)
                bond_set = bond_set.union(neighbouring_bonds)
                self.bonds[atom][i + 1] = bond_set

                feature = sorted(list(bond_set) + list(self.seen_atoms[atom][i + 1]), key=lambda x: (x.nr, x.type))

                new_features.append((Feature(new_identifier, feature), new_identifier, atom))

            for new_feature, identifier, atom in new_features:
                # TODO: Make a better feature representation, perhaps as SMILES string
                self.fingerprint.add(identifier)
                self.hash_to_feature[identifier] = new_feature


def build_ecfp_bitvector(structures, depth=2, bits=1024):
    fingerprints = []
    identifier_to_feature = {}

    for structure in structures:
        ecfp = ECFP(structure, iterations=depth)
        fingerprints.append(ecfp.fingerprint)
        identifier_to_feature.update(ecfp.hash_to_feature)

    identifier_to_count = {}

    for fingerprint in fingerprints:
        for identifier in fingerprint:
            if identifier not in identifier_to_count:
                identifier_to_count[identifier] = 0
            identifier_to_count[identifier] += 1

    identifiers_and_counts = sorted(list(identifier_to_count.items()), key=lambda x: x[1], reverse=True)
    bitvector_identifiers = [x[0] for x in identifiers_and_counts[:bits]]
    bitvector_mapping = {}
    for identifier in bitvector_identifiers:
        bitvector_mapping[identifier] = identifier_to_feature[identifier]

    return bitvector_identifiers, bitvector_mapping, fingerprints




