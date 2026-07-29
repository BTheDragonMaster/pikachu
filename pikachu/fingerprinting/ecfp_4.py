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
        self.disambiguated_bond_stereo = {}

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

                for bond in bonds:
                    if bond.chiral:
                        self.disambiguated_bond_stereo[bond] = False

    def ecfp(self):
        # Nr of iterations determines the radius of the fingerprinting
        for i in range(self.iterations):
            new_features = []

            # Precompute, once per iteration, which stereo bonds become
            # resolved at this iteration. This must happen before the atom
            # loop below: both atoms on a bond need to fold in the same bit
            # during this same iteration, so we can't flip
            # disambiguated_bond_stereo until after both have been processed.

            bond_stereo_bits = {}
            for bond, already_disambiguated in self.disambiguated_bond_stereo.items():
                if not already_disambiguated and self.stereo_bond_resolved(bond, i):
                    bond_stereo_bits[bond] = self.get_stereo_bit(bond, i)

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

                # Fold in cis/trans bits for any stereo double bonds attached to this
                # atom that became resolved this iteration. Sorted by bond.nr so the
                # order bits are appended in is deterministic (an atom could in
                # principle sit on more than one stereo double bond).
                stereo_bond = next(
                    (bond for bond in atom.get_non_hydrogen_bonds() if bond.chiral),
                    None
                )
                if stereo_bond is not None and stereo_bond in bond_stereo_bits:
                    array.append(bond_stereo_bits[stereo_bond])

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

            for bond in bond_stereo_bits:
                self.disambiguated_bond_stereo[bond] = True

            for new_feature, identifier, atom in new_features:
                # TODO: Make a better feature representation, perhaps as SMILES string
                self.fingerprint.add(identifier)
                self.hash_to_feature[identifier] = new_feature

    @staticmethod
    def get_stereo_substituents(atom, partner):
        """
        Return the neighbours of 'atom' on a stereo double bond, excluding the
        double bond partner itself. Normally one or two such substituents exist
        (a missing second slot means an implicit H).
        """
        return [n for n in atom.neighbours if n != partner]

    def stereo_bond_resolved(self, bond, i):
        """
        Check whether both sides of a stereo double bond currently have
        distinguishable substituents at iteration i, i.e. whether we can
        unambiguously identify which substituent is which when consulting
        bond.chiral_dict.
        """
        atom_1 = bond.atom_1
        atom_2 = bond.atom_2

        for atom, partner in ((atom_1, atom_2), (atom_2, atom_1)):
            # Get the neighbours for a stereo bond
            substituents = self.get_stereo_substituents(atom, partner)

            if len(substituents) < 2:
                # Only one real substituent -> nothing to disambiguate
                continue

            sub_identifiers = []
            for substituent in substituents:
                if substituent.type == 'H':
                    sub_identifiers.append('dummy')
                else:
                    sub_identifiers.append(self.identifiers[substituent][i])

            # Can't tell the two substituents on this side apart yet
            if len(sub_identifiers) != len(set(sub_identifiers)):
                return False

        return True

    def get_canonical_stereo_substituent(self, atom, partner, i):
        """
        Return the substituent on this side of a stereo double bond that
        should be used to query bond.chiral_dict, chosen by a rule that is
        invariant across structures (unlike atom.nr, which just reflects
        SMILES parsing order).

        Rule: prefer identifier-based comparison, since self.identifiers[atom][i]
        is a structural hash - the same local chemical environment gets the
        same identifier no matter what molecule it's embedded in or how it
        was numbered. Hydrogen substituents need no comparison: if only one
        heavy substituent exists on this side, it's the only candidate and is
        trivially canonical (the "other slot" being an implicit/explicit H is
        not a real ambiguity to resolve).
        """
        substituents = self.get_stereo_substituents(atom, partner)

        if not substituents:
            raise ValueError("Stereo bonds must have at least one neighbour")

        keyed_substituents = []
        for substituent in substituents:
            if substituent.type == 'H':
                # Sentinel key: sorts after every real identifier, so H is
                # only ever picked as canonical when it's the sole substituent
                # on this side (the tuple's first element, 1, always loses to
                # any real identifier's first element, 0)
                key = (1, 0)
            else:
                key = (0, self.identifiers[substituent][i])
            keyed_substituents.append((key, substituent))

        # Pick the substituent with the lowest key. Real identifiers always
        # beat the H sentinel; between two real identifiers, the structurally
        # "smaller" one wins - a choice that's reproducible across structures
        # since self.identifiers[atom][i] depends only on local environment
        return min(keyed_substituents, key=lambda pair: pair[0])[1]

    def get_stereo_bit(self, bond, i):
        """
        Given a bond with defined cis/trans stereochemistry that is resolved
        at iteration i, return 0 or 1 encoding that relationship.
        """
        substituent_1 = self.get_canonical_stereo_substituent(bond.atom_1, bond.atom_2, i)
        substituent_2 = self.get_canonical_stereo_substituent(bond.atom_2, bond.atom_1, i)

        if substituent_1 is None or substituent_2 is None:
            return None

        relationship = bond.chiral_dict[substituent_1][substituent_2]

        if relationship == 'cis':
            return 1
        elif relationship == 'trans':
            return 0
        else:
            raise ValueError(f"Expected 'cis' or 'trans' in chiral_dict, got {relationship}")


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




