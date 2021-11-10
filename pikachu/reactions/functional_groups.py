from pikachu.smiles.smiles import Smiles


def find_bonds(bond_neighbourhood, structure):
    """
    Return list of bond occurrences in structure from a bond neighbourhood

    Input:
    bond_neighbourhood: BondDefiner object

    Output:
    bonds: list of Bond objects

    """

    locations = structure.find_substructures(bond_neighbourhood.structure)
    bonds = []
    for match in locations:
        atom_1 = match.atoms[bond_neighbourhood.atom_1]
        atom_2 = match.atoms[bond_neighbourhood.atom_2]
        bond = structure.bond_lookup[atom_1][atom_2]
        bonds.append(bond)

    return bonds


def find_atoms(atom_neighbourhood, structure):
    """
    Return list of atom occurrences in structure from an atom neighbourhood

    Input:
    atom_neighbourhood: GroupDefiner object

    Output:
    atoms: list of Atom objects

    """
    locations = structure.find_substructures(atom_neighbourhood.structure)
    atoms = []
    for match in locations:
        atom = match.atoms[atom_neighbourhood.atom_1]
        atoms.append(atom)

    return atoms


class BondDefiner:
    """
    Class to store an object that defines a bond type

    Attributes
    ----------
    name: str, name of the bond type
    smiles: Smiles object indicating the structure of the bond
    atom_1: int, index of atom 1 in the bond in the smiles string
    atom_2: int, index of atom 2 in the bond in the smiles string
    """
    def __init__(self, name, smiles, atom_nr_1, atom_nr_2):
        self.name = name
        self.smiles = Smiles(smiles)
        self.structure = self.smiles.smiles_to_structure()
        self.find_atoms(atom_nr_1, atom_nr_2)
        self.bond = self.structure.bond_lookup[self.atom_1][self.atom_2]

    def __repr__(self):
        return self.name

    def find_atoms(self, atom_nr_1, atom_nr_2):
        self.atom_1 = None
        self.atom_2 = None
        for atom in self.structure.graph:
            if atom.nr == atom_nr_1:
                self.atom_1 = atom
            elif atom.nr == atom_nr_2:
                self.atom_2 = atom
            if self.atom_1 and self.atom_2:
                break
        if not self.atom_1 or not self.atom_2:
            raise Exception("Can't find atoms adjacent to bond.")


class GroupDefiner:
    def __init__(self, name, smiles, atom_nr_1):
        self.name = name
        self.smiles = Smiles(smiles)
        self.structure = self.smiles.smiles_to_structure()
        self.find_atom(atom_nr_1)

    def __repr__(self):
        return self.name

    def find_atom(self, atom_nr_1):

        self.atom_1 = None

        for atom in self.structure.graph:
            if atom.nr == atom_nr_1:
                self.atom_1 = atom

            if self.atom_1:
                break

        if not self.atom_1:
            raise Exception("Can't find atoms adjacent to bond.")