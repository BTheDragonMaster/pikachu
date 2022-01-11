from pikachu.smiles.smiles import Smiles
from pikachu.chem.structure import Structure
from pikachu.chem.bond import Bond


class IndexTracker:

    def __init__(self):
        self.atoms = {}
        self.bonds = {}

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
        atom_1 = None
        atom_2 = None

        if bond_neighbourhood.atom_1.type != 'H':
            atom_1 = match.atoms[bond_neighbourhood.atom_1]
        elif bond_neighbourhood.atom_2.type != 'H':
            atom_2 = match.atoms[bond_neighbourhood.atom_2]
            if atom_2.has_neighbour('H'):
                atom_1 = atom_2.get_neighbour('H')
            else:
                continue

        if bond_neighbourhood.atom_2.type != 'H':
            atom_2 = match.atoms[bond_neighbourhood.atom_2]
        elif bond_neighbourhood.atom_1.type != 'H':
            atom_1 = match.atoms[bond_neighbourhood.atom_1]
            if atom_1.has_neighbour('H'):
                atom_2 = atom_1.get_neighbour('H')
            else:
                continue
        if atom_1 and atom_2:
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
        if atom_neighbourhood.atom_1.type != 'H':
            atom = match.atoms[atom_neighbourhood.atom_1]
            atoms.append(atom)
        else:
            neighbour = atom_neighbourhood.atom_1.neighbours[0]
            if neighbour.type != 'H':
                neighbouring_atom = match.atoms[neighbour]
                atom = neighbouring_atom.get_neighbour('H')
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


def combine_structures(structures):
    """
    Returns a single combined Structure object from multiple input Structures
    Structures: tuple of to be combined Structure objects
    """
    index_tracker = IndexTracker()

    new_graph = {}
    new_bonds = {}
    struct1, struct2 = structures
    bonds1 = struct1.bonds
    bonds2 = struct2.bonds
    new_bonds_1 = {}

    #make sure you dont add atoms with the same atom.nr
    atom_nrs = []
    for atom in struct2.graph:
        atom_nrs.append(atom.nr)
    index = 0
    for atom in struct1.graph:
        while index in atom_nrs:
            if index not in atom_nrs:
                atom.nr = index
                break
            index += 1
        else:
            atom.nr = index
            for other_atom in struct1.graph:
                if atom in other_atom.neighbours:
                    atom.nr = index
            atom_nrs.append(index)
            index += 1

    #refresh bond_lookup and structure
    struct1.make_bond_lookup()
    struct1.refresh_structure()

    # make sure you add no double bond nrs
    start = 0
    bonds = []
    bond_nrs = []
    for atom in struct2.graph:
        for bond in atom.bonds:
            if bond.nr not in bond_nrs:
                bond_nrs.append(bond.nr)
    new_nrs = []

    #Bond objects are identified only by Bond.nr, reset each Bond.nr first to
    #nr that is not used in both to-be combined structures, to make sure
    #no bonds are skipped in the loop
    for atom in struct1.graph:
        for bond in atom.bonds[:]:
            bond.nr = 9999999
    for atom in struct1.graph:
        for bond in atom.bonds[:]:
            if bond not in bonds:
                while start in bond_nrs:
                    if start not in bond_nrs and start not in new_nrs:
                        bond.nr = start
                        bonds.append(bond)
                        bond_nrs.append(start)
                        new_nrs.append(start)
                        break
                    start += 1
                else:
                    if start not in new_nrs:
                        bond.nr = start
                        bonds.append(bond)
                        bond_nrs.append(start)
                        new_nrs.append(start)
                        start += 1

    # Refresh second structure
    struct2.get_connectivities()
    struct2.set_connectivities()
    struct2.set_atom_neighbours()
    struct2.make_bond_lookup()

    #added 29-11-21: bond nrs are not refreshed as keys in structure.bonds dict
    struct1.bonds = {}
    for bond_nr, bond in bonds1.items():
        updated_bond_nr = bond.nr
        struct1.bonds[updated_bond_nr] = bond
    bonds1 = struct1.bonds

    #In order to prevent aliasing issues, make all bonds new Bond instances for
    #the bonds in structure 1:
    for nr, bond in bonds1.items():
        atom1 = bond.atom_1
        atom2 = bond.atom_2
        bond_type = bond.type
        bond_summary = bond.bond_summary
        bond_electrons = bond.electrons
        bond_nr = bond.nr
        new_bond = Bond(atom1, atom2, bond_type, bond_nr)
        new_bond.bond_summary = bond_summary
        new_bonds_1[bond_nr] = new_bond
        new_bond.electrons = bond_electrons

    #Make new structure object of the combined structure
    bonds2.update(new_bonds_1)
    for structure in structures:
        new_graph.update(structure.graph)
        new_bonds.update(structure.bonds)
    new_structure = Structure(new_graph, bonds2)
    new_structure.make_bond_lookup()
    new_structure.set_connectivities()
    new_structure.refresh_structure()
    new_structure.find_cycles()
    for bond_nr, bond in new_structure.bonds.items():
        bond.set_bond_summary()
    return new_structure