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

    return list(set(bonds))


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

    return list(set(atoms))


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

    struct1, struct2 = structures

    # make sure you dont add atoms with the same atom.nr
    atom_nrs = []
    for atom in struct2.graph:
        atom_nrs.append(atom.nr)

    index = 0

    for atom_1 in struct1.graph:
        while index in atom_nrs:
            index += 1

        atom_1.nr = index

        atom_nrs.append(index)
        index += 1

    struct1.refresh_structure()

    # make sure you add no double bond nrs

    bond_nrs = []

    # inventorise bond numbers of structure 2

    for atom in struct2.graph:
        for bond in atom.bonds:
            if bond.nr not in bond_nrs:
                bond_nrs.append(bond.nr)

    new_bond_dict = {}

    bond_idx = 0

    for bond_nr, bond in struct1.bonds.items():
        while bond_idx in bond_nrs:
            bond_idx += 1

        bond.nr = bond_idx
        bond_nrs.append(bond_idx)
        new_bond_dict[bond_idx] = bond
        bond_idx += 1

    struct1.bonds = new_bond_dict

    annotations_1 = struct1.annotations
    annotations_2 = struct2.annotations

    new_annotations = annotations_1.union(annotations_2)

    new_graph = {}
    new_bonds = {}

    for structure in [struct2, struct1]:
        new_graph.update(structure.graph)
        new_bonds.update(structure.bonds)

    new_structure = Structure(new_graph, new_bonds)
    
    for new_annotation, default in new_annotations:
        new_structure.add_attribute(new_annotation, default)
        
    new_structure.make_bond_lookup()
    new_structure.refresh_structure(find_cycles=True)

    return new_structure
