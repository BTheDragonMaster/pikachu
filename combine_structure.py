from pikachu import *
from groupdefiner import *
from reactions import *
from graph_to_smiles import *
from drawing import *

def correct_structures_index(structures):

    for i in range(len(structures) - 1):
        structure = structures[i]

        upvalue = structure.find_highest_atom_nr()
        upvalue_bond = structure.find_highest_bond_nr()

        structure_to_mend = structures[i + 1]

        for atom in structure_to_mend.graph:
                atom.nr += upvalue + 1
        new_graph = {}
        new_bonds = {}

        for atom, atoms in structure_to_mend.graph.items():

            new_graph[atom] = atoms

        for bond_nr, bond in structure_to_mend.bonds.items():
            bond.nr += upvalue_bond + 1
            new_bonds[bond.nr] = bond

        structure_to_mend.graph = new_graph
        structure_to_mend.bonds = new_bonds
        structure_to_mend.make_bond_lookup()
        structure_to_mend.set_atom_neighbours()

    return structures


def combine_graphs(structures):
    new_graph = {}
    new_bonds = {}
    for group in structures:
        new_graph.update(group.graph)
        new_bonds.update(group.bonds)

    structure = Structure(new_graph, new_bonds)
    structure.make_bond_lookup()
    structure.set_atom_neighbours()
    return structure


if __name__ == "__main__":
    string1 = "C(CCN)C[C@@H](C(=O)O)N"
    string2 = "C1CCN[C@@H](C1)O"
    string3 = "C1CCN[C@@H](C1)O"

    smiles_1 = Smiles(string1)
    smiles_2 = Smiles(string2)
    smiles_3 = Smiles(string3)

    structure_1 = smiles_1.smiles_to_structure()
    structure_2 = smiles_2.smiles_to_structure()
    structure_3 = smiles_3.smiles_to_structure()

    structures = (structure_1, structure_2, structure_3)
    structures = correct_structures_index(structures)
    for structure in structures:
        print(structure.graph)
        print(structure.bonds)

    GraphToSmiles(structures[1])


