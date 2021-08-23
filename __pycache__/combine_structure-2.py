from pikachu import *
from groupdefiner import *
from reactions import *

def correct_structures_index(structures):
    #print(len(structures))
    #corrected_structures = copy.deepcopy(structures)
    for i in range(len(structures)-1):
        structure = structures[i]
        upvalue = structure.find_highest_atom_nr()
        structure_to_mend = structures[i+1]
        for atom in structure_to_mend.graph:
                #print(atom.nr)
                atom.nr += upvalue+1
                #print(atom.nr)
        new_graph = {}
        for atom, atoms in structure_to_mend.graph.items():
            #print(atom)
            #print(atoms)

            new_graph[atom] = atoms

            structure_to_mend.graph = new_graph
            structure_to_mend.make_bond_lookup()

    return structures

def combine_graphs(structures):
    product = Smiles('C').smiles_to_structure()
    for group in structures:
        for atom, atoms in group.graph.items():
            product.graph[atom] = atoms
            product.make_bond_lookup()
    return product


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

    structures = (structure_1, structure_2)
    print(structures[1].graph)
    structures = correct_structures_index(structures)
    print(structures[1].graph)

