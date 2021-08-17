from pikachu import *

def combine_structure(structure_1, structure_2):

    for atom in structure_2.graph:
        #print(atom.nr)
        atom.nr += len(structure_1.graph)
        #print(atom.nr)

    return structure_1, structure_2



if __name__ == "__main__":
    string1 = "C1CCN[C@@H](C1)O"
    string2 = "C1CCN[C@@H](C1)O"
    smiles_1 = Smiles(string1)
    #print(smiles_1)
    smiles_2 = Smiles(string2)
    structure_1 = smiles_1.smiles_to_structure()
    #print(structure_1)
    structure_2 = smiles_2.smiles_to_structure()
    print(structure_2.graph)
    #print(structure_1)
    combine_structure(structure_1, structure_2)
    print(structure_2.graph)