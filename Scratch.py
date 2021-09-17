from pikachu import *
from drawing import *
import copy



structure = Smiles('CCCC').smiles_to_structure()


for atom in structure.graph:
    if atom.nr == 0:
        the_bonds = copy.copy(atom.bonds)
print(the_bonds)
print(len(the_bonds))

for bond in the_bonds:
    print(bond)
    structure.break_bond(bond)

print(structure.graph) #break_bond only breaks bond between C0-C1 and C0-H5?



