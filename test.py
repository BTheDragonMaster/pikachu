from __future__ import absolute_import

import pikachu.pikachu as pikachu
from pikachu.fingerprinting.ecfp_4 import ECFP

smiles = r'CCCCCCCCCC(=O)N[C@@H](CC1=CNC2=CC=CC=C21)C(=O)N[C@@H](CC(=O)N)C(=O)N[C@@H](CC(=O)O)C(=O)N[C@H]3[C@H](OC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@H](NC(=O)CNC(=O)[C@@H](NC(=O)[C@H](NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)CNC3=O)CCCN)CC(=O)O)C)CC(=O)O)CO)[C@H](C)CC(=O)O)CC(=O)C4=CC=CC=C4N)C'
smiles = r'CCCC(=O)N'
smiles_2 = r'C([C@@H](C(=O)O)N)O'
structure = pikachu.read_smiles(smiles)

#structure_2 = pikachu.read_smiles(smiles_2)
# print(pikachu.to_smiles(structure))

# pikachu.draw_smiles(smiles)
# pikachu.draw_structure(structure)

ecfp = ECFP(structure)
print(ecfp.fingerprint)