from __future__ import absolute_import
from sys import argv
from pprint import pprint

import pikachu.general as pikachu
from pikachu.fingerprinting.ecfp_4 import ECFP, build_ecfp_bitvector
from pikachu.parsers.np_atlas import smiles_from_npatlas_tabular

def build_npatlas_bitvector(npatlas_file):
    npaid_to_smiles = smiles_from_npatlas_tabular(npatlas_file)
    structures = []
    counter = 0
    for smiles in npaid_to_smiles.values():
        counter += 1

        try:

            structure = pikachu.read_smiles(smiles)
            structures.append(structure)
            print('Done', smiles)
        except Exception as e:
            print(e, smiles)

        if counter > 100:
            break

    structures = [pikachu.read_smiles(r'CCCC(=O)NCCCC(=O)N')] * 500 + [pikachu.read_smiles(r'C(=O)(N)CCC')] * 500

    bitvector, mapping = build_ecfp_bitvector(structures)
    return bitvector, mapping


smiles = r'CCCCCCCCCC(=O)N[C@@H](CC1=CNC2=CC=CC=C21)C(=O)N[C@@H](CC(=O)N)C(=O)N[C@@H](CC(=O)O)C(=O)N[C@H]3[C@H](OC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@H](NC(=O)CNC(=O)[C@@H](NC(=O)[C@H](NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)CNC3=O)CCCN)CC(=O)O)C)CC(=O)O)CO)[C@H](C)CC(=O)O)CC(=O)C4=CC=CC=C4N)C'
#smiles = r'CCCC(=O)N'
#smiles_2 = r'C([C@@H](C(=O)O)N)O'
structure = pikachu.read_smiles(smiles)

#structure_2 = pikachu.read_smiles(smiles_2)
# print(pikachu.to_smiles(structure))

# pikachu.draw_smiles(smiles)
# pikachu.draw_structure(structure)

# ecfp = ECFP(structure)
# print(ecfp.fingerprint)
#
# bitvector, mapping = build_npatlas_bitvector(argv[1])
# print(bitvector)
# pprint(mapping)



if __name__ == "__main__":
    structure = pikachu.read_smiles(r"C1CC(=O)[C@H]2[C@@H]([C@H]1C[C@@H](C(=O)O)N)O2")
    for cycle in structure.cycles.unique_cycles:
        print(cycle)
    pikachu.draw_smiles(r"C1CC(=O)[C@H]2[C@@H]([C@H]1C[C@@H](C(=O)O)N)O2")
    # structure = pikachu.read_smiles(r"C12=NNC=C1C3N(C(C2)C4=C3C=CC=C4)[SH5]")
    # pikachu.draw_smiles(r"C[C@H]1[C@H]([C@@](C[C@@H](O1)O[C@@H]2[C@H]([C@@H]([C@H](O[C@H]2OC3=C4C=C5C=C3OC6=C(C=C(C=C6)[C@H]([C@H](C(=O)N[C@H](C(=O)N[C@H]5C(=O)N[C@@H]7C8=CC(=C(C=C8)O)C9=C(C=C(C=C9O)O)[C@H](NC(=O)[C@H]([C@@H](C1=CC(=C(O4)C=C1)Cl)O)NC7=O)C(=O)O)CC(=O)N)NC(=O)[C@@H](CC(C)C)NC)O)Cl)CO)O)O)(C)N)O")