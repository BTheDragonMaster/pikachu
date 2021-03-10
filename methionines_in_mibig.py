#!/usr/bin/env python
from sys import argv
import pikachu

MET_SMILES = "CSCCC(C=O)N"
MET_STRUCTURE = pikachu.read_smiles(MET_SMILES)

def parse_compound_file(compound_dir):
    with open(compound_dir, 'r') as compound_file:
        bgc_to_compound_to_smiles = {}
        for line in compound_file:
            line = line.strip()
            bgc, compound_name, compound_smiles = line.split('\t')

            if not bgc in bgc_to_compound_to_smiles:
                bgc_to_compound_to_smiles[bgc] = {}

            bgc_to_compound_to_smiles[bgc][compound_name] = compound_smiles


    return bgc_to_compound_to_smiles


def find_methionines(bgc_to_compound_to_smiles):
    for bgc, compound_to_smiles in bgc_to_compound_to_smiles.items():
        for compound, smiles in compound_to_smiles.items():
            print(smiles)
            structure = pikachu.read_smiles(smiles)
            matches = structure.find_substructures(MET_STRUCTURE, check_chiral_centres=False, check_chiral_double_bonds=False)
            print(matches)

if __name__ == "__main__":
    compound_dict = parse_compound_file(argv[1])
    find_methionines(compound_dict)