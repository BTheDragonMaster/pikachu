import time
import os
from sys import argv
from random import shuffle, seed

from pikachu.general import read_smiles, position_smiles

from rdkit.Chem import MolFromSmiles

import timeout_decorator


def substructure_matching_speed_pikachu(smiles, subsmiles):
    start_time = time.time()
    print("Starting substructure search..")

    for s in smiles:
        structure = read_smiles(s)
        for sub in subsmiles:
            substructure = read_smiles(sub)
            matches = structure.find_substructures(substructure)
            has_substructure.append(len(matches))

    time_1 = time.time()
    print(f'Time spent by PIKAChU finding {len(subsmiles)} substructures in {len(smiles)} structures: {time_1 - start_time}')

    return has_substructure


def drawing_speed_pikachu(smiles):

    start_time = time.time()
    for i, s in enumerate(smiles):
        position_smiles(s)
    time_1 = time.time()
    print(f'Time spent by PIKAChU drawing {len(smiles)} SMILES: {time_1 - start_time}')


def reading_speed_pikachu(smiles):
    start_time = time.time()
    for s in smiles:
        x = read_smiles(s)
    time_1 = time.time()
    print(f'Time spent by PIKAChU reading {len(smiles)} SMILES: {time_1 - start_time}')


def substructure_matching_speed_rdkit(smiles, subsmiles):
    start_time = time.time()
    has_substructure = []
    for s in smiles:
        structure = MolFromSmiles(s)
        for sub in subsmiles:
            substructure = MolFromSmiles(sub)
            matches = structure.GetSubstructMatches(substructure, useChirality=True)
            has_substructure.append(len(matches))

    time_1 = time.time()
    print(f'Time spent by RDKit finding {len(subsmiles)} substructures in {len(smiles)} structures: {time_1 - start_time}')

    return has_substructure


def drawing_speed_rdkit(smiles):
    start_time = time.time()
    for s in smiles:
        x = MolFromSmiles(s)

    time_1 = time.time()
    print(f'Time spent by RDKit drawing {len(smiles)} SMILES: {time_1 - start_time}')


def reading_speed_rdkit(smiles):
    start_time = time.time()
    for s in smiles:
        x = MolFromSmiles(s)

    time_1 = time.time()
    print(f'Time spent by RDKit reading {len(smiles)} SMILES: {time_1 - start_time}')


def compare_substructure_matching_outcomes(pikachu_list, rdkit_list):
    correct = 0
    incorrect = 0
    mistake_indices = []

    for i, pikachu_entry in enumerate(pikachu_list):
        rdkit_entry = rdkit_list[i]
        if pikachu_entry == rdkit_entry:
            correct += 1
        else:
            incorrect += 1
            mistake_indices.append(i)

    return correct, incorrect, mistake_indices


def read_smiles_file(smiles_file):
    smiles_strings = []
    with open(smiles_file, 'r') as s_file:
        for line in s_file:
            smiles = line.strip()
            if smiles:
                smiles_strings.append(smiles)
                
    seed(11)
    shuffle(smiles_strings)
    return smiles_strings


if __name__ == "__main__":
    smiles_file = argv[1]
    supersmiles_file = argv[2]
    subsmiles_file = argv[3]

    smiles_strings = read_smiles_file(smiles_file)
    supersmiles_strings = read_smiles_file(supersmiles_file)
    subsmiles_strings = read_smiles_file(subsmiles_file)

    drawing_speed_pikachu(smiles_strings)
    reading_speed_pikachu(smiles_strings)

    drawing_speed_rdkit(smiles_strings)
    reading_speed_rdkit(smiles_strings)

    pikachu_list = substructure_matching_speed_pikachu(supersmiles_strings, subsmiles_strings)
    rdkit_list = substructure_matching_speed_pikachu(supersmiles_strings, subsmiles_strings)

    correct, incorrect, mistake_indices = compare_substructure_matching_outcomes()

    


