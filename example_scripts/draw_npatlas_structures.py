#!/usr/bin/env python

import os
from sys import argv
from pikachu.general import read_smiles, svg_from_smiles
import timeout_decorator


def parse_npatlas_smiles(npatlas_file):
    smiles = []
    with open(npatlas_file, 'r') as npatlas:
        for line in npatlas:
            line = line.strip()
            if line:
                smiles.append(line)

    return smiles

@timeout_decorator.timeout(20)
def draw(smiles, drawing_file):
    svg_from_smiles(smiles, drawing_file)


def draw_npatlas(npatlas_file, drawing_dir, failed_smiles_dir):
    smiles_strings = parse_npatlas_smiles(npatlas_file)
    if not os.path.exists(drawing_dir):
        os.mkdir(drawing_dir)
    if not os.path.exists(failed_smiles_dir):
        os.mkdir(failed_smiles_dir)
    failed_smiles_file = os.path.join(failed_smiles_dir, 'failed_smiles.txt')
    failed_drawing_file = os.path.join(failed_smiles_dir, 'failed_drawings.txt')
    failed_drawings = open(failed_drawing_file, 'w')
    with open(failed_smiles_file, 'w') as failed_smiles:
        for i, smiles in enumerate(smiles_strings):
            try:
                structure = read_smiles(smiles)
                if not structure:
                    print(smiles)
                    failed_smiles.write(f'{smiles}\tSmilesError\n')
                else:
                    drawing_file = os.path.join(drawing_dir, f'{i}.svg')
                    try:
                        draw(smiles, drawing_file)
                    except Exception as e:
                        print("Drawing failure:", smiles)
                        failed_drawings.write(f'{smiles}\t{e}\n')

            except Exception as e:
                print(smiles)
                failed_smiles.write(f'{smiles}\t{e}\n')

            print(f"Handled smiles number {i}: {smiles}.")

    failed_drawings.close()


if __name__ == "__main__":
    draw_npatlas(argv[1], argv[2], argv[3])
