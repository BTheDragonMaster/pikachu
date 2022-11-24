from sys import argv

from pikachu.general import read_smiles, draw_smiles


def draw_nitrogen_smiles(smiles_file):

    with open(smiles_file, "r") as smi:
        for line in smi:
            smiles = line.strip()
            draw_smiles(smiles)


if __name__ == "__main__":

    draw_nitrogen_smiles(argv[1])
