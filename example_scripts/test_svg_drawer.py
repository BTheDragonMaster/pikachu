from pikachu.general import svg_from_smiles
from sys import argv


def test_drawer(smiles_string):
    svg_from_smiles(smiles_string, "test.svg")


if __name__ == "__main__":
    test_drawer(argv[1])
