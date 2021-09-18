from smiles import Smiles
from errors import SmilesError
from graph_to_smiles import GraphToSmiles
from drawing import Drawer


def read_smiles(smiles_string):
    if smiles_string:
        try:
            smiles = Smiles(smiles_string)
            structure = smiles.smiles_to_structure()
            return structure
        except SmilesError as e:
            print(f'Error parsing "{smiles_string}": {e.message}')
            return


def to_smiles(structure, kekule=False):
    if kekule:
        structure = structure.kekulise()

    return GraphToSmiles(structure).smiles


def draw_structure(structure):
    structure = structure.kekulise()
    Drawer(structure)


def draw_smiles(smiles):
    structure = read_smiles(smiles)
    structure = structure.kekulise()
    Drawer(structure)