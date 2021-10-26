#!/usr/bin/env python

from pikachu.smiles.smiles import Smiles
from pikachu.errors import SmilesError, ColourError
from pikachu.smiles.graph_to_smiles import GraphToSmiles
from pikachu.drawing.drawing import Drawer
from pikachu.drawing.colours import *


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
    drawer = Drawer(structure)
    drawer.show_molecule()

def draw_svg_from_smiles(smiles, svg_out):
    structure = read_smiles(smiles)
    structure = structure.kekulise()
    drawer = Drawer(structure)
    drawer.save_svg(svg_out)

def highlight_subsmiles_single(substructure_smiles, parent_smiles, colour=RASPBERRY,
                               check_chiral_centres=True,
                               check_bond_chirality=True,
                               visualisation='show',
                               out_file=None):
    child_structure = read_smiles(substructure_smiles)
    parent_structure = read_smiles(parent_smiles)
    parent_structure.colour_substructure_single(child_structure, colour=colour,
                                                check_chiral_centres=check_chiral_centres,
                                                check_bond_chirality=check_bond_chirality)

    drawer = Drawer(parent_structure)
    if visualisation == 'show':
        drawer.show_molecule()
    elif visualisation == 'svg':
        assert out_file
        drawer.save_svg(out_file)

def highlight_subsmiles_all(substructure_smiles, parent_smiles, colour=RASPBERRY,
                            check_chiral_centres=True,
                            check_bond_chirality=True,
                            visualisation='show',
                            out_file=None):
    child_structure = read_smiles(substructure_smiles)
    parent_structure = read_smiles(parent_smiles)
    parent_structure.colour_substructure_all(child_structure, colour=colour,
                                             check_chiral_centres=check_chiral_centres,
                                             check_bond_chirality=check_bond_chirality)

    drawer = Drawer(parent_structure)
    if visualisation == 'show':
        drawer.show_molecule()
    elif visualisation == 'svg':
        assert out_file
        drawer.save_svg(out_file)


def highlight_subsmiles_multiple(substructure_smiles_list, parent_smiles, colours=None,
                                 check_chiral_centres=True,
                                 check_bond_chirality=True,
                                 visualisation='show',
                                 out_file=None):
    parent_structure = read_smiles(parent_smiles)

    smiles_nr = len(substructure_smiles_list)

    if not colours:
        colour_list = RANDOM_PALETTE_2[:smiles_nr]

    else:
        colour_list = []
        for colour in colours:
            hex_colour = get_hex(colour)
            colour_list.append(hex_colour)

        colour_list = colour_list[:smiles_nr]

    try:
        assert len(colour_list) == smiles_nr
    except AssertionError:
        raise ColourError('too few colours')

    for i, smiles in enumerate(substructure_smiles_list):
        print("Smiles:", smiles)
        child_structure = read_smiles(smiles)
        colour = colour_list[i]
        parent_structure.colour_substructure_all(child_structure, colour=colour,
                                                 check_chiral_centres=check_chiral_centres,
                                                 check_bond_chirality=check_bond_chirality)

    drawer = Drawer(parent_structure)
    if visualisation == 'show':
        drawer.show_molecule()
    elif visualisation == 'svg':
        assert out_file
        drawer.save_svg(out_file)
