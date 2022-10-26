#!/usr/bin/env python

import time
import os

import matplotlib.pyplot as plt

from pikachu.smiles.smiles import Smiles
from pikachu.errors import StructureError, ColourError
from pikachu.smiles.graph_to_smiles import GraphToSmiles
from pikachu.drawing.drawing import Drawer, Options, draw_multiple
from pikachu.drawing.colours import *
from pikachu.chem.molfile.write_molfile import MolFileWriter
from pikachu.chem.structure import Structure


def smiles_from_file(smiles_file, read_all=False):
    if not read_all:
        with open(smiles_file, 'r') as smiles:
            smiles_string = smiles.readline().strip()

        return smiles_string
    else:
        smiles_strings = []
        with open(smiles_file, 'r') as smiles:
            for line in smiles:
                smiles_string = line.strip()
                smiles_strings.append(smiles_string)

        return smiles_strings
    

def read_smiles(smiles_string: str) -> Structure:
    """
    Return structure object from SMILES string

    Input:
    smiles_string: str, SMILES string

    Output:
    Structure object if correct SMILES string was parsed, None otherwise
    """

    if not isinstance(smiles_string, str):
        raise TypeError(f"Expected str, got {type(smiles_string)}")
    elif not smiles_string:
        raise ValueError(f"Can't pass empty string as SMILES.")
    else:
        smiles = Smiles(smiles_string)
        structure = smiles.smiles_to_structure()
        if not structure:
            raise ValueError(f"Could not produce structure for SMILES: {smiles_string}.")
        return structure


def structure_to_smiles(structure, kekule=False):
    """
    Return SMILES string from structure object

    Input:

    structure: Structure object
    kekule: bool, return kekulised SMILES string if True, unkekulised SMILES string if False

    Output:
    str, SMILES string

    """
    if kekule:
        structure = structure.kekulise()

    return GraphToSmiles(structure).smiles


def draw_structure(structure, options=None):
    """
    Display structure from structure object

    Input:
    structure: Structure object

    """

    if not options:
        options = Options()
    drawer = Drawer(structure, options=options)

    drawer.show_molecule()


def position_smiles(smiles):
    """
    Return structure with stored atom coordinates

    Input
    ----------
    smiles: str, SMILES string

    Output
    ----------
    structure: Structure object


    """
    structure = read_smiles(smiles)
    if '.' in smiles:
        drawer = draw_multiple(structure, coords_only=True)
    else:
        drawer = Drawer(structure, coords_only=True)

    return drawer


def draw_smiles(smiles, options=None):
    """
    Display structure from SMILES string

    Input:
    smiles: str, SMILES string

    """

    if not options:
        options = Options()

    structure = read_smiles(smiles)
    if '.' in smiles:
        drawer = draw_multiple(structure, options=options)

    else:

        drawer = Drawer(structure, options=options)
    
    drawer.show_molecule()


def smiles_to_molfile(smiles, molfile, options=None):
    if not options:
        options = Options()

    structure = read_smiles(smiles)
    if '.' in smiles:
        MolFileWriter(structure, molfile, drawing_options=options, multiple=True).write_mol_file()
    else:
        MolFileWriter(structure, molfile, drawing_options=options).write_mol_file()


def svg_from_smiles_timed(smiles, svg_out):
    start_time = time.time()
    print("Start")
    time_1 = time.time()
    print(time_1 - start_time)
    structure = read_smiles(smiles)
    print("reading smiles")
    time_2 = time.time()
    print(time_2 - time_1)
    structure = structure.kekulise()
    print("Kekulising")
    time_3 = time.time()
    print(time_3 - time_2)
    drawer = Drawer(structure)
    print("Drawing")
    time_4 = time.time()
    print(time_4 - time_3)
    drawer.save_svg(svg_out)
    print("Saving")
    time_5 = time.time()
    print(time_5 - time_4)


def svg_from_structure(structure, svg_out, options=None):
    """
    Save structure drawing of Structure object to .svg

    Input:
    structure: Structure object
    svg_out: str, output file name, should end in .svg

    """
    if not options:
        options = Options()

    drawer = Drawer(structure, options=options)
    drawer.save_svg(svg_out)


def svg_string_from_structure(structure, options=None):
    if not options:
        options = Options()

    drawer = Drawer(structure, options=options)
    svg_string = drawer.save_svg_string()

    return svg_string


def png_from_structure(structure, png_out, options=None):
    """
    Save structure drawing of Structure object to .png

    Input:
    structure: Structure object
    png_out: str, output file name, should end in .png

    """
    if not options:
        options = Options()

    drawer = Drawer(structure, options=options)
    drawer.save_png(png_out)


def svg_from_smiles(smiles, svg_out, options=None):
    """
    Save structure drawing of SMILES string to .svg

    Input:
    smiles: str, SMILES string
    svg_out: str, output file name, should end in .svg

    """
    structure = read_smiles(smiles)
    if not options:
        options = Options()

    if '.' in smiles:
        drawer = draw_multiple(structure, options=options)
    else:
        drawer = Drawer(structure, options=options, coords_only=True)
    drawer.write_svg(svg_out)


def png_from_smiles(smiles, png_out, options=None):
    """
    Save structure drawing of SMILES string to .png

    Input:
    smiles: str, SMILES string
    png_out: str, output file name, should end in .png

    """
    structure = read_smiles(smiles)
    if not options:
        options = Options()
    drawer = Drawer(structure, options=options)
    drawer.save_png(png_out)


def highlight_substructure(substructure_smiles, parent_smiles, search_mode='all',
                           colour=None,
                           check_chiral_centres=True,
                           check_bond_chirality=True,
                           visualisation='show',
                           out_file=None):
    """
    Find occurrences of (a) substructure(s) in a parent structure and highlight it in a drawing

    Input:
    substructure_smiles: str, SMILES string of substructure, OR list of str, with each str a SMILES string
    parent_smiles: str, SMILES string of superstructure
    search_mode: str, 'single', 'multiple' or 'all. If single, highlight only the first detected instance of a
        substructure. If 'all', highlight all instances of a substructure. If 'multiple', highlight all instances of
        all substructures, assigning one colour per substructure.
    colour: str, hex colour code, ie #ffffff, colour in which substructure will be highlighted, OR list of str,
        with each str a colour.
        Default: None (RASPBERRY for single/ all matching, RANDOM_PALETTE_2 for multiple matching
    check_chiral_centres: bool, if True, only matches substructure to superstructure if stereochemistry
        of all stereocentres match; if False, matches substructure to superstructure regardless of
        stereochemistry of stereocentres.
    check_bond_chirality: bool, if True, only matches substructure to superstructure if stereochemistry
        of all stereobonds match; if False, matches substructure to superstructure regardless of
        stereochemistry of stereobonds.
    visualisation: str, 'show', 'png', or 'svg'. If 'png' or 'svg', out_file is required.
    out_file: str, output file of png or svg drawing

    """
    assert search_mode in {'all', 'single', 'multiple'}

    if search_mode == 'all' or search_mode == 'single':
        assert type(substructure_smiles) == str
        if colour:
            assert type(colour) in {str}
        else:
            colour = RASPBERRY
    elif search_mode == 'multiple':
        assert type(substructure_smiles) in {list, tuple, set}
        assert type(colour) in {list, tuple, set}

    if search_mode == 'all':
        highlight_subsmiles_all(substructure_smiles, parent_smiles, colour=colour,
                                check_chiral_centres=check_chiral_centres,
                                check_bond_chirality=check_bond_chirality,
                                visualisation=visualisation,
                                out_file=out_file)
    elif search_mode == 'multiple':
        highlight_subsmiles_multiple(substructure_smiles, parent_smiles, colours=colour,
                                     check_chiral_centres=check_chiral_centres,
                                     check_bond_chirality=check_bond_chirality,
                                     visualisation=visualisation,
                                     out_file=out_file)
    elif search_mode == 'single':
        highlight_subsmiles_single(substructure_smiles, parent_smiles, colour=colour,
                                   check_chiral_centres=check_chiral_centres,
                                   check_bond_chirality=check_bond_chirality,
                                   visualisation=visualisation,
                                   out_file=out_file)


def highlight_subsmiles_single(substructure_smiles, parent_smiles, colour=RASPBERRY,
                               check_chiral_centres=True,
                               check_bond_chirality=True,
                               visualisation='show',
                               out_file=None):
    """
    Draw structure with a single occurrence of substructure_smiles highlighted with colour

    Input:
    substructure_smiles: str, SMILES string of substructure
    parent_smiles: str, SMILES string of superstructure
    colour: str, hex colour code, ie #ffffff, colour in which substructure will be highlighted
        Default: inbuilt colour raspberry
    check_chiral_centres: bool, if True, only matches substructure to superstructure if stereochemistry
        of all stereocentres match; if False, matches substructure to superstructure regardless of
        stereochemistry of stereocentres.
    check_bond_chirality: bool, if True, only matches substructure to superstructure if stereochemistry
        of all stereobonds match; if False, matches substructure to superstructure regardless of
        stereochemistry of stereobonds.
    visualisation: str, 'show', 'png', or 'svg'. If 'png' or 'svg', out_file is required.
    out_file: str, output file of png or svg drawing

    """
    child_structure = read_smiles(substructure_smiles)
    parent_structure = read_smiles(parent_smiles)

    if not colour.startswith('#'):
        colour = get_hex(colour)
        
    parent_structure.colour_substructure_single(child_structure, colour=colour,
                                                check_chiral_centres=check_chiral_centres,
                                                check_bond_chirality=check_bond_chirality)

    drawer = Drawer(parent_structure)
    if visualisation == 'show':
        drawer.show_molecule()
    elif visualisation == 'svg':
        assert out_file
        drawer.save_svg(out_file)
    elif visualisation == 'png':
        assert out_file
        drawer.save_png(out_file)


def highlight_subsmiles_all(substructure_smiles, parent_smiles, colour=RASPBERRY,
                            check_chiral_centres=True,
                            check_bond_chirality=True,
                            visualisation='show',
                            out_file=None):
    """
    Draw structure with all occurrences of substructure_smiles highlighted with colour

    Input:
    substructure_smiles: str, SMILES string of substructure
    parent_smiles: str, SMILES string of superstructure
    colour: str, hex colour code, ie #ffffff, colour in which substructure will be highlighted.
        Default: inbuilt colour raspberry
    check_chiral_centres: bool, if True, only matches substructure to superstructure if stereochemistry
        of all stereocentres match; if False, matches substructure to superstructure regardless of
        stereochemistry of stereocentres.
    check_bond_chirality: bool, if True, only matches substructure to superstructure if stereochemistry
        of all stereobonds match; if False, matches substructure to superstructure regardless of
        stereochemistry of stereobonds.
    visualisation: str, 'show', 'png', or 'svg'. If 'png' or 'svg', out_file is required.
    out_file: str, output file of png or svg drawing

    """
    child_structure = read_smiles(substructure_smiles)
    parent_structure = read_smiles(parent_smiles)

    if not colour.startswith('#'):
        colour = get_hex(colour)

    parent_structure.colour_substructure_all(child_structure, colour=colour,
                                             check_chiral_centres=check_chiral_centres,
                                             check_bond_chirality=check_bond_chirality)

    drawer = Drawer(parent_structure)
    if visualisation == 'show':
        drawer.show_molecule()
    elif visualisation == 'svg':
        assert out_file
        drawer.save_svg(out_file)
    elif visualisation == 'png':
        assert out_file
        drawer.save_png(out_file)


def highlight_subsmiles_multiple(substructure_smiles_list, parent_smiles, colours=None,
                                 check_chiral_centres=True,
                                 check_bond_chirality=True,
                                 visualisation='show',
                                 out_file=None):
    """
    Draw structure with all occurrences of all substructure_smiles highlighted in different colours

    Input:
    substructure_smiles_list: list of str, with each str a SMILES string of substructure. Length must be shorter
        than or equal to the length of colours.
    parent_smiles: str, SMILES string of superstructure
    colours: list of str, with each str a hex colour code, ie #ffffff, colours in which substructures will be
        highlighted in order of occurrence. Length must be longer than or equal to the length of
        substructure_smiles_list
    check_chiral_centres: bool, if True, only matches substructure to superstructure if stereochemistry
        of all stereocentres match; if False, matches substructure to superstructure regardless of
        stereochemistry of stereocentres.
    check_bond_chirality: bool, if True, only matches substructure to superstructure if stereochemistry
        of all stereobonds match; if False, matches substructure to superstructure regardless of
        stereochemistry of stereobonds.
    visualisation: str, 'show', 'png', or 'svg'. If 'png' or 'svg', out_file is required.
    out_file: str, output file of png or svg drawing

    """
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
    elif visualisation == 'png':
        assert out_file
        drawer.save_png(out_file)
