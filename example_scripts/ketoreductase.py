#!/usr/bin/env python

from pikachu.reactions.functional_groups import BondDefiner
from pikachu.general import read_smiles, draw_structure

RECENT_ELONGATION = BondDefiner('recent_elongation', 'O=CCC(=O)S', 0, 1)

def carbonyl_to_hydroxyl(double_bond):
    """Alters the double bond in a carbonyl group to a single bond

    double_bond: Bond object of the bond between the C and O atom of a carbonyl
    group
    """
    #Assert the correct bond is parsed
    assert double_bond.type == 'double'
    atom_types_in_bond = []
    for atom in double_bond.neighbours:
        atom_types_in_bond.append(atom.type)
    assert 'O' in atom_types_in_bond and 'C' in atom_types_in_bond

    # Define the two electrons in the pi bond (inside the db) that needs
    # to be broken
    electrons_in_db = double_bond.electrons
    for electron in electrons_in_db:
        if electron.atom.type == 'O' and electron.orbital_type == 'p':
            electron_1 = electron
        elif electron.atom.type == 'C' and electron.orbital_type == 'p':
            electron_2 = electron

    #Remove the pi electrons from their respective orbital
    orbital_1 = electron_1.orbital
    orbital_2 = electron_2.orbital
    orbital_1.remove_electron(electron_2)
    orbital_2.remove_electron(electron_1)

    #Set bond type to single
    for bond in double_bond.atom_1.bonds:
        if bond == double_bond:
            bond.type = 'single'
    for bond in double_bond.atom_2.bonds:
        if bond == double_bond:
            bond.type = 'single'

    #Remove pi electrons from the Bond between the C and O atom
    for electron in double_bond.electrons[:]:
        if electron == electron_1 or electron == electron_2:
            double_bond.electrons.remove(electron)

    #Change hybridisation of both C and O atoms to sp3
    atom_1, atom_2 = double_bond.neighbours
    atom_1.valence_shell.dehybridise()
    atom_1.valence_shell.hybridise('sp3')
    atom_2.valence_shell.dehybridise()
    atom_2.valence_shell.hybridise('sp3')

    #Change bond_type of Bond instance to single
    double_bond.type = 'single'

    new_single_bond = double_bond
    new_single_bond.set_bond_summary()
    return new_single_bond


def find_betaketon(structure):
    """
    Returns the bond between the carbonyl oxygen and -carbon of the beta ketone
    group

    structure: Structure object
    """
    locations = structure.find_substructures(RECENT_ELONGATION.structure)
    bonds = []
    for match in locations:
        atom_1 = match.atoms[RECENT_ELONGATION.atom_1]
        atom_2 = match.atoms[RECENT_ELONGATION.atom_2]
        bond = structure.bond_lookup[atom_1][atom_2]
        bonds.append(bond)
    return bonds


def ketoreductase(chain_intermediate):
    """
    Performs the ketoreductase reaction on the PKS chain intermediate, returns
    the reaction product as a Structure object

    chain_intermediate: Structure object of a PKS chain intermediate just
    after an elongation step using (methyl)malonyl-CoA
    """
    #Reset all colours to black:
    for atom in chain_intermediate.graph:
        atom.draw.colour = 'black'
    for bond_nr, bond in chain_intermediate.bonds.items():
        bond.set_bond_summary()

    #Identify beta-ketone bond, identify O- and C-atom participating in bond
    beta_ketone_bond = find_betaketon(chain_intermediate)
    for bond in beta_ketone_bond:
        for atom in bond.neighbours:
            if atom.type == 'O':
                carbonyl_oxygen = atom
            elif atom.type == 'C':
                carbonyl_carbon = atom
            else:
                raise Exception('Cannot find atoms in beta ketone bond')

    #Change carbonyl bond to single bond
    for bond in beta_ketone_bond:
        new_single_bond = carbonyl_to_hydroxyl(bond)


    #Add H atom to form hydroxyl group and another H to the C
    chain_intermediate.add_atom('H', [carbonyl_oxygen])
    chain_intermediate.add_atom('H', [carbonyl_carbon])
    carbonyl_carbon.chiral = 'counterclockwise'

    #Set bond summary for newly formed bond (cannot do from struct.bonds?)
    for atom in chain_intermediate.graph:
        if atom == carbonyl_carbon:
            for bond in atom.bonds:
                for neighbour in bond.neighbours:
                    if neighbour.type == 'O':
                        the_bond = bond

    the_bond.set_bond_summary()

    # Refresh structure :)
    chain_intermediate.set_atom_neighbours()
    chain_intermediate.find_cycles()
    for atom in chain_intermediate.graph:
        atom.set_connectivity()
    for bond_nr, bond in chain_intermediate.bonds.items():
        bond.set_bond_summary()

    #Add colouring to the tailored group
    for atom in new_single_bond.neighbours:
        atom.draw.colour = 'red'

    return chain_intermediate

if __name__ == "__main__":
    structure = read_smiles('SC(=O)CC(=O)CCCCC')
    draw_structure(structure)
    reduced_structure = ketoreductase(structure)
    draw_structure(reduced_structure)
