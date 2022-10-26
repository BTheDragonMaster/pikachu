from pikachu.general import read_smiles, draw_structure
from pikachu.reactions.functional_groups import find_bonds, BondDefiner, combine_structures, GroupDefiner, find_atoms
from pikachu.chem.atom_properties import ATOM_PROPERTIES


def internal_condensation(structure, oh_bond, h_bond):
    # Ensure that the defined OH-bond actually has an -OH leaving group
    o_atom = None
    o_neighbour = None

    for atom in oh_bond.neighbours:
        if not o_atom and atom.type == 'O' and atom.has_neighbour('H'):
            o_atom = atom
        else:
            o_neighbour = atom

    if not o_atom:
        raise Exception(f"The selected bond {oh_bond} does not attach to an -OH leaving group.")

    # Ensure that the defined H-bond actually has an -H leaving group

    h_atom = None
    h_neighbour = None

    for atom in h_bond.neighbours:
        if not h_atom and atom.type == 'H':
            h_atom = atom
        else:
            h_neighbour = atom

    if not h_atom:
        raise Exception(f"The selected bond {h_bond} does not attach to an -H leaving group.")

    # Break the bonds between the leaving groups and the rest of the product

    structure.break_bond(oh_bond)
    structure.break_bond(h_bond)

    # Retrieve the atoms that have to form bonds: h_atom with o_atom to form water, and h_neighbour with o_neighbour
    # to form the product

    h_atom = structure.get_atom(h_atom)
    o_atom = structure.get_atom(o_atom)
    h_neighbour = structure.get_atom(h_neighbour)
    o_neighbour = structure.get_atom(o_neighbour)

    # Create the bonds

    structure.make_bond(h_atom, o_atom, structure.find_next_bond_nr())
    structure.make_bond(h_neighbour, o_neighbour, structure.find_next_bond_nr())

    # Put the water and the product into different Structure instances

    structures = structure.split_disconnected_structures()

    water = None
    product = None

    # Find out which of the structures is your product and which is your water

    for structure in structures:
        if h_atom in structure.graph:
            water = structure
        elif h_neighbour in structure.graph:
            structure.refresh_structure(find_cycles=True)
            product = structure

    return [product, water]


def hydrolysis(structure, bond):
    """
    Return hydrolysed product from a structure and the bond to be hydrolysed

    Input
    ----------
    structure: Structure object, structure to be hydrolysed
    bond: Bond object, bond to be hydrolysed

    Returns
    -------
    list of [h_structure, oh_structure], with h_structure the product that gained an 'H' atom from the hydrolysis,
        and oh_structure the product that gained an 'OH' group.

    """

    # Create the water molecule required for hydrolysis
    water = read_smiles("O")
    oxygen = water.atoms[0]
    water_bond = water.bonds[0]
    hydrogen = water_bond.get_connected_atom(oxygen)

    # Determine which atom is more electronegative, which will therefore bond with 'H' rather than 'OH'

    atom_1 = bond.atom_1
    atom_2 = bond.atom_2

    electronegativity_1 = ATOM_PROPERTIES.element_to_electronegativity[atom_1.type]
    electronegativity_2 = ATOM_PROPERTIES.element_to_electronegativity[atom_2.type]

    h_acceptor = atom_1
    oh_acceptor = atom_2

    if electronegativity_2 > electronegativity_1:
        h_acceptor = atom_2
        oh_acceptor = atom_1

    # Break water up into H and OH groups

    hydrolysed_structure = combine_structures([structure, water])

    hydrolysed_structure.break_bond(water_bond)
    hydrolysed_structure.break_bond(bond)

    hydrolysed_structure.make_bond(hydrogen, h_acceptor, hydrolysed_structure.find_next_bond_nr())
    hydrolysed_structure.make_bond(oxygen, oh_acceptor, hydrolysed_structure.find_next_bond_nr())

    structures = hydrolysed_structure.split_disconnected_structures()
    if len(structures) == 1:
        structures[0].refresh_structure(find_cycles=True)

    return structures


def condensation(structure_1, structure_2, oh_bond, h_bond):
    """
    Returns condensed product of structure 1 and structure 2, and a water molecule

    Input
    ----------
    structure_1: Structure object, structure containing -OH leaving group
    structure_2: Structure object, structure containing -H leaving group
    oh_bond: Bond object, bond attached to -OH leaving group
    h_bond: Bond object, bond attached to -H leaving group

    Returns
    -------
    list of [product, water], both Structure objects, with product the condensed product of structure 1 and structure 2
        and water a water molecule.

    """

    # Ensure that the defined OH-bond actually has an -OH leaving group
    o_atom = None
    o_neighbour = None

    for atom in oh_bond.neighbours:
        if not o_atom and atom.type == 'O' and atom.has_neighbour('H'):
            o_atom = atom
        else:
            o_neighbour = atom

    if not o_atom:
        raise Exception(f"The selected bond {oh_bond} does not attach to an -OH leaving group.")

    # Ensure that the defined H-bond actually has an -H leaving group

    h_atom = None
    h_neighbour = None

    for atom in h_bond.neighbours:
        if not h_atom and atom.type == 'H':
            h_atom = atom
        else:
            h_neighbour = atom

    if not h_atom:
        raise Exception(f"The selected bond {h_bond} does not attach to an -H leaving group.")

    # Break the bonds between the leaving groups and the rest of the product

    structure_1.break_bond(oh_bond)
    structure_2.break_bond(h_bond)

    # Put the two structures into the same object

    structure = combine_structures([structure_1, structure_2])

    # Retrieve the atoms that have to form bonds: h_atom with o_atom to form water, and h_neighbour with o_neighbour
    # to form the product

    h_atom = structure.get_atom(h_atom)
    o_atom = structure.get_atom(o_atom)
    h_neighbour = structure.get_atom(h_neighbour)
    o_neighbour = structure.get_atom(o_neighbour)

    # Create the bonds

    structure.make_bond(h_atom, o_atom, structure.find_next_bond_nr())
    structure.make_bond(h_neighbour, o_neighbour, structure.find_next_bond_nr())

    # Put the water and the product into different Structure instances

    structures = structure.split_disconnected_structures()

    water = None
    product = None

    # Find out which of the structures is your product and which is your water

    for structure in structures:
        if h_atom in structure.graph:
            water = structure
        elif h_neighbour in structure.graph:
            structure.refresh_structure(find_cycles=True)
            product = structure

    return [product, water]

#
# if __name__ == "__main__":
#     smiles = "NCC(=O)NCC(=O)O"
#     structure = read_smiles(smiles)
#     peptide_bond = BondDefiner("peptide bond", "C(=O)N", 0, 2)
#     peptide_bonds = find_bonds(peptide_bond, structure)
#     structures = hydrolysis(structure, peptide_bonds[0])
#     glycine = "NCC(=O)O"
#     alanine = "NC(C)C(=O)O"
#
#
#
#     c_oh_bond = BondDefiner("c_oh_bond", "C(=O)O", 0, 2)
#     # n_h_bond = BondDefiner("n_h_bond", "[H]NCC(=O)", 0, 1)
#     n_h_bond = BondDefiner("n_h_bond", "NCC(=O)", 0, 4)
#
#     structure_1 = read_smiles(glycine)
#     structure_2 = read_smiles(alanine)
#
#     c_oh_glycine = find_bonds(c_oh_bond, structure_1)[0]
#     n_h_alanine = find_bonds(n_h_bond, structure_2)[0]
#
#     product, water = condensation(structure_1, structure_2, c_oh_glycine, n_h_alanine)
#
#     draw_structure(product)
