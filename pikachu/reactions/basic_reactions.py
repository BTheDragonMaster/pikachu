from pikachu.general import read_smiles, draw_structure
from pikachu.reactions.functional_groups import find_bonds, BondDefiner, combine_structures
from pikachu.chem.atom_properties import ATOM_PROPERTIES


def hydrolysis(structure, bond):

    # Create the water molecule required for hydrolysis
    water = read_smiles("O")
    water_bond = water.bonds[0]

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
        
    water.break_bond(water_bond)
    h_and_oh_1, h_and_oh_2 = water.split_disconnected_structures()
    
    h = h_and_oh_1
    oh = h_and_oh_2
    
    if len(h_and_oh_1.graph) == 2:
        oh = h_and_oh_1
        h = h_and_oh_2

    # Begin hydrolysis by breaking the desired bond
        
    structure.break_bond(bond)
    product_1, product_2 = structure.split_disconnected_structures()

    # Determine which structure will bond with OH and which with H
    
    oh_acceptor_structure = product_1
    h_acceptor_structure = product_2
    
    if h_acceptor in product_1.graph:
        h_acceptor_structure = product_1
        oh_acceptor_structure = product_2

    # Put structures into a single structure object. This also updates the atom numbers of the original structures

    h_structure = combine_structures([h, h_acceptor_structure])
    oh_structure = combine_structures([oh, oh_acceptor_structure])

    # Obtain the reacting oxygen, hydrogen, OH acceptor and H acceptor atoms from the combined structure

    oxygen = None
    for atom in oh.graph:
        if atom.type == 'O':
            oxygen = oh_structure.get_atom(atom)
            break

    hydrogen = h_structure.get_atom(list(h.graph.keys())[0])

    h_acceptor = h_structure.get_atom(h_acceptor)
    oh_acceptor = oh_structure.get_atom(oh_acceptor)

    # Create the required bonds

    h_structure.make_bond(hydrogen, h_acceptor, h_structure.find_next_bond_nr())
    oh_structure.make_bond(oxygen, oh_acceptor, oh_structure.find_next_bond_nr())

    # Return the hydrolysed product

    return [h_structure, oh_structure]


if __name__ == "__main__":
    smiles = "NCC(=O)NCC(=O)O"
    structure = read_smiles(smiles)
    peptide_bond = BondDefiner("peptide bond", "C(=O)N", 0, 2)
    peptide_bonds = find_bonds(peptide_bond, structure)
    structures = hydrolysis(structure, peptide_bonds[0])

    for s in structures:
        draw_structure(s)
