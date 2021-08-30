from pikachu import *
from combine_structure import *
from drawing import *
from graph_to_smiles import *

def remove_PPi(ATP):
    reactive_P = find_group(ATP, AMP_PHOSPHOR)
    #print(reactive_P)
    active_P = reactive_P[0]
    PPi_Oxygen = find_group(ATP, PPi_O)
    #print(PPi_Oxygen)
    for atom in PPi_Oxygen:
        for neighbour in atom.neighbours:
            if neighbour.type == 'C':
                PPi_Oxygen.remove(atom)
                #print(PPi_Oxygen)
    active_Oxygen = PPi_Oxygen[0]

    ATP.break_bond_between_atoms(active_P, active_Oxygen)
    products = ATP.split_disconnected_structures()

    if len(products[0].bonds) > len(products[1].bonds):
        AMP_product = products[0]
        PPi_product = products[1]
        #print('this is AMP')

    elif len(products[1].bonds) > len(products[0].bonds):
        AMP_product = products[1]
        #print('no this is AMP')

    else:
        print('ERROR ERROR BEEP BEEP')

    #print(AMP_product.graph)

    return AMP_product

def prepare_oxygen(AA):
    #print(AA.graph)
    peptide_oxygens = find_group(AA, PEPTIDE_OXYGEN)

    if len(peptide_oxygens) > 1:
        print('MULTIPLE PEPTIDE BONDS AVAILABLE')

    oxygen = peptide_oxygens[0]

    for neighbour in oxygen.neighbours:
        if neighbour.type == 'H':
            #print(neighbour)
            hydrogen = neighbour
            AA.break_bond_between_atoms(oxygen, hydrogen)
            AA.remove_atom(hydrogen)
            break
    #print(oxygen)

    #print(AA.graph)

    return AA, oxygen


def do_adenylation(ATP, AA):
    AMP = remove_PPi(ATP)

    #AMP = activate_AMP(AMP_structure)

    #print(AMP.graph)

    active_AA, active_O = prepare_oxygen(AA)
    #print(active_O)
    structures = [AMP, active_AA]
    correct_structures_index(structures)

    reactive_P = find_group(AMP, AMP_PHOSPHOR)
    active_P = reactive_P[0]

    product = combine_graphs(structures)

    #DEZE DOET HET
    products = product.split_disconnected_structures()

    for mol in products:
        #mol.print_graph()
        #Drawer(mol)
        pass

   # print(product.graph)
   # print(active_O)
   # print(active_P)

    product.refresh_structure()
    next_nr = product.find_next_bond_nr()
    product.make_bond(active_O, active_P, next_nr)
    product.set_atom_neighbours()
    product.make_bond_lookup()
    product.refresh_structure()

    #DEZE CRASHT

  #  product.print_graph()
  #  product.print_bonds()
  #  Drawer(product)

    x = product.split_disconnected_structures()

    print(x)
    for mol in x:
        mol.print_graph()
        Drawer(mol)


    #print(product.graph)

    return product

def activate_AMP(AMP):
    for atom in AMP.graph:
        if atom.type == 'P':
            for neighbour in atom.neighbours:
                if neighbour.type == 'O':
                    for nbour in neighbour.neighbours:
                        if nbour.type == 'H':
                            AMP.break_bond_between_atoms(atom, neighbour)
                            AMP.remove_atom(neighbour)
                            AMP.remove_atom(nbour)
                            #print(AMP.graph)
                            return AMP


if __name__ == "__main__":
    from graph_to_smiles import *
    ATP_string = "C1=NC(=C2C(=N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N"
    Serine_string = "C([C@@H](C(=O)O)N)O"
    AMP_string = 'C1=NC(=C2C(=N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)O)O)O)N'

    ATP_smiles = Smiles(ATP_string)
    Serine_smiles = Smiles(Serine_string)
    AMP_smiles = Smiles(AMP_string)

    ATP = ATP_smiles.smiles_to_structure()
    Serine = Serine_smiles.smiles_to_structure()
    AMP = AMP_smiles.smiles_to_structure()

    #GraphToSmiles(ATP)
    #Drawer(ATP)
    #Drawer(Serine)
    #print(ATP.graph)
    #print(Serine.graph)
    #remove_PPi(ATP)

    #prepare_oxygen(Serine)
    #print(AMP.graph)
    product = do_adenylation(ATP, Serine)
    #Drawer(product)
    #print(product)


  #  products = product.split_disconnected_structures()
   # for product in products:
  #      product.print_graph()
    #print(products)
    #product.print_graph()
    #products[0].print_graph()
    #products[1].print_graph()
    #Drawer(product)
    #Drawer(products[0])

