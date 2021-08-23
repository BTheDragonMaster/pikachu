from pikachu import *
from combine_structure import *
from drawing import *

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
    AMP = copy.deepcopy(ATP)
    AMP.break_bond_between_atoms(active_P, active_Oxygen)
    to_remove = []
    for atom in AMP.graph:
        #print(atom.nr)
        if atom.nr > 22:
            #print(atom)
            to_remove.append(atom)
    #print(to_remove)
    #for atom in to_remove:
        #for neighbour in atom.neighbours:
            #print(atom, neighbour)
            #AMP.break_any_bond(atom, neighbour)

        #AMP.remove_atom(atom)
    #print(AMP.graph)
    return AMP

def prepare_oxygen(AA):
    #print(AA.graph)
    hydroxyl_groups = find_group(AA, HYDROXYL)
    #print(hydroxyl_groups)
    to_remove = []
    for atom in hydroxyl_groups:
        #print('atom is', atom)
        #print(atom.neighbours)
        for neighbour in atom.neighbours:
            #print('neighbour is', neighbour)
            for nbour in neighbour.neighbours:
                #print('nbour is', nbour)
                if nbour.type == 'O' and nbour != atom :
                    #print('REMOVE THIS')
                    to_remove.append(atom)
    #print(hydroxyl_groups)
    #print(to_remove)
    for A in to_remove:
        if A in hydroxyl_groups:
            hydroxyl_groups.remove(A)
    #print(hydroxyl_groups)
    oxygen = hydroxyl_groups[0]
    for neighbour in oxygen.neighbours:
        if neighbour.type == 'H':
            hydrogen = neighbour
            AA.break_bond_between_atoms(oxygen, hydrogen)
            AA.remove_atom(hydrogen)
            break
    #print(oxygen)

    return AA, oxygen


def do_adenylation(ATP, AA, AMP_structure):
    AMP = remove_PPi(ATP)

    #AMP = activate_AMP(AMP_structure)

    #print(AMP.graph)

    active_AA, active_O = prepare_oxygen(AA)
    #print(active_O)
    structures = [AMP, active_AA]
    correct_structures_index(structures)
    #print(AMP.graph)
    #print(active_AA.graph)

    reactive_P = find_group(AMP, AMP_PHOSPHOR)
    active_P = reactive_P[0]

    product = combine_graphs(structures)
    #print(product.graph)

    #print(active_O)
    #print(reactive_P)
    next_nr = product.find_next_bond_nr()
    product.make_bond(active_O, active_P, next_nr)
    product.set_atom_neighbours()
    product.make_bond_lookup()

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
    ATP_string = "C1=NC(=C2C(=N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N"
    Serine_string = "C([C@@H](C(=O)O)N)O"
    AMP_string = 'C1=NC(=C2C(=N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)O)O)O)N'

    ATP_smiles = Smiles(ATP_string)
    Serine_smiles = Smiles(Serine_string)
    AMP_smiles = Smiles(AMP_string)

    ATP = ATP_smiles.smiles_to_structure()
    Serine = Serine_smiles.smiles_to_structure()
    AMP = AMP_smiles.smiles_to_structure()

    #print(ATP.graph)
    #print(Serine.graph)
    #remove_PPi(ATP)

    #prepare_oxygen(Serine)

    products = do_adenylation(ATP, Serine, AMP).split_disconnected_structures()
    Drawer(products[1].kekulise())
    for product in products:
        product.print_graph()
        GraphToSmiles(product)
        product = product.kekulise()
        Drawer(product)



