from smiles import *
from structures import *
import time

class BondDefiner():
    def __init__(self, name, smiles, atom_nr_1, atom_nr_2):
        self.name = name
        self.smiles = Smiles(smiles)
        self.structure = self.smiles.smiles_to_structure()
        self.find_atoms(atom_nr_1, atom_nr_2)
        self.bond = self.structure.bond_lookup[self.atom_1][self.atom_2]
        
    def __repr__(self):
        return self.name
        

    def find_atoms(self, atom_nr_1, atom_nr_2):

        self.atom_1 = None
        self.atom_2 = None
        
        for atom in self.structure.graph:
            if atom.nr == atom_nr_1:
                self.atom_1 = atom
            elif atom.nr == atom_nr_2:
                self.atom_2 = atom

            if self.atom_1 and self.atom_2:
                break

        if not self.atom_1 or not self.atom_2:
            raise Exception("Can't find atoms adjacent to bond.")

PEPTIDEBOND = BondDefiner('peptide_bond', 'C(=O)NC', 0, 2)
ESTERCOCBOND = BondDefiner('ester_coc_bond', 'C(=O)OC', 2, 3)
BENZENEBENZENEBOND = BondDefiner('benzene_benzene_bond', 'c1ccccc1-c2ccccc2', 5, 6)


def do_hydrolysis(product, bond, heteroatom, other_atom):
    """Break indicated bond in indicated structure, adding an -H to heteroatom and an -OH to other_atom

    Auxillary function to 'hydrolyse'.
    """
    
    product.break_bond_nr(bond)

    oxygen = product.make_atom('O', [other_atom])
    hydrogen_1 = product.make_atom('H', [oxygen])
    hydrogen_2 = product.make_atom('H', [heteroatom])
    
    
def hydrolyse(bond, structure, oxygen_recipient = None):
    assert bond.type == 'single'

    hydrogen_recipient = None

    if oxygen_recipient:
        for atom in bond.neighbours:
            if atom != oxygen_recipient:
                hydrogen_recipient = atom

        do_hydrolysis(structure, bond, hydrogen_recipient, oxygen_recipient)

    else:
    
        heteroatoms = []
        carbons = []
        
        for atom in bond.neighbours:
            if atom.type != 'C':
                heteroatoms.append(atom)
            else:
                carbons.append(atom)

        if len(heteroatoms) == 0:
            do_hydrolysis(structure, bond, bond.neighbours[0], oxygen_recipient)

        elif len(heteroatoms) == 2:
            non_oxygens = []
            oxygens = []
            for atom in heteroatoms:
                if atom.type != 'O':
                    non_oxygens.append(atom)
                else:
                    oxygens.append(atom)

            if len(non_oxygens) == 2:
                do_hydrolysis(structure, bond, non_oxygens[0], non_oxygens[1])
                
            elif len(non_oxygens) == 1:
                do_hydrolysis(structure, bond, oxygens[0], non_oxygens[0])
            else:
                raise Exception(f"Not a valid bond between {heteroatoms[0]} and {heteroatoms[1]}")
                

        elif len(heteroatoms) == 1:
            do_hydrolysis(structure, bond, heteroatoms[0], carbons[0])


def do_condensation(structure, hydrogen_donor, hydrogen, oxygen_donor, oxygen):

    structure.break_bond_atoms(hydrogen_donor, hydrogen)
    structure.break_bond_atoms(oxygen_donor, oxygen)
    
    for atom in oxygen.neighbours:
        if atom.type == 'H':
            hydrogen_2 = atom

    structure.remove_atom(hydrogen)
    structure.remove_atom(hydrogen_2)
    structure.remove_atom(oxygen)

    bond_nr = structure.find_next_bond_nr()
    structure.make_bond(hydrogen_donor, oxygen_donor, bond_nr)

def condensation(structure, atom_1, atom_2):
    atom_1_hydrogens = []
    atom_1_oxygens = []

    atom_2_hydrogens = []
    atom_2_oxygens = []

    products = []

    for neighbour in atom_1.neighbours:
        if neighbour.type == 'H':
            atom_1_hydrogens.append(neighbour)
        elif neighbour.type == 'O':
            atom_1_oxygens.append(neighbour)

    for neighbour in atom_2.neighbours:
        if neighbour.type == 'H':
            atom_2_hydrogens.append(neighbour)
        elif neighbour.type == 'O':
            atom_2_oxygens.append(neighbour)

    atom_1_reactive_oxygen = None
    atom_2_reactive_oxygen = None

    for oxygen in atom_1_oxygens:
        for neighbour in oxygen.neighbours:
            if neighbour.type == 'H':
                atom_1_reactive_oxygen = oxygen
                break

    for oxygen in atom_2_oxygens:
        for neighbour in oxygen.neighbours:
            if neighbour.type == 'H':
                atom_2_reactive_oxygen = oxygen
                break

    if not atom_1_hydrogens and not atom_2_hydrogens:
        raise Exception("Can't condensate these two atoms! No hydrogens!")

    if not atom_2_hydrogens and not atom_2_reactive_oxygen:
        raise Exception(f"Can't condensate these two atoms! Atom {atom_2} cannot react!")

    if not atom_1_hydrogens and not atom_1_reactive_oxygen:
        raise Exception(f"Can't condensate these two atoms! Atom {atom_1} cannot react!")
        
    if not atom_1_reactive_oxygen and not atom_2_reactive_oxygen:
        raise Exception("Can't condensate these two atoms! No reactive oxygens!")

    if atom_1_reactive_oxygen:
        if atom_2_hydrogens:
            product = copy.deepcopy(structure)
            do_condensation(product, atom_2, atom_2_hydrogens[0], atom_1, atom_1_reactive_oxygen)
            products.append(product)

    if atom_2_reactive_oxygen:
        if atom_1_hydrogens:
            product = copy.deepcopy(structure)
            do_condensation(product, atom_1, atom_1_hydrogens[0], atom_2, atom_2_reactive_oxygen)
            products.append(product)


    return products
                
                    
                    
            
        
            
            
        

class BondStructures():
    peptide_bond = {'smiles': 'C(=O)NC',
                    'structure': Smiles('C(=O)NC').smiles_to_structure(),
                    'atom_1': 'C_0',
                    'atom_2': 'N_2'}
    

BONDSTRUCTURES = BondStructures()

def find_bonds(structure, bond_type):
    locations = structure.substructure_matching(bond_type.structure)
    bonds = []
    for match in locations:
        atom_1 = match[bond_type.atom_1]
        atom_2 = match[bond_type.atom_2]
        bond = structure.bond_lookup[atom_1][atom_2]
        bonds.append(bond)

    return bonds


    

if __name__ == "__main__":
    string = "CCCCCCCCCC(=O)NC1C(O)C(O)C(CO)OC1Oc2c3Oc4ccc(CC5NC(=O)C(N)c6ccc(O)c(Oc7cc(O)cc(c7)C(NC5=O)C(=O)NC8C(=O)NC9C(=O)NC(C(OC%10OC(CO)C(O)C(O)C%10NC(C)=O)c%11ccc(Oc2cc8c3)c(Cl)c%11)C(=O)NC(C(O)=O)c%12cc(O)cc(OC%13OC(CO)C(O)C(O)C%13O)c%12c%14cc9ccc%14O)c6)cc4Cl"
    string = "CCCCCCCCCC(=O)N[C@@H](CC1=CNC2=CC=CC=C21)C(=O)N[C@@H](CC(=O)N)C(=O)N[C@@H](CC(=O)O)C(=O)N[C@H]3[C@H](OC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@H](NC(=O)CNC(=O)[C@@H](NC(=O)[C@H](NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)CNC3=O)CCCN)CC(=O)O)C)CC(=O)O)CO)[C@H](C)CC(=O)O)CC(=O)C4=CC=CC=C4N)C"
    smiles_1 = Smiles(string)
    structure_1 = smiles_1.smiles_to_structure()
    

    peptide_bonds = find_bonds(structure_1, PEPTIDEBOND)
    coc_ester_bonds = find_bonds(structure_1, ESTERCOCBOND)
    

    for carbon in structure_1.graph:
        print(carbon)
        for shell in carbon.shells:
            for orbital_set in carbon.shells[shell].orbital_sets:
                for orbital in carbon.shells[shell].orbital_sets[orbital_set].orbitals:
                    print(orbital)
                    print(orbital.electrons)

    print(peptide_bonds)
    bonds = coc_ester_bonds + peptide_bonds
    product = copy.deepcopy(structure_1)

    for bond in bonds:
        hydrolyse(bond, product)


    structures = product.split_disconnected_structures()
    print(len(structures))
    for structure in structures:
        pprint(structure.graph)

        for carbon in structure.graph:
            print(carbon)
            for shell in carbon.shells:
                for orbital_set in carbon.shells[shell].orbital_sets:
                    for orbital in carbon.shells[shell].orbital_sets[orbital_set].orbitals:
                        if orbital.electron_nr != 2:
                            print(orbital.electrons)
 #                       print(orbital)
#                        print(orbital.electrons)
    
        
