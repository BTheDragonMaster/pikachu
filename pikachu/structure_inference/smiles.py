#!/usr/bin/env python

from pikachu.structure_inference.structures import Structure
from pikachu.structure_inference.atom import Atom
from pikachu.structure_inference.bond import Bond

from typing import *
from dataclasses import dataclass
import sys
"""
Smiles class
"""

sys.setrecursionlimit(100000)

def calc_charge(sign, value):
    if sign == '+':
        return value
    elif sign == '-':
        return value * -1
    else:
        raise Exception("Wrong character to indicate charge!")

def parse_explicit(component):
    skip = False
    informative = component[1:-1]
    
    charges = []
    hydrogen = None
    numbers = []
    element = []
    chirals = []

    hydrogens = 0
    chiral = None
    charge = 0


    for i, character in enumerate(informative):
        if skip:
            skip = False
            continue
        if character.isupper():
            
            if character == 'H':
                hydrogen = i
            else:
                try:
                    if informative[i + 1].islower():
                        element.append(i)
                        element.append(i + 1)
                        skip = True
                    else:
                        element.append(i)
                except IndexError:
                    element.append(i)
        elif character.islower():
            element.append(i)
        elif character.isdigit():
            numbers.append(i)
        elif character == '+' or character == '-':
            charges.append(i)
        elif character == '@':
            chirals.append(i)

    element = ''.join([informative[x] for x in element])

    #Parsing the charge

    if len(charges) == 1:
        index = charges[0]
        charge_type = informative[index]
        try:
            if (index + 1) in numbers:
                charge_value = int(informative[index + 1])
                charge = calc_charge(charge_type, charge_value)
            else:
                charge = calc_charge(charge_type, 1)
                
        except IndexError:
            charge = calc_charge(charge_type, 1)

    elif len(charges) == 0:
        charge = 0
    else:
        charge_type = informative[charges[0]]
        charge_value = len(charges)
        charge = calc_charge(charge_type, charge_value)

    #Parsing an explicit hydrogen

    if hydrogen != None:
        try:
            if (hydrogen + 1) in numbers:
                hydrogens = int(informative[hydrogen + 1])
            else:
                hydrogens = 1
        except IndexError:
            hydrogens = 1
    else:
        hydrogens = 0

    #Parsing chirality

    if len(chirals) == 1:
        chiral = 'counterclockwise'
    elif len(chirals) == 2:
        chiral = 'clockwise'
    else:
        chiral = None

    return element, chiral, charge, hydrogens


def make_character_dict() -> Dict[str, str]:
    """Create dict of {character: label, ->} to label smiles characters
    """
    character_dict = {}
    atoms = ["C", "O", "N", "S", "B", "P", "F", "I", "c", "n", "o", '*',
             'Cl', 'Br', 'p', 'b', 'p', 's']
    cyclic = list(range(1,100))
    

    for atom in atoms:
        character_dict[atom] = "atom"
    for number in cyclic:
        character_dict[str(number)] = "cyclic"

    character_dict["="] = "double_bond"
    character_dict["("] = "branch_start"
    character_dict[")"] = "branch_end"
    character_dict['\\'] = 'chiral_double_bond'
    character_dict['/'] = 'chiral_double_bond'
    character_dict['#'] = 'triple_bond'
    character_dict['$'] = 'quadruple_bond'
    character_dict[':'] = 'aromatic'
    character_dict['.'] = 'split'
    character_dict['-'] = 'single_bond'
    character_dict[':'] = 'aromatic_bond'
    
    return character_dict

class Smiles():
    character_dict = make_character_dict()
    two_atom_dict = {'B': {'r'}, 'C': {'l'}}
    
            
    def __init__(self, string: str) -> None:
        self.smiles = string
        self.get_components()
        

    def get_components(self):
        self.components = []

        skip = False
        double_digits = False
        square_brackets = False
        
        for i, character in enumerate(self.smiles):
            if skip:
                skip = False
            elif square_brackets:
                component += character
                if character == ']':
                    square_brackets = False
                    self.components.append(component)
                    component = ''
            elif double_digits:
                try:
                    next_character = self.smiles[i + 1]
                    
                    if next_character not in {'0', '1', '2', '3', '4',
                                              '5', '6', '7', '8', '9'}:
                        double_digits = False
                        self.components.append(component + character)
                        component = ''
                    else:
                        component += character
                except IndexError:
                    self.components.append(component + character)
            else:
                
                if character in self.two_atom_dict:
                    try:
                        next_character = self.smiles[i + 1]
                        if next_character in self.two_atom_dict[character]:
                            self.components.append(character + next_character)
                            skip = True
                        else:
                            self.components.append(character)
                            
                        
                    except IndexError:
                        self.components.append(character)
                elif character == '[':
                    square_brackets = True
                    component = character
                elif character == '%':
                    double_digits = True
                    component = ''
                else:
                    self.components.append(character)
            

    def smiles_to_structure(self):

        structure = Structure()

        #Keeps track of the layer of brackets the atom is in

        branch_level = 0

        #Keeps track of which cyclic atoms have been encountered
        cyclic_dict = {}

        #Keeps track of the last atom encountered per branch level
        last_atoms_dict = {0: None}

        #Keeps track of the last double bond encountered, and the last bond chirality marker associated with it
        last_double_chiral_dict = {0: None}
        last_double_bond_dict = {0: None}
        double_chiral_active_dict = {0: False}
        chiral_atoms_1_double = {0: None}


        #Keeps track of the nature of the bond
        bond_type = 'single'

        #Keeps track of chiral centres
        chiral_dict = {}
        chiral_double_bond_dict = {}
        current_chiral_atom = None
        

        explicit = False

        atom_nr = -1
        bond_nr = -1
        
        for i, component in enumerate(self.components):
            if component[0] == '[':
                label = 'atom'
                explicit = True
            else:
                label = self.character_dict[component]

            if label == 'split':
                #Starts disconnected structure; set everything back to default
                branch_level = 0
                cyclic_dict = {}
                last_atoms_dict = {0: None}
                double = False
                triple = False
                chiral_dict = {}

            elif label == "atom":

                if not explicit:
                    element = component
                    chiral = None
                    charge = 0
                    hydrogens = 0
                else:
                    element, chiral, charge, hydrogens = parse_explicit(component)
                    explicit = False

                #Determine atom
                
                atom_nr += 1

                #Turn atom into an atom object
                atom_2 = Atom(element, atom_nr, chiral, charge)

                #Keep track of atoms surrounding chiral double bonds
                try:
                    if double_chiral_active_dict[branch_level]:
                        last_double_bond, last_double_bond_index = last_double_bond_dict[branch_level]

                        chiral_double_bond_dict[last_double_bond]['atom 2'] = atom_2
                        double_chiral_active_dict[branch_level] = False
                except KeyError:
                    pass
                    
                for i in range(hydrogens):
                    atom_nr += 1
                    bond_nr += 1
                    hydrogen = Atom('H', atom_nr, None, 0)
                    structure.add_bond(atom_2, hydrogen, 'single', bond_nr)

                
                        
                    
                atom_1 = self.get_last_atom(branch_level, last_atoms_dict)

                previous_atom_branch_level = branch_level

                #Go back through the branches to identify the atom that the
                #current atom is attached to; call this atom_1
                
                while previous_atom_branch_level > 0 and not atom_1:
                    previous_atom_branch_level -= 1
                    atom_1 = last_atoms_dict[previous_atom_branch_level]

                if atom_1:
                    bond_nr += 1
                    
                    if atom_1.type.islower() and atom_2.type.islower():
                        if not bond_type == 'explicit_single':
                            bond_type = 'aromatic'
                        else:
                            bond_type = 'single'
                        
                    structure.add_bond(atom_1, atom_2, bond_type, bond_nr)
                    if bond_type == 'double' or bond_type == 'triple':
                        last_double_bond_dict[branch_level] = (bond_nr, i)
                    

                    if atom_1.chiral:
                        #Add current atom to list of residues
                        chiral_dict[atom_1].append(atom_2)                          
                            
                    if atom_2.chiral:
                        chiral_dict[atom_2] = [atom_1]
                        
                        if hydrogens == 1:
                            chiral_dict[atom_2].append(hydrogen)

                    bond_type = 'single'

                            
                            

                else:
                    
                    #This happens only if atom_2 is completely disconnected
                    #of any atoms already in the graph
                    structure.add_atom(atom_2)
                    if atom_2.chiral:
                        chiral_dict[atom_2] = []
                        if hydrogens == 1:
                            chiral_dict[atom_2].append(hydrogen)

                #Set atom_2 as the last atom in the current branch level
                self.track_last_atoms_per_branch(atom_2, branch_level,
                                                 last_atoms_dict)

                
            elif label == "single_bond":
                bond_type = 'explicit_single'

            elif label == 'aromatic_bond':
                bond_type = 'aromatic'
                
            elif label == "double_bond":
                bond_type = 'double'
                
            elif label == "triple_bond":
                bond_type = 'triple'
                
            elif label == "quadruple_bond":
                bond_type = 'quadruple'

            #If there are brackets: go up or down a branch level
            elif label == "branch_start":
                branch_level += 1
                last_atoms_dict[branch_level] = None

            elif label == "branch_end":
                last_atoms_dict[branch_level] = None
                branch_level -= 1

            elif label == "cyclic":

                cycle_nr = int(component)
                atom = self.get_last_atom(branch_level, last_atoms_dict)
                

                #If we encounter a number that hasn't been encountered before
                
                if self.new_cycle(cyclic_dict, cycle_nr):
                    self.start_cycle(cycle_nr, atom, cyclic_dict)
                    
                    if atom in chiral_dict:
                        self.add_cycle_placeholder(chiral_dict, atom, cycle_nr)

                #Otherwise look up the atom that the cycle closes on
                else:
                    bond_nr += 1
                    
                        
                    atom_1, atom_2 = self.end_cycle(cycle_nr, atom,
                                                    cyclic_dict)
                    
                    if atom_1.type.islower() and atom_2.type.islower():
                        if not bond_type == 'explicit_single':
                            bond_type = 'aromatic'
                        else:
                            bond_type = 'single'
                        
                    structure.add_bond(atom_1, atom_2, bond_type, bond_nr)
                    
                    
                    if atom_1 in chiral_dict:
                        self.replace_cycle_placeholder(chiral_dict, atom_1, atom_2, cycle_nr)
                    if atom_2 in chiral_dict:
                        chiral_dict[atom_2].append(atom_1)

                    bond_type = 'single'
            elif label == 'chiral_double_bond':
                if branch_level in last_double_chiral_dict and last_double_chiral_dict[branch_level]:
                    last_chiral_sign, last_chiral_index = last_double_chiral_dict[branch_level]
                    last_double_bond, last_double_bond_index = last_double_bond_dict[branch_level]
                    if last_chiral_sign and (last_chiral_index < last_double_bond_index):
                        if last_chiral_sign == component:
                            orientation = 'trans'
                        else:
                            orientation = 'cis'
                            
                        chiral_double_bond_dict[last_double_bond] = {}
                        chiral_double_bond_dict[last_double_bond]['orientation'] = orientation
                        chiral_double_bond_dict[last_double_bond]['atom 1'] = chiral_atoms_1_double[branch_level]

                        double_chiral_active_dict[branch_level] = True


                last_double_chiral_dict[branch_level] = (component, i)
                chiral_atoms_1_double[branch_level] = last_atoms_dict[branch_level]


        for atom in chiral_dict:
            order = chiral_dict[atom]

            #Save chirality in the atom object
            current_chirality = atom.chiral
            new_chirality = self.determine_chirality(order, current_chirality)
            
            
            atom.chiral = new_chirality

        structure.refine_structure()

        for double_bond in chiral_double_bond_dict:
            bond = structure.bonds[double_bond]
            orientation = chiral_double_bond_dict[double_bond]['orientation']
            atom_1 = chiral_double_bond_dict[double_bond]['atom 1']
            if 'atom 2' in chiral_double_bond_dict[double_bond]:
                atom_2 = chiral_double_bond_dict[double_bond]['atom 2']
                bond.set_double_bond_chirality(atom_1, atom_2, orientation)

        

        return structure
    

    def kekulize(self) -> str:
        """Return kekulized smiles from smiles.

        Input:
        smiles: str, smiles string

        Output:
        K_smiles: str, smiles string, kekulized version of smiles
        """
        
        return K_smiles

 

    def draw_smiles(self, img_dir: str, highlight_map: Dict[int, List[float]], bond_list: List[int], bond_colour: List[float], ID: str) -> None:
        """Create image from smiles string

        Input:
        img_dir: str, directory of output image
        highlight_map: dict of {atom_nr: [float_1, float_2, float_3], ->}, with
            atom_nr int and the combination of floats representing an RGB colour
        bond_list: list of int, with each int a bond index of the molecule
        bond_colour: list of [float_1, float_2, float_3], representing an
            RGB colour
        ID: str, name of the molecule
        """

        mol = Chem.MolFromSmiles(self.smiles)
        mol = Chem.Mol(mol.ToBinary())

        if not mol.GetNumConformers():
            rdDepictor.Compute2DCoords(mol)

        colormap = dict(zip(bond_list, [bond_colour] * len(bond_list)))
        
        print(colormap)        
        
        size = (500, 500)
 #       drawer = rdMolDraw2D.MolDraw2DSVG(size[0], size[1])
#        drawer.DrawMolecule(mol, highlightAtoms = list(highlight_map.keys()),
#                            highlightAtomColors = highlight_map,
#                            highlightBonds = bond_list,
#                            highlightBondColors=colormap)
#        drawer.FinishDrawing()
        
 #       img = Chem.Draw.MolToImage(mol, size=size, kekulize=True,
#                                   wedge_bonds = True, fitImage = True,
#                                   highlightMap = highlight_map,
#                                   highlightBonds = bond_list,
#                                   highlightColor = bond_colour)

        img = Chem.Draw.MolToFile(mol, "%s/%s.svg" % (img_dir, ID),
                                   size=size, kekulize=True,
                                   wedge_bonds = True, fitImage = True,
                                   highlightMap = highlight_map,
                                   highlightColor = bond_colour,
                                   highlightBonds = bond_list)

        
 #       img.save("%s/%s.png" % (img_dir, ID))
#        svg = drawer.GetDrawingText()
 #       with open("%s/%s.svg" % (img_dir, ID), 'w') as f:
#            f.write(svg)
        

    #=========================================================================
    #Auxillary functions to smiles_to_structure
    #=========================================================================
        

    def add_chiral_atom(self, chiral_dict: Dict['Atom', Dict[str, Any]], last_atom: 'Atom', current_atom: 'Atom') -> None:
        """Place current_atom in one of the four bond slots of last_atom

        Input:
        chiral_dict: dict of {atom: {'direction': direction, 'order':
            [atom_1, atom_2, atom_3, atom_4]}}, with atom Atom Object,
            direction int, and order a list of Atom Object or int or None
        last_atom: Atom Object, chiral atom
        current_atom: Atom Object, neighbour of chiral atom

        """


        chiral_dict[last_atom].append(current_atom)

    def add_cycle_placeholder(self, chiral_dict, atom, cycle_nr):
        chiral_dict[atom].append(cycle_nr)

    def replace_cycle_placeholder(self, chiral_dict, chiral_atom, current_atom, cycle_nr):
        for i, atom in enumerate(chiral_dict[chiral_atom]):
            if type(atom) == int:
                if atom == cycle_nr:
                    chiral_dict[chiral_atom][i] = current_atom
        

    def get_chiral_permutations(self, order):
        permutations = [tuple(order)]
        permutations.append((order[0], order[3], order[1], order[2]))
        permutations.append((order[0], order[2], order[3], order[1]))
        permutations.append((order[1], order[0], order[3], order[2]))
        permutations.append((order[1], order[2], order[0], order[3]))
        permutations.append((order[1], order[3], order[2], order[0]))
        permutations.append((order[2], order[0], order[1], order[3]))
        permutations.append((order[2], order[3], order[0], order[1]))
        permutations.append((order[2], order[1], order[3], order[0]))
        permutations.append((order[3], order[0], order[2], order[1]))
        permutations.append((order[3], order[1], order[0], order[2]))
        permutations.append((order[3], order[2], order[1], order[0]))

        return permutations
        
        
        
    def determine_chirality(self, order, chirality):
        original_order = tuple(order)
        chiral_permutations = self.get_chiral_permutations(original_order)
        order.sort(key = lambda x: x.nr)
        new_order = tuple(order)
        if new_order in chiral_permutations:
            if chirality == 'counterclockwise':
                new_chirality = 'counterclockwise'
            else:
                new_chirality = 'clockwise'
        else:
            if chirality == 'counterclockwise':
                new_chirality = 'clockwise'
            else:
                new_chirality = 'counterclockwise'

        return new_chirality



    def new_cycle(self, cyclic_dict: Dict[int, 'Atom'], cycle_nr: int) -> bool:
        """Return bool, True if a new cycle is recorded, False if not

        Input:
        cyclic_dict: dict of {cycle_nr: atom, ->}, with cycle_nr int and atom
            an Atom Object
        cycle_nr: int, nr of the current cycle

        Output:
        bool: True if an atom with cycle_nr at position 0 does not yet exist in
            cyclic_atoms, False if it does
        """
        if cycle_nr in cyclic_dict:
            return False
        else:
            return True

    def start_cycle(self, cycle_nr: int, atom: 'Atom', cyclic_dict: Dict[int, 'Atom']) -> None:
        """Add a new atom and corresponding cycle number to cyclic dict

        Input:
        cycle_nr: int, nr of the current cycle
        atom: Atom Object
        cyclic_dict: dict of {cycle_nr: atom, ->}, with cycle_nr int and atom
            an Atom Object

        """
        cyclic_dict[cycle_nr] = atom

    def end_cycle(self, cycle_nr: int, atom: 'Atom', cyclic_dict: Dict[int, 'Atom']) -> Tuple['Atom']:
        """Return pair of atoms that close a cycle

        Input:
        cycle_nr: int, nr of the current cycle
        atom: Atom Object
        cyclic_dict: dict of {cycle_nr: atom, ->}, with cycle_nr int and atom
            an Atom Object

        Output:
        atom_pair: tuple of two atoms, with each atom an Atom Object
        
        """
        atom_old = cyclic_dict[cycle_nr]
        atom_pair = (atom_old, atom)
        
        del cyclic_dict[cycle_nr]

        #Implement atom pair here?
        
        return atom_pair

    def track_last_atoms_per_branch(self, new_atom: 'Atom', current_level: int, last_atoms_dict: Dict[int, 'Atom']) -> None:
        """Update which atom was the last to occur at the current branch level

        Input:
        new_atom: Atom Object
        current_level: int, current branch level, level to which new_atom is
            added
        last_atoms_dict: dict of {int: atom, ->}, with int representing a branch
            level, and atom representing the last atom that occurs in that
            branch. 
        """
        last_atoms_dict[current_level] = new_atom

    def get_last_atom(self, current_level: int, last_atoms_dict: Dict[int, 'Atom']) -> 'Atom':
        """Return the last atom in the current branch level

        Input:
        current_level: int, current branch level, level from which the last atom
            is to be extracted
        last_atoms_dict: dict of {int: atom, ->}, with int representing a branch
            level, and atom representing the last atom that occurs in that
            branch.

        Output:
        last_atom: Atom Object, last atom that was encountered in the
            current_level branch
        """
        last_atom = last_atoms_dict[current_level]
        return last_atom
        
    def check_character(self, character: str) -> str:
        """Return the label for a character to determine its type

        Input:
        character: str, substring of smiles string of length 1

        Output:
        label: str, descriptor of character
        """

        try:
            #Consider making character_dict a class attribute defined before
            #__init__
            label = character_dict[character]
        except KeyError:
            label = "Unknown"
        return label

    def update_structure(self, atom_1: 'Atom', atom_2: 'Atom', structure_graph: Dict['Atom', List['Atom']], bond_type: str) -> None:
        """Add an atom to the structure graph

        Input:
        atom_1: Atom Object
        atom_2: Atom Object
        structure_graph: dict of {atom: [atom, ->], ->}, with each atom an
        Atom Object
        """
        if atom_1 in structure_graph:
            structure_graph[atom_1] += [atom_2]
        else:
            structure_graph[atom_1] = [atom_2]

        if atom_2 in structure_graph:
            structure_graph[atom_2] += [atom_1]
        else:
            structure_graph[atom_2] = [atom_1]

def main():
    sys.setrecursionlimit(100000)
    string = "CC(C)(O1)C[C@@H](O)[C@@]1(O2)[C@@H](C)[C@@H]3CC=C4[C@]3(C2)C(=O)C[C@H]5[C@H]4CC[C@@H](C6)[C@]5(C)Cc(n7)c6nc(C[C@@]89(C))c7C[C@@H]8CC[C@@H]%10[C@@H]9C[C@@H](O)[C@@]%11(C)C%10=C[C@H](O%12)[C@]%11(O)[C@H](C)[C@]%12(O%13)[C@H](O)C[C@@]%13(C)CO"
 #   string = '[Co+2].[Fe+++].[O-].[OH2].[OH-].[CH4]'
 #   string = "c1ccccc1-c2ccccc2"
    string = "CCCCCCCCCC(=O)NC1C(O)C(O)C(CO)OC1Oc2c3Oc4ccc(CC5NC(=O)C(N)c6ccc(O)c(Oc7cc(O)cc(c7)C(NC5=O)C(=O)NC8C(=O)NC9C(=O)NC(C(OC%10OC(CO)C(O)C(O)C%10NC(C)=O)c%11ccc(Oc2cc8c3)c(Cl)c%11)C(=O)NC(C(O)=O)c%12cc(O)cc(OC%13OC(CO)C(O)C(O)C%13O)c%12c%14cc9ccc%14O)c6)cc4Cl"
    string = "CCCCCCCCCC(=O)N[C@@H](CC1=CNC2=CC=CC=C21)C(=O)N[C@H](CC(=O)N)C(=O)N[C@@H](CC(=O)O)C(=O)N[C@H]3[C@H](OC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@H](NC(=O)CNC(=O)[C@@H](NC(=O)[C@H](NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)CNC3=O)CCCN)CC(=O)O)C)CC(=O)O)CO)[C@H](C)CC(=O)O)CC(=O)C4=CC=CC=C4N)C"
#    string = "CC1CCC/C(C)=C1/C=C/C(C)=C/C=C/C(C)=C/C=C/C=C(C)/C=C/C=C(C)/C=C/C2=C(C)/CCCC2(C)C"
#    string = 'CC'
#    string = 'O=P([O-])([O-])[O-]'
 #   string = 'C#N'
 #   string = 'FS(F)(F)F'
#    string_2 = 'c1ccccc1'
#    string = 'C(C[C@@H](C(=O))N)CN'
 #   string_2 = 'C(C[C@@H](C(=O))N)CN'
#    string = 'C([C@@H](C(=O)O)N)'
#    string_2 = 'CC[C@@H](C)[C@@H](C(=O))N'
    string_2 = 'CC(C)[C@H](C(=O))N'
#    string_2 = 'C(=O)NC'
 #   string = 'FS(F)(F)(F)'
 #   string = 'C/C=C/C'
 #   string_2 = 'C/C=C/C'
    smiles_1 = Smiles(string)
    smiles_2 = Smiles(string_2)
    structure_1 = smiles_1.smiles_to_structure()
#    pprint(structure_1.bonds)
#    pprint(structure_1.graph)
#    print(structure_1.find_peptide_bonds())
    
    structure_2 = smiles_2.smiles_to_structure()

    pprint(structure_1.substructure_matching(structure_2, False))
    pprint(structure_1.find_aromatic_systems())
#    print('reverse')
 #   print(structure_2.substructure_matching(structure_1))
    
 ##   structure_1.print_graph()
#    structure_2.print_bonds()

 #   structure_1.bonds[1].break_bond()
#    for carbon in structure_1.graph:
#        print(carbon)
#        for shell in carbon.shells:
#            for orbital_set in carbon.shells[shell].orbital_sets:
#                for orbital in carbon.shells[shell].orbital_sets[orbital_set].orbitals:
#                    print(orbital)
#                    print(orbital.electrons)
#    aromatics = structure.find_aromatic_graphs()
#    for aromatic in aromatics:
#        aromatic.print_graph()
    


if __name__ == "__main__":
    main()
    


