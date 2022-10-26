from hashlib import sha256
from pikachu.chem.atom import Atom
from pikachu.chem.rings.find_cycles import find_cycles, in_cycle

from pprint import pprint


class InChI:
    def __init__(self, inchi):
        self.inchi = inchi
        self.layers = {'version': '',
                       'formula': '',
                       'connectivity': '',
                       'hydrogens': '',
                       'charge': '',
                       'protonation': '',
                       'double bond stereochemistry': '',
                       'tetrahedral stereochemistry': '',
                       'allene stereochemistry': '',
                       'stereochemistry information': '',
                       'isotopic layer': ''}

        self.graph = {}
        self.atoms = []
        self.hydrogens = []
        self.atom_to_hydrogens = {}
        self.cycles = None
        self.atom_to_cyclic = {}

    def set_cycles(self):
        self.cycles = find_cycles(self.graph)
        for atom in self.atoms:
            if in_cycle(self.cycles, atom):
                self.atom_to_cyclic[atom] = True
            else:
                self.atom_to_cyclic[atom] = False

    def inchi_to_structure(self):
        self.get_layers()
        self.parse_formula_layer()
        self.parse_connectivity_layer()
        self.set_cycles()

        print(self.atoms)
        print(self.hydrogens)
        print(self.graph)
        self.parse_hydrogen_layer()
        print(self.atom_to_cyclic)

    def get_layers(self):
        layers = self.inchi.split(r'/')
        version = layers[0].split('InChI=')[-1]
        version = version.split('S')[0]
        self.layers['version'] = version
        self.layers['formula'] = layers[1]

        fixed_h_layer = False

        for layer in layers[2:]:
            prefix = layer[0]

            if prefix == 'c':
                self.layers['connectivity'] = layer[1:]
            elif prefix == 'h':
                self.layers['hydrogens'] = layer[1:]
        print(self.layers)

    def parse_formula_layer(self):

        hydrogen_nr = 0
        atom_type_list = []

        atom = []
        digit = []
        formula_layer = self.layers['formula']
        for i, symbol in enumerate(formula_layer):
            previous_symbol = None
            if i != 0:
                previous_symbol = formula_layer[i - 1]

            # Start of a new atom type

            if symbol.isalpha() and symbol.isupper():
                print(previous_symbol)
                # Nothing to record yet if it is the first atom type of the list

                if previous_symbol:
                    atom_type = ''.join(atom)

                    # If there is more than one of the previous atom type

                    if previous_symbol.isdigit():
                        atom_number = int(''.join(digit))
                        if atom_type != 'H':
                            for i in range(atom_number):
                                atom_type_list.append(atom_type)
                        else:
                            hydrogen_nr = atom_number

                    # Happens if there is only one of the previous atom type

                    else:
                        if atom_type != 'H':
                            atom_type_list.append(atom_type)
                        else:
                            hydrogen_nr = 1

                    atom = []
                    digit = []

                atom.append(symbol)

            # For atom types like Cl an Br

            elif symbol.isalpha() and symbol.islower():
                atom.append(symbol)

            # Indicates number of atoms of a certain type.

            elif symbol.isdigit():
                digit.append(symbol)

            # If there is more than one of the previous atom type

        atom_type = ''.join(atom)

        # More than one of the last atom type

        if formula_layer[-1].isdigit():

            atom_number = int(''.join(digit))
            if atom_type != 'H':
                for i in range(atom_number):
                    atom_type_list.append(atom_type)
            else:
                hydrogen_nr = atom_number

        # Happens if there is only one of the last atom type

        else:
            if atom_type != 'H':
                atom_type_list.append(atom_type)
            else:
                hydrogen_nr = 1

        for i in range(hydrogen_nr):
            atom_type_list.append('H')

        for i, atom_type in enumerate(atom_type_list):
            atom = Atom(atom_type, i, None, 0, False)
            if atom_type != 'H':
                self.atoms.append(atom)
            else:
                self.hydrogens.append(atom)

    def get_connectivity_components(self):
        components = []
        digit = []
        for i, symbol in enumerate(self.layers['connectivity']):
            if symbol.isdigit():
                digit.append(symbol)
            else:
                digit_int = int(''.join(digit))
                components.append(digit_int)
                components.append(symbol)
                digit = []

        if self.layers['connectivity'][-1].isdigit():
            digit_int = int(''.join(digit))
            components.append(digit_int)

        return components

    def parse_connectivity_layer(self):

        for atom in self.atoms:
            self.graph[atom] = []

        components = self.get_connectivity_components()

        previous_component = None
        branch_to_previous_atom = {}

        current_branch = 0

        branch_to_previous_atom[0] = None

        for component in components:

            # Entering a new branch

            if component == '(':
                current_branch += 1
                if current_branch not in branch_to_previous_atom:
                    branch_to_previous_atom[current_branch] = None

            # Terminating a branch

            elif component == ')':
                current_branch -= 1

            elif type(component) == int:
                print(component)

                # Python uses 0-indexing, InChI doesn't
                atom = self.atoms[component - 1]

                # Previous atom is in the current branch

                if previous_component == '-' or previous_component == ')':
                    previous_atom = branch_to_previous_atom[current_branch]

                # Previous atom is in the previous branch
                elif previous_component == ',' or previous_component == '(':
                    previous_atom = branch_to_previous_atom[current_branch - 1]

                # This is the first atom; there is no previous atom
                else:
                    previous_atom = None

                if previous_atom:
                    print(previous_atom)
                    self.graph[atom].append(previous_atom)
                    self.graph[previous_atom].append(atom)

                branch_to_previous_atom[current_branch] = atom

            previous_component = component

    def get_hydrogen_components(self):
        hydrogen_layer = self.layers['hydrogens']
        components = []
        component = []
        between_brackets = False
        hydrogen_active = False
        for i, symbol in enumerate(hydrogen_layer):

            if symbol == '(':
                if component:
                    components.append(component)
                component = [symbol]
                between_brackets = True

            elif symbol == ')':
                component.append(symbol)
                between_brackets = False
                components.append(component)
                component = []

            elif symbol == 'H':
                hydrogen_active = True
                component.append(symbol)

            elif symbol.isdigit() or symbol == '-':
                component.append(symbol)

            elif symbol == ',':
                component.append(symbol)
                if hydrogen_active and not between_brackets:
                    components.append(component)
                    component = []

                hydrogen_active = False

        if component:
            components.append(component)

        for i, component in enumerate(components):
            components[i] = ''.join(component)

        return components

    def parse_hydrogen_component(self, component):

        hydrogen_position = None
        last_symbol = component[-1]

        if last_symbol == ')':
            component_type = 'shared'
            hydrogen_position = 1

        else:
            component_type = 'single'
            if last_symbol == 'H':
                hydrogen_position = len(component) - 1

            else:
                for i in range(len(component) - 1, -1, -1):
                    symbol = component[i]
                    if symbol == 'H':
                        hydrogen_position = i
                        break

        digit = []
        for symbol in component[hydrogen_position + 1:]:
            if symbol.isdigit():
                digit.append(symbol)
            else:
                break

        if digit:
            hydrogen_nr = int(''.join(digit))
        else:
            hydrogen_nr = 1

        atoms = []
        atom = []

        hydrogen_active = False
        range_active = False

        for symbol in component:
            if symbol.isdigit() and not hydrogen_active:
                atom.append(symbol)
            elif symbol == '-':

                current_atom = int(''.join(atom))
                atoms.append(current_atom)
                atom = []

                range_active = True
            elif symbol == 'H':

                hydrogen_active = True
                if atom:
                    current_atom = int(''.join(atom))
                    if not range_active:
                        atoms.append(current_atom)
                    else:
                        previous_atom = atoms[-1]
                        for i in range(previous_atom + 1, current_atom + 1):
                            atoms.append(i)
                    atom = []

                range_active = False

            elif symbol == ',':

                if atom and not hydrogen_active:
                    current_atom = int(''.join(atom))
                    if not range_active:
                        atoms.append(current_atom)
                    else:
                        previous_atom = atoms[-1]
                        for i in range(previous_atom + 1, current_atom + 1):
                            atoms.append(i)

                    atom = []

                hydrogen_active = False
                range_active = False

            else:
                if atom:
                    current_atom = int(''.join(atom))
                    if not range_active:
                        atoms.append(current_atom)
                    else:
                        previous_atom = atoms[-1]
                        for i in range(previous_atom + 1, current_atom + 1):
                            atoms.append(i)
                    atom = []

        if atom:
            print("Shouldn't happen!")
            current_atom = int(''.join(atom))
            atoms.append(current_atom)

        atom_objects = []
        for atom in atoms:
            for atom_object in self.atoms:
                if atom_object.nr + 1 == atom:
                    atom_objects.append(atom_object)

        return atom_objects, hydrogen_nr, component_type

    def add_hydrogens_to_graph(self):
        for atom, hydrogen_info in atom_to_hydrogens.items():
            if hydrogen_info:
                component_type = hydrogen_info[0]
                hydrogen_nr = hydrogen_info[1]
                if component_type == 'single':
                    for i in range(hydrogen_nr):
                        hydrogen = self.hydrogens.pop()
                        self.graph[atom].append(hydrogen)
                        self.graph[hydrogen] = [atom]

        for atom, hydrogen_info in atom_to_hydrogens.items():
            if hydrogen_info:
                component_type = hydrogen_info[0]
                hydrogen_nr = hydrogen_info[1]
                if component_type == 'shared':
                    atoms = hydrogen_info[2]



    def parse_hydrogen_layer(self):
        for atom in self.atoms:
            self.atom_to_hydrogens[atom] = []

        components = self.get_hydrogen_components()
        print(components)
        for component in components:
            atoms, hydrogen_nr, component_type = self.parse_hydrogen_component(component)
            for atom in atoms:
                if component_type == 'single':
                    self.atom_to_hydrogens[atom].append((component_type, hydrogen_nr))
                else:
                    self.atom_to_hydrogens[atom].append((component_type, hydrogen_nr, atoms))

        pprint(self.atom_to_hydrogens)
        for atom, hydrogens in self.atom_to_hydrogens.items():
            print(atom)
            print(self.graph[atom])
            print(hydrogens)
            print('\n')













