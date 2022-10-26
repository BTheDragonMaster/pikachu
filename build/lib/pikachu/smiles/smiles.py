from typing import *
from pikachu.chem.structure import Structure
from pikachu.errors import StructureError
from pikachu.chem.atom import Atom
# from pikachu.chem.aromatic_system import AromaticSystem


def read_smiles(smiles_string):
    if smiles_string:
        try:
            smiles = Smiles(smiles_string)
            structure = smiles.smiles_to_structure()
            return structure
        except StructureError as e:
            print(f'Error parsing "{smiles_string}": {e.message}')
            return


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

    for i, character in enumerate(informative):
        last = informative[i-1]
        if len(informative) > 2 and i > 1:
            second_last = informative[i-2]
        else:
            second_last = ''
        if skip:
            skip = False
            continue
        if character.isupper():
            if character == 'H':
                if not (len(informative) >= 2 and informative[1] in {'o', 'e', 'f', 's'}):
                    hydrogen = i
                else:
                    element.append(i)
                    element.append(i + 1)
                    skip = True

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
        # Add indices of R groups to the element symbol
        elif character.isdigit() and informative[i-1] in ['R', 'X', 'Z']:
            element.append(i)
        elif (character + last).isdigit() and second_last in ['R', 'X', 'Z']:
            element.append(i)
        elif character.isdigit():
            numbers.append(i)
        elif character == '+' or character == '-':
            charges.append(i)
        elif character == '@':
            chirals.append(i)
        elif character == '*':
            element.append(i)

    element = ''.join([informative[x] for x in element])

    # Parsing the charge

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

    # Parsing an explicit hydrogen

    if hydrogen is not None:
        try:
            if (hydrogen + 1) in numbers:
                hydrogens = int(informative[hydrogen + 1])
            else:
                hydrogens = 1
        except IndexError:
            hydrogens = 1
    else:
        hydrogens = 0

    # Parsing chirality

    if len(chirals) == 1:
        chiral = 'counterclockwise'
    elif len(chirals) == 2:
        chiral = 'clockwise'
    else:
        chiral = None

    # dealing with lone explicit hydrogens

    if not element and hydrogen == 0:
        element = 'H'
        hydrogens = 0

    return element, chiral, charge, hydrogens


def make_character_dict() -> Dict[str, str]:
    """Create dict of {character: label, ->} to label smiles characters
    """
    character_dict = {}
    atoms = ["C", "O", "N", "S", "B", "P", "F", "I", "c", "n", "o", r'*',
             'Cl', 'Br', 'p', 'b', 'p', 's']
    cyclic = list(range(1, 100))

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
    character_dict['.'] = 'split'
    character_dict['-'] = 'single_bond'
    character_dict[':'] = 'aromatic_bond'

    return character_dict


class Smiles:
    character_dict = make_character_dict()
    two_atom_dict = {'B': {'r'}, 'C': {'l'}}

    def __init__(self, string: str) -> None:
        self.smiles = string
        self.components = []
        self.get_components()

    def get_components(self):

        skip = False
        double_digits = False
        square_brackets = False
        component = ''

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

                assert character in {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'}
                component += character
                if len(component) == 2:
                    self.components.append(component)
                    double_digits = False
                    component = ''

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

        bond_labels = {'single_bond',
                       'double_bond',
                       'triple_bond',
                       'quadruple_bond',
                       'chiral_double_bond',
                       'aromatic_bond'}

        structure = Structure()

        # Keeps track of the layer of brackets the atom is in

        branch_level = 0

        # Keeps track of which cyclic atoms have been encountered
        cyclic_dict = {}

        # Keeps track of the last atom encountered per branch level
        last_atoms_dict = {0: None}

        # Keeps track of the nature of the bond
        bond_type = 'single'
        bond_chiral_symbol = None

        # Keeps track of chiral centres
        chiral_dict = {}
        cycle_to_chiral_symbol = {}

        explicit = False
        pyrrole = False
        furan = False
        thiophene = False

        atom_nr = -1
        bond_nr = -1

        for i, component in enumerate(self.components):
            if component[0] == '[':
                label = 'atom'
                explicit = True
            else:
                label = self.character_dict[component]

            if label in bond_labels:
                try:
                    next_component = self.components[i + 1]
                    if next_component[0] == '[':
                        next_label = 'atom'
                    else:
                        next_label = self.character_dict[next_component]

                    if next_label != 'atom' and next_label != 'cyclic':
                        raise StructureError('bond')
                except IndexError:
                    raise StructureError('bond')

            if label == 'split':
                # Starts disconnected structure; set everything back to default
                branch_level = 0
                last_atoms_dict = {0: None}

            elif label == "atom":

                if not explicit:
                    element = component
                    chiral = None
                    charge = 0
                    hydrogens = 0
                else:
                    element, chiral, charge, hydrogens = parse_explicit(component)
                    if element == 'n' and hydrogens == 1:
                        pyrrole = True
                    explicit = False

                if element.islower():
                    aromatic = True
                    element = element.upper()
                    if element == 'O':
                        furan = True
                    elif element == 'S':
                        thiophene = True
                else:
                    aromatic = False

                # Determine atom

                atom_nr += 1

                # Turn atom into an atom object
                atom_2 = Atom(element, atom_nr, chiral, charge, aromatic)
                if pyrrole:
                    atom_2.pyrrole = True
                    pyrrole = False

                if furan:
                    atom_2.furan = True
                    furan = False

                if thiophene:
                    atom_2.thiophene = True
                    thiophene = False

                hydrogen = None

                for j in range(hydrogens):
                    atom_nr += 1
                    bond_nr += 1
                    hydrogen = Atom('H', atom_nr, None, 0, False)
                    structure.add_bond(atom_2, hydrogen, 'single', bond_nr)

                atom_1 = self.get_last_atom(branch_level, last_atoms_dict)

                previous_atom_branch_level = branch_level

                # Go back through the branches to identify the atom that the
                # current atom is attached to; call this atom_1

                while previous_atom_branch_level > 0 and not atom_1:
                    previous_atom_branch_level -= 1
                    atom_1 = last_atoms_dict[previous_atom_branch_level]

                if atom_1:
                    bond_nr += 1

                    if atom_1.aromatic and atom_2.aromatic:
                        if not bond_type == 'explicit_single':
                            bond_type = 'aromatic'
                        else:
                            bond_type = 'single'

                    if bond_type == 'single_chiral':
                        bond_type = 'single'

                    structure.add_bond(atom_1, atom_2, bond_type, bond_nr, bond_chiral_symbol)

                    if atom_1.chiral:
                        # Add current atom to list of residues
                        chiral_dict[atom_1].append(atom_2)

                    if atom_2.chiral:
                        chiral_dict[atom_2] = [atom_1]

                        if hydrogens == 1:
                            chiral_dict[atom_2].append(hydrogen)

                    bond_type = 'single'
                    bond_chiral_symbol = None

                else:

                    # This happens only if atom_2 is completely disconnected
                    # from any atoms already in the graph
                    if hydrogens == 0:
                        structure.add_disconnected_atom(atom_2)

                    if atom_2.chiral:
                        chiral_dict[atom_2] = []
                        if hydrogens == 1:
                            chiral_dict[atom_2].append(hydrogen)

                    # if atom_2.aromatic:
                    #     aromatic_system = AromaticSystem(aromatic_system_id)
                    #     aromatic_system.add_atom(atom_2)
                    #     aromatic_system_id += 1

                # Set atom_2 as the last atom in the current branch level
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

            # If there are brackets: go up or down a branch level

            elif label == "branch_start":
                branch_level += 1
                last_atoms_dict[branch_level] = None

            elif label == "branch_end":
                last_atoms_dict[branch_level] = None
                branch_level -= 1

            elif label == "cyclic":

                cycle_nr = int(component)
                atom = self.get_last_atom(branch_level, last_atoms_dict)

                # If we encounter a number that hasn't been encountered before, start a new cycle

                if self.is_new_cycle(cyclic_dict, cycle_nr):
                    self.start_cycle(cycle_nr, atom, cyclic_dict, bond_type)

                    if atom in chiral_dict:
                        self.add_cycle_placeholder(chiral_dict, atom, cycle_nr)

                    if bond_type == 'single_chiral':
                        cycle_to_chiral_symbol[cycle_nr] = bond_chiral_symbol

                # Otherwise look up the atom that the cycle closes on

                else:
                    bond_nr += 1

                    atom_1, atom_2, old_bond_type = self.end_cycle(cycle_nr, atom, cyclic_dict)

                    if atom_1.aromatic and atom_2.aromatic:
                        if not bond_type == 'explicit_single':
                            bond_type = 'aromatic'
                            # atom_1.aromatic_system.add_atom(atom_2)
                        else:
                            bond_type = 'single'
                            # aromatic_system = AromaticSystem(aromatic_system_id)
                            # aromatic_system.add_atom(atom_2)
                            # aromatic_system_id += 1
                    #
                    # elif atom_2.aromatic:
                    #     aromatic_system = AromaticSystem(aromatic_system_id)
                    #     aromatic_system.add_atom(atom_2)
                    #     aromatic_system_id += 1

                    if bond_type == 'single_chiral':
                        bond_type = 'single'

                        # We have to flip the symbol here, as the previous atom occurred before the double bond

                        if bond_chiral_symbol == '/':
                            bond_chiral_symbol = '\\'
                        elif bond_chiral_symbol == '\\':
                            bond_chiral_symbol = '/'

                        structure.add_bond(atom_1, atom_2, bond_type, bond_nr, bond_chiral_symbol)
                        bond_chiral_symbol = None

                    else:
                        if old_bond_type != 'single':
                            if old_bond_type == 'single_chiral':

                                bond_chiral_symbol = cycle_to_chiral_symbol[cycle_nr]

                                # if bond_chiral_symbol == '/':
                                #     bond_chiral_symbol = '\\'
                                # elif bond_chiral_symbol == '\\':
                                #     bond_chiral_symbol = '/'

                                structure.add_bond(atom_1, atom_2, 'single', bond_nr, bond_chiral_symbol)
                                bond_chiral_symbol = None
                            else:
                                structure.add_bond(atom_1, atom_2, old_bond_type, bond_nr)
                        else:
                            structure.add_bond(atom_1, atom_2, bond_type, bond_nr)

                    if atom_1 in chiral_dict:
                        self.replace_cycle_placeholder(chiral_dict, atom_1, atom_2, cycle_nr)

                    if atom_2 in chiral_dict:
                        chiral_dict[atom_2].append(atom_1)

                bond_type = 'single'

            elif label == 'chiral_double_bond':
                bond_type = 'single_chiral'
                bond_chiral_symbol = component

        structure.refine_structure()
        structure.set_double_bond_chirality()

        for atom in chiral_dict:
            order = chiral_dict[atom]

            # Save chirality in the atom object
            current_chirality = atom.chiral
            new_chirality = self.determine_chirality(order, current_chirality, atom)

            atom.chiral = new_chirality

        return structure

    # =========================================================================
    # Auxillary functions to smiles_to_structure
    # =========================================================================

    @staticmethod
    def add_chiral_atom(chiral_dict: Dict['Atom', Dict[str, Any]], last_atom: 'Atom',
                        current_atom: 'Atom') -> None:
        """Place current_atom in one of the four bond slots of last_atom

        Input:
        chiral_dict: dict of {atom: {'direction': direction, 'order':
            [atom_1, atom_2, atom_3, atom_4]}}, with atom Atom Object,
            direction int, and order a list of Atom Object or int or None
        last_atom: Atom Object, chiral atom
        current_atom: Atom Object, neighbour of chiral atom

        """

        chiral_dict[last_atom].append(current_atom)

    @staticmethod
    def add_cycle_placeholder(chiral_dict, atom, cycle_nr):
        chiral_dict[atom].append(cycle_nr)

    @staticmethod
    def replace_cycle_placeholder(chiral_dict, chiral_atom, current_atom, cycle_nr):
        for i, atom in enumerate(chiral_dict[chiral_atom]):
            if type(atom) == int:
                if atom == cycle_nr:
                    chiral_dict[chiral_atom][i] = current_atom

    @staticmethod
    def get_chiral_permutations(order):
        permutations = [tuple(order),
                        (order[0], order[3], order[1], order[2]),
                        (order[0], order[2], order[3], order[1]),
                        (order[1], order[0], order[3], order[2]),
                        (order[1], order[2], order[0], order[3]),
                        (order[1], order[3], order[2], order[0]),
                        (order[2], order[0], order[1], order[3]),
                        (order[2], order[3], order[0], order[1]),
                        (order[2], order[1], order[3], order[0]),
                        (order[3], order[0], order[2], order[1]),
                        (order[3], order[1], order[0], order[2]),
                        (order[3], order[2], order[1], order[0])]

        return permutations

    def determine_chirality(self, order, chirality, atom):
        if len(order) != 4:
            lone_pairs = atom.lone_pairs

            try:
                assert len(order) + len(lone_pairs) == 4
            except AssertionError:
                raise StructureError('chiral centre')

            original_order = order + lone_pairs

        else:
            original_order = order[:]

        chiral_permutations = self.get_chiral_permutations(original_order)
        original_order.sort(key=lambda x: x.nr)
        new_order = tuple(original_order)
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

    @staticmethod
    def is_new_cycle(cyclic_dict: Dict[int, 'Atom'], cycle_nr: int) -> bool:
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

    @staticmethod
    def start_cycle(cycle_nr: int, atom: 'Atom',
                    cyclic_dict: Dict[int, Tuple['Atom', str]], bond_type: str) -> None:
        """Add a new atom and corresponding cycle number to cyclic dict

        Input:
        cycle_nr: int, nr of the current cycle
        atom: Atom Object
        cyclic_dict: dict of {cycle_nr: atom, ->}, with cycle_nr int and atom
            an Atom Object

        """
        cyclic_dict[cycle_nr] = (atom, bond_type)

    @staticmethod
    def end_cycle(cycle_nr: int, atom: 'Atom',
                  cyclic_dict: Dict[int, 'Atom']) -> Tuple[Union['Atom', None], 'Atom', str]:
        """Return pair of atoms that close a cycle

        Input:
        cycle_nr: int, nr of the current cycle
        atom: Atom Object
        cyclic_dict: dict of {cycle_nr: atom, ->}, with cycle_nr int and atom
            an Atom Object

        Output:
        atom_pair: tuple of two atoms, with each atom an Atom Object

        """
        atom_old, bond_type = cyclic_dict[cycle_nr]

        del cyclic_dict[cycle_nr]

        return atom_old, atom, bond_type

    @staticmethod
    def track_last_atoms_per_branch(new_atom: 'Atom', current_level: int,
                                    last_atoms_dict: Dict[int, Union['Atom', None]]) -> None:
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

    @staticmethod
    def get_last_atom(current_level: int, last_atoms_dict: Dict[int, Union['Atom', None]]) -> 'Atom':
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
