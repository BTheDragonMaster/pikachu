import copy
from pprint import pprint
import sys

from pikachu.chem.bond_properties import BOND_PROPERTIES
from pikachu.chem.atom_properties import ATOM_PROPERTIES
from pikachu.errors import StructureError, KekulisationError
from pikachu.chem.atom import Atom
from pikachu.chem.bond import Bond
from pikachu.chem.kekulisation import Match
from pikachu.chem.substructure_matching import check_same_chirality, compare_all_matches, SubstructureMatch, \
    find_substructures
from pikachu.chem.rings.ring_identification import is_aromatic
import pikachu.chem.rings.find_cycles as find_cycles
from pikachu.chem.aromatic_system import AromaticSystem

sys.setrecursionlimit(1000000)


class Structure:
    """
    Class to store a molecular structure

    Attributes
    ----------
    graph: dict of {atom: [atom, ->], ->}
        with each atom (Atom object) pointing to a list of atoms
        covalently bonded to that atom
    bonds: dict of {bond nr: bond, ->}
        with bond nr (int) an arbitrarily chosen and unique number for a bond,
        and bond (Bond object) the bond associated with that number
    bond_lookup: dict of {atom: {atom: bond, ->}, ->}
        with each pair of atoms (Atom object) pointing to the bond (Bond)
        between those atoms
    """

    def __init__(self, graph=None, bonds=None, bond_lookup=None):

        if graph:
            self.graph = graph
            self.set_atom_neighbours()
            self.set_atoms()
        else:
            self.graph = {}
            self.atoms = {}
        if bonds:
            self.bonds = bonds
            self.make_bond_lookup()
        else:
            self.bonds = {}

        if bond_lookup:
            self.bond_lookup = bond_lookup
        else:
            self.make_bond_lookup()

        self.annotations = set()

        self.aromatic_cycles = []
        self.aromatic_systems = []

    def get_atoms(self):
        atoms = []
        for atom in self.graph:
            atoms.append(atom)

        return atoms

    def bond_exists(self, atom_1, atom_2):
        if atom_1 in self.bond_lookup:
            if atom_2 in self.bond_lookup[atom_1]:
                return True
            return False
        return False

    def get_subtree(self, atom, parent_atom):
        pass

    def set_atoms(self):
        self.atoms = {}
        for atom in self.graph:
            self.atoms[atom.nr] = atom
            
    def deepcopy(self):

        new_graph = {}
        new_bonds = {}
        new_atoms = {}

        for atom_nr, atom in self.atoms.items():
            new_atoms[atom_nr] = atom.copy()
        
        for atom_1, atoms in self.graph.items():
            new_atom_1 = new_atoms[atom_1.nr]
            new_graph[new_atom_1] = []
            for atom_2 in atoms:
                new_atom_2 = new_atoms[atom_2.nr]
                new_graph[new_atom_1].append(new_atom_2)
                
        for bond_nr, bond in self.bonds.items():
            new_atom_1 = new_atoms[bond.atom_1.nr]
            new_atom_2 = new_atoms[bond.atom_2.nr]
            new_bond = Bond(new_atom_1, new_atom_2, bond.type, bond.nr)
            new_bond.chiral = bond.chiral
            new_bond.chiral_symbol = bond.chiral_symbol

            for atom_1, atoms_and_chirality in bond.chiral_dict.items():
                new_1 = None
                if type(atom_1) == Atom:
                    new_1 = new_atoms[atom_1.nr]
                else:
                    pass

                assert new_1

                new_bond.chiral_dict[new_1] = {}

                for atom_2, chirality in atoms_and_chirality.items():
                    if type(atom_2) == Atom:
                        new_2 = new_atoms[atom_2.nr]
                        new_bond.chiral_dict[new_1][new_2] = chirality
                    else:
                        pass

            new_bonds[bond_nr] = new_bond
            if new_bond not in new_atom_1.bonds:
                new_atom_1.bonds.append(new_bond)
            if new_bond not in new_atom_2.bonds:
                new_atom_2.bonds.append(new_bond)
                
        new_structure = Structure(new_graph, new_bonds)

        new_structure.add_shells_non_hydrogens()
        new_structure.add_shells()

        new_structure.form_pi_bonds()
        new_structure.hybridise_atoms()
        new_structure.promote_pi_bonds()

        new_structure.find_cycles()
        new_structure.aromatic_cycles = new_structure.find_aromatic_cycles()
        new_structure.aromatic_systems = new_structure.find_aromatic_systems()

        new_structure.form_sigma_bonds()
        new_structure.drop_electrons()
        new_structure.make_lone_pairs()

        new_structure.set_atom_neighbours()
        new_structure.make_bond_lookup()
        new_structure.set_connectivities()

        for annotation, default in self.annotations:
            new_structure.annotations.add((annotation, default))

        return new_structure

    def get_priority_groups(self, priority_groups):
        ordered_priorities = sorted(priorities, reverse=True)
        priority_groups = {0: [],
                           1: [],
                           2: [],
                           3: []}

        previous_priority = None
        previous_priority_idx = 0

        for i, priority in enumerate(ordered_priorities):
            if priority != previous_priority:
                previous_priority = priority
                previous_priority_idx = i
                priority_groups[i] = [priority]
            else:
                priority_groups[previous_priority_idx].append(priority)


    def get_absolute_chirality(self, chiral_center):
        assert chiral_center.chiral

        chiral_atom = self.get_atom(chiral_center)
        priorities = []
        next_atoms = []
        masked_atoms = {chiral_atom}
        for atom in chiral_atom.neighbours:
            priorities.append((ATOM_PROPERTIES.element_to_atomic_nr[atom.type], 1))
            next_atoms.append(atom)



        unresolved_priority_groups = []



        while len(set(priorities)) != len(priorities):


            for priority in priorities:
                pass


    def copy(self):
        new_graph = {}
        new_bonds = {}

        for atom_1, atoms in self.graph.items():
            new_graph[atom_1] = []
            for atom_2 in atoms:
                new_graph[atom_1].append(self.atoms[atom_2.nr])

        for bond_nr, bond in self.bonds.items():
            new_bonds[bond_nr] = bond

        new_structure = Structure(new_graph, new_bonds)
        new_structure.cycles = self.cycles

        return new_structure

    def refresh_structure(self, find_cycles=False):
        new_graph = {}
        self.set_atoms()

        for atom_1, atoms in self.graph.items():
            new_graph[atom_1] = []
            for atom_2 in atoms:
                new_graph[atom_1].append(self.atoms[atom_2.nr])

        for bond_nr, bond in self.bonds.items():
            bond.atom_1 = self.atoms[bond.atom_1.nr]
            bond.atom_2 = self.atoms[bond.atom_2.nr]
            if bond not in bond.atom_1.bonds:
                bond.atom_1.bonds.append(bond)
            if bond not in bond.atom_2.bonds:
                bond.atom_2.bonds.append(bond)

        self.graph = new_graph

        self.make_bond_lookup()
        self.set_atom_neighbours()
        self.set_connectivities()
        for bond in self.bonds.values():
            bond.set_bond_summary()

        if find_cycles:
            self.find_cycles()

    def get_next_in_ring(self, ring, current_atom, previous_atom):
        neighbours = self.graph[current_atom]

        for neighbour in neighbours:
            for member in ring.members:
                if neighbour == member:
                    if previous_atom != neighbour:
                        return neighbour

        return None

    def colour_substructure_single(self, substructure, colour="hot pink", check_chiral_centres=True,
                                   check_bond_chirality=True):
        matches = self.find_substructures(substructure,
                                          check_chiral_centres=check_chiral_centres,
                                          check_chiral_double_bonds=check_bond_chirality)

        if matches:
            match = matches[0]
            for parent_atom in match.atoms.values():
                for atom in self.graph:
                    if atom == parent_atom:
                        atom.draw.colour = colour

    def colour_substructure_all(self, substructure, colour="hot pink", check_chiral_centres=True,
                                check_bond_chirality=True):
        matches = self.find_substructures(substructure,
                                          check_chiral_centres=check_chiral_centres,
                                          check_chiral_double_bonds=check_bond_chirality)

        for match in matches:
            for parent_atom in match.atoms.values():
                for atom in self.graph:
                    if atom == parent_atom:
                        atom.draw.colour = colour

    def to_dash_molecule2d_input(self):
        nodes = []
        links = []

        kekulised_structure = self.kekulise()

        for atom in kekulised_structure.graph:
            atom_dict = {'id': atom.nr,
                         'atom': atom.type}
            nodes.append(atom_dict)

        for bond_nr, bond in kekulised_structure.bonds.items():
            assert bond.type != 'aromatic'
            bond_dict = {'id': bond_nr,
                         'source': bond.atom_1.nr,
                         'target': bond.atom_2.nr,
                         'bond': BOND_PROPERTIES.bond_type_to_weight[bond.type]}
            links.append(bond_dict)

        dash_molecule2d_input = {'nodes': nodes, 'links': links}
        return dash_molecule2d_input

    def get_drawn_atoms(self):
        atoms = []
        for atom in self.graph:
            if atom.draw.is_drawn:
                atoms.append(atom)
        return atoms

    def get_drawn_bonds(self):
        bonds = []

        for bond_nr, bond in self.bonds.items():
            if bond.atom_1.draw.is_drawn and bond.atom_2.draw.is_drawn:
                bonds.append(bond)

        return bonds

    def set_double_bond_chirality(self):

        for bond_nr, bond in self.bonds.items():

            # iterate over all double bonds

            if bond.type == 'double' or bond.type == 'triple':

                # double bonds neighboured by three bonds on each atom, e.g. a C=C bond

                if len(bond.atom_1.bonds) + len(bond.atom_1.lone_pairs) == 3 and \
                        len(bond.atom_2.bonds) + len(bond.atom_2.lone_pairs) == 3 and \
                        len(bond.atom_1.lone_pairs) < 2 and len(bond.atom_2.lone_pairs) < 2:

                    # define atoms adjacent to the atoms involved in the double bond
                    # also keep track of the chiral symbol that defines these bonds

                    atom_1_1 = None
                    atom_1_2 = None
                    chiral_1_1 = None
                    chiral_1_2 = None

                    atom_2_1 = None
                    atom_2_2 = None
                    chiral_2_1 = None
                    chiral_2_2 = None

                    if bond.atom_1.lone_pairs:
                        atom_1_1 = bond.atom_1.lone_pairs[0]
                    if bond.atom_2.lone_pairs:
                        atom_2_1 = bond.atom_2.lone_pairs[0]

                    # Check bonds adjacent to the first atom

                    for bond_1 in bond.atom_1.bonds:

                        if bond_1.type == 'single':

                            # Looks at the bonds between the atom adjacent to the stereobond and its neighbours
                            if bond.atom_1 == bond_1.atom_1:
                                if bond_1.chiral_symbol == '/':
                                    direction = 'up'
                                elif bond_1.chiral_symbol == '\\':
                                    direction = 'down'
                                else:
                                    direction = None

                                # First time it runs through this it will define atom_1_1

                                if not atom_1_1:
                                    atom_1_1 = bond_1.atom_2
                                    chiral_1_1 = direction

                                # Second time it runs through this, it will define atom_1_2
                                else:
                                    atom_1_2 = bond_1.atom_2
                                    chiral_1_2 = direction

                            elif bond.atom_1 == bond_1.atom_2:

                                if bond_1.chiral_symbol == '/':
                                    direction = 'down'
                                elif bond_1.chiral_symbol == '\\':
                                    direction = 'up'
                                else:
                                    direction = None

                                if not atom_1_1:
                                    atom_1_1 = bond_1.atom_1
                                    chiral_1_1 = direction
                                else:
                                    atom_1_2 = bond_1.atom_1
                                    chiral_1_2 = direction

                    for bond_2 in bond.atom_2.bonds:

                        if bond_2.type == 'single':
                            if bond.atom_2 == bond_2.atom_1:
                                if bond_2.chiral_symbol == '/':
                                    direction = 'up'
                                elif bond_2.chiral_symbol == '\\':
                                    direction = 'down'
                                else:
                                    direction = None

                                if not atom_2_1:
                                    atom_2_1 = bond_2.atom_2
                                    chiral_2_1 = direction
                                else:
                                    atom_2_2 = bond_2.atom_2
                                    chiral_2_2 = direction

                            elif bond.atom_2 == bond_2.atom_2:
                                if bond_2.chiral_symbol == '/':
                                    direction = 'down'
                                elif bond_2.chiral_symbol == '\\':
                                    direction = 'up'
                                else:
                                    direction = None

                                if not atom_2_1:
                                    atom_2_1 = bond_2.atom_1
                                    chiral_2_1 = direction
                                else:
                                    atom_2_2 = bond_2.atom_1
                                    chiral_2_2 = direction

                    chiral_1 = False
                    chiral_2 = False

                    if chiral_1_1 or chiral_1_2:
                        chiral_1 = True

                    if chiral_2_1 or chiral_2_2:
                        chiral_2 = True

                    if chiral_1 and chiral_2:
                        if chiral_1_1 == chiral_1_2:
                            raise StructureError('chiral double bond')
                        if chiral_2_2 == chiral_2_1:
                            raise StructureError('chiral double bond')

                        if chiral_1_1:
                            first_atom = atom_1_1
                            first_other_atom = atom_1_2
                            first_chiral_symbol = chiral_1_1

                            if not chiral_1_2 and type(atom_1_2) == Atom:
                                # Make sure where chiral symbols are not defined, they are added

                                if (atom_1_1.nr > bond.atom_1.nr and atom_1_2.nr > bond.atom_1.nr) or \
                                        (atom_1_1.nr < bond.atom_1.nr and atom_1_2.nr < bond.atom_1.nr):
                                    if self.bond_lookup[bond.atom_1][atom_1_1].chiral_symbol == '/':
                                        self.bond_lookup[bond.atom_1][atom_1_2].chiral_symbol = '\\'
                                    else:
                                        self.bond_lookup[bond.atom_1][atom_1_2].chiral_symbol = '/'

                                else:
                                    if self.bond_lookup[bond.atom_1][atom_1_1].chiral_symbol == '/':
                                        self.bond_lookup[bond.atom_1][atom_1_2].chiral_symbol = '/'
                                    else:
                                        self.bond_lookup[bond.atom_1][atom_1_2].chiral_symbol = '\\'

                        else:
                            first_atom = atom_1_2
                            first_other_atom = atom_1_1
                            first_chiral_symbol = chiral_1_2

                            # Make sure where chiral symbols are not defined, they are added

                            if type(atom_1_1) == Atom:

                                if (atom_1_1.nr > bond.atom_1.nr and atom_1_2.nr > bond.atom_1.nr) or \
                                   (atom_1_1.nr < bond.atom_1.nr and atom_1_2.nr < bond.atom_1.nr):
                                    if self.bond_lookup[bond.atom_1][atom_1_2].chiral_symbol == '/':
                                        self.bond_lookup[bond.atom_1][atom_1_1].chiral_symbol = '\\'
                                    else:
                                        self.bond_lookup[bond.atom_1][atom_1_1].chiral_symbol = '/'

                                else:
                                    if self.bond_lookup[bond.atom_1][atom_1_2].chiral_symbol == '/':
                                        self.bond_lookup[bond.atom_1][atom_1_1].chiral_symbol = '/'
                                    else:
                                        self.bond_lookup[bond.atom_1][atom_1_1].chiral_symbol = '\\'

                        if chiral_2_1:
                            second_atom = atom_2_1
                            second_other_atom = atom_2_2
                            second_chiral_symbol = chiral_2_1

                            if not chiral_2_2 and type(atom_2_2) == Atom:

                                # Make sure where chiral symbols are not defined, they are added

                                if (atom_2_1.nr > bond.atom_2.nr and atom_2_2.nr > bond.atom_2.nr) or \
                                        (atom_2_1.nr < bond.atom_2.nr and atom_2_2.nr < bond.atom_2.nr):
                                    if self.bond_lookup[bond.atom_2][atom_2_1].chiral_symbol == '/':
                                        self.bond_lookup[bond.atom_2][atom_2_2].chiral_symbol = '\\'
                                    else:
                                        self.bond_lookup[bond.atom_2][atom_2_2].chiral_symbol = '/'

                                else:
                                    if self.bond_lookup[bond.atom_2][atom_2_1].chiral_symbol == '/':
                                        self.bond_lookup[bond.atom_2][atom_2_2].chiral_symbol = '/'
                                    else:
                                        self.bond_lookup[bond.atom_2][atom_2_2].chiral_symbol = '\\'

                        else:
                            second_atom = atom_2_2
                            second_other_atom = atom_2_1
                            second_chiral_symbol = chiral_2_2

                            # Make sure where chiral symbols are not defined, they are added

                            if type(atom_2_1) == Atom:
                                if (atom_2_1.nr > bond.atom_2.nr and atom_2_2.nr > bond.atom_2.nr) or \
                                        (atom_2_1.nr < bond.atom_2.nr and atom_2_2.nr < bond.atom_2.nr):
                                    if self.bond_lookup[bond.atom_2][atom_2_2].chiral_symbol == '/':
                                        self.bond_lookup[bond.atom_2][atom_2_1].chiral_symbol = '\\'
                                    else:
                                        self.bond_lookup[bond.atom_2][atom_2_1].chiral_symbol = '/'

                                else:
                                    if self.bond_lookup[bond.atom_2][atom_2_2].chiral_symbol == '/':
                                        self.bond_lookup[bond.atom_2][atom_2_1].chiral_symbol = '/'
                                    else:
                                        self.bond_lookup[bond.atom_2][atom_2_1].chiral_symbol = '\\'

                        if type(first_atom) == Atom:
                            bond.chiral_dict[first_atom] = {}
                        if type(first_other_atom) == Atom:
                            bond.chiral_dict[first_other_atom] = {}
                        if type(second_atom) == Atom:
                            bond.chiral_dict[second_atom] = {}
                        if type(second_other_atom) == Atom:
                            bond.chiral_dict[second_other_atom] = {}

                        if first_chiral_symbol == second_chiral_symbol:
                            if type(first_atom) == Atom and type(second_atom) == Atom:
                                bond.chiral_dict[first_atom][second_atom] = 'cis'
                                bond.chiral_dict[second_atom][first_atom] = 'cis'

                            if type(first_other_atom) == Atom and type(second_other_atom) == Atom:
                                bond.chiral_dict[first_other_atom][second_other_atom] = 'cis'
                                bond.chiral_dict[second_other_atom][first_other_atom] = 'cis'

                            if type(first_atom) == Atom and type(second_other_atom) == Atom:
                                bond.chiral_dict[first_atom][second_other_atom] = 'trans'
                                bond.chiral_dict[second_other_atom][first_atom] = 'trans'

                            if type(first_other_atom) == Atom and type(second_atom) == Atom:
                                bond.chiral_dict[first_other_atom][second_atom] = 'trans'
                                bond.chiral_dict[second_atom][first_other_atom] = 'trans'

                        else:
                            if type(first_atom) == Atom and type(second_atom) == Atom:
                                bond.chiral_dict[first_atom][second_atom] = 'trans'
                                bond.chiral_dict[second_atom][first_atom] = 'trans'

                            if type(first_other_atom) == Atom and type(second_other_atom) == Atom:
                                bond.chiral_dict[first_other_atom][second_other_atom] = 'trans'
                                bond.chiral_dict[second_other_atom][first_other_atom] = 'trans'

                            if type(first_atom) == Atom and type(second_other_atom) == Atom:
                                bond.chiral_dict[first_atom][second_other_atom] = 'cis'
                                bond.chiral_dict[second_other_atom][first_atom] = 'cis'

                            if type(first_other_atom) == Atom and type(second_atom) == Atom:
                                bond.chiral_dict[first_other_atom][second_atom] = 'cis'
                                bond.chiral_dict[second_atom][first_other_atom] = 'cis'

                        bond.chiral = True

    def find_next_atom_nr(self):
        """
        Return the next available integer to label an atom

        Output
        ------
        next_atom_nr: int
        """
        atom_nrs = []

        for atom in self.graph:
            atom_nrs.append(atom.nr)

        atom_nrs.sort()

        next_atom_nr = None

        for i, atom_nr in enumerate(atom_nrs):
            if i != atom_nr:
                next_atom_nr = i
                break

        if next_atom_nr is None:
            next_atom_nr = atom_nrs[-1] + 1

        return next_atom_nr

    def find_highest_atom_nr(self):
        highest_atom_nr = -1
        for atom in self.graph:
            if atom.nr > highest_atom_nr:
                highest_atom_nr = atom.nr

        return highest_atom_nr

    def find_highest_bond_nr(self):
        highest_bond_nr = -1
        for bond_nr in self.bonds:
            if bond_nr > highest_bond_nr:
                highest_bond_nr = bond_nr

        return highest_bond_nr

    def get_steric_atoms(self):
        """
        Return list of atoms that are chiral centres

        Output
        ------
        steric atoms: list of [atom, ->]
            with each atom (Atom object) a chiral centre in the structure
        """
        steric_atoms = []
        for atom in self.graph:
            if atom.chiral:
                steric_atoms.append(atom)

        return steric_atoms

    # TODO:
    # call find_next_bond_nr in make_bond

    def find_next_bond_nr(self):
        """
        Return the next available integer to label a bond

        Output
        ------
        next_bond_nr: int
        """
        bond_nrs = list(self.bonds.keys())
        bond_nrs.sort()

        next_bond_nr = None

        for i, bond_nr in enumerate(bond_nrs):
            if i != bond_nr:
                next_bond_nr = i
                break

        if next_bond_nr is None:
            next_bond_nr = bond_nrs[-1] + 1

        return next_bond_nr

    def make_bond_lookup(self):
        """
        Update bond_lookup from current collection of bonds
        """
        self.bond_lookup = {}
        for bond in self.bonds:
            atom_1 = self.bonds[bond].atom_1
            atom_2 = self.bonds[bond].atom_2
            if atom_1 not in self.bond_lookup:
                self.bond_lookup[atom_1] = {}
            if atom_2 not in self.bond_lookup:
                self.bond_lookup[atom_2] = {}
            self.bond_lookup[atom_1][atom_2] = self.bonds[bond]
            self.bond_lookup[atom_2][atom_1] = self.bonds[bond]

    def make_atom_index(self):
        """
        Return dict of {atom nr: atom}

        Output
        ------
        atom_nr_to_atom: dict of {atom nr: atom}
            with atom nr (int) the unique number associated with an atom based
            on the order they occurred in in the input SMILES, and atom
            (Atom object) the atom associated with that number
        """
        atom_nr_to_atom = {}
        for atom in self.graph:
            atom_nr_to_atom[atom.nr] = atom
        return atom_nr_to_atom

    def add_atom(self, atom_type, neighbours, chiral=None, charge=0, aromatic=False):
        """
        Add a new atom to an already existing structure

        Input
        -----
        atom_type: str, an abbreviation of an element of the periodic table
        neighbours: list of Atom objects
        chiral: str or None, str indicating the chirality (clockwise or
            counterclockwise) if the new atom is a chiral centre, None
            if the new atom is not a chiral centre
        charge: int, charge of the new atom
        """

        next_atom_nr = self.find_next_atom_nr()
        atom = Atom(atom_type, next_atom_nr, chiral, charge, aromatic)
        atom.add_electron_shells()

        for annotation, default in self.annotations:
            atom.annotations.add_annotation(annotation, default)

        for i, neighbour in enumerate(neighbours):
            next_bond_nr = self.find_next_bond_nr()
            self.make_bond(atom, neighbour, next_bond_nr)

        atom.hybridise()

        self.atoms[atom.nr] = atom

        return atom

    def find_cycles(self):
        """
        Find all cycles in a structure and store them
        """
        self.cycles = find_cycles.Cycles(self)
        self.sssr = self.cycles.find_sssr()

    @staticmethod
    def promote_lone_pairs_in_aromatic_cycles(cycles):
        for cycle in cycles:
            for atom in cycle:
                if atom.hybridisation == 'sp3':
                    atom.promote_lone_pair_to_p_orbital()
                    if atom.type == 'N':
                        atom.pyrrole = True
                    elif atom.type == 'S':
                        atom.thiophene = True
                    elif atom.type == 'O':
                        atom.furan = True

    def get_bounding_box(self):
        min_x = 1000000000
        max_x = -1000000000
        min_y = 1000000000
        max_y = -1000000000

        for atom in self.graph:
            if atom.draw.is_drawn:
                if atom.draw.position.x < min_x:
                    min_x = atom.draw.position.x
                if atom.draw.position.x > max_x:
                    max_x = atom.draw.position.x
                if atom.draw.position.y < min_y:
                    min_y = atom.draw.position.y
                if atom.draw.position.y > max_y:
                    max_y = atom.draw.position.y
                    
        return min_x, min_y, max_x, max_y

    def find_aromatic_cycles(self):
        """
        Returns cycles that are aromatic

        Output
        ------
        aromatic_cycles: list of [[Atom, ->], ->], with each list of
            atoms the atoms that comprise an aromatic cycle

        """
        aromatic_cycles = set()
        previous_nr_aromatic_cycles = -1
        current_nr_aromatic_cycles = 0

        while previous_nr_aromatic_cycles != current_nr_aromatic_cycles:
            previous_nr_aromatic_cycles = current_nr_aromatic_cycles
            for cycle in self.sssr:
                if tuple(cycle) not in aromatic_cycles and is_aromatic(cycle):
                    self.make_cycle_aromatic(cycle)
                    aromatic_cycles.add(tuple(cycle))

            current_nr_aromatic_cycles = len(aromatic_cycles)

        self.promote_lone_pairs_in_aromatic_cycles(aromatic_cycles)

        return aromatic_cycles

    def find_aromatic_systems(self):
        """
        Return ring systems that are aromatic

        Output
        ------
        aromatic_ring_systems: list of [[Atom, ->], ->], with each list of
            atoms the atoms that comprise an aromatic ring system
        """
        previous_system_nr = -1
        current_system_nr = 1

        aromatic_systems = [list(aromatic_cycle) for aromatic_cycle in self.aromatic_cycles]

        while current_system_nr != previous_system_nr:
            previous_system_nr = current_system_nr
            indices_to_remove = None
            new_system = None
            for i, system_1 in enumerate(aromatic_systems):
                system_found = False
                for j, system_2 in enumerate(aromatic_systems):
                    if i != j:
                        if len(set(system_1).intersection(set(system_2))) >= 2:
                            indices_to_remove = [i, j]
                            new_system = list(set(system_1[:] + system_2[:]))
                            system_found = True
                            break
                if system_found:
                    break

            if new_system:
                indices_to_remove.sort(reverse=True)
                for index in indices_to_remove:
                    aromatic_systems.pop(index)

                aromatic_systems.append(new_system)
                
            current_system_nr = len(aromatic_systems)

        aromatic_ring_systems = []

        for i, aromatic_system in enumerate(aromatic_systems):
            system = AromaticSystem(i, aromatic_system)
            aromatic_ring_systems.append(system)

        return aromatic_ring_systems

    def find_double_bond_sequences(self):
        double_bond_fragments = []
        for bond in self.bonds.values():
            if bond.type == 'single':
                stereobond_1 = None
                stereobond_2 = None
                for bond_1 in bond.atom_1.bonds:
                    if bond_1.chiral:
                        stereobond_1 = bond_1
                        break

                for bond_2 in bond.atom_2.bonds:
                    if bond_2.chiral:
                        stereobond_2 = bond_2
                        break

                if stereobond_1 and stereobond_2:
                    double_bond_fragments.append([stereobond_1, stereobond_2])

        previous_fragment_nr = -1
        fragment_nr = len(double_bond_fragments)

        while fragment_nr != previous_fragment_nr:
            previous_fragment_nr = fragment_nr

            indices_to_remove = None
            new_fragment = None

            for i, fragment_1 in enumerate(double_bond_fragments):
                found_fragment = False
                for j, fragment_2 in enumerate(double_bond_fragments):
                    if i != j:
                        if fragment_1[-1] == fragment_2[0]:
                            new_fragment = fragment_1[:] + fragment_2[1:]
                        elif fragment_1[-1] == fragment_2[-1]:
                            new_fragment = fragment_1[:] + list(reversed(fragment_2[:-1]))
                        elif fragment_1[0] == fragment_2[0]:
                            new_fragment = list(reversed(fragment_2[1:])) + fragment_1[:]
                        elif fragment_1[0] == fragment_2[-1]:
                            new_fragment = fragment_2[:-1] + fragment_1[:]

                        if new_fragment:
                            indices_to_remove = [i, j]
                            found_fragment = True
                            break
                if found_fragment:
                    break

            if indices_to_remove:
                assert new_fragment
                # sort indices from large to small to make sure you first remove the last index
                indices_to_remove.sort(reverse=True)
                double_bond_fragments.pop(indices_to_remove[0])
                double_bond_fragments.pop(indices_to_remove[1])
                double_bond_fragments.append(new_fragment)
                fragment_nr = len(double_bond_fragments)

        return double_bond_fragments

    @staticmethod
    def make_cycle_aromatic(cycle):
        for atom_1 in cycle:
            atom_1.aromatic = True
            for atom_2 in cycle:
                if atom_1 != atom_2:
                    bond = atom_1.get_bond(atom_2)
                    if bond and bond.type != 'aromatic':
                        bond.make_aromatic()

    def refine_structure(self):
        """
        """

        self.add_shells()
        self.add_hydrogens()
        self.add_shells()

        self.sort_by_nr()

        self.form_pi_bonds()
        self.hybridise_atoms()
        self.promote_pi_bonds()

        self.find_cycles()
        self.aromatic_cycles = self.find_aromatic_cycles()
        self.aromatic_systems = self.find_aromatic_systems()

        self.form_sigma_bonds()
        self.drop_electrons()
        self.make_lone_pairs()

        self.set_atom_neighbours()

        self.make_bond_lookup()

        self.set_connectivities()
        self.set_atoms()

    def make_lone_pairs(self):
        for atom in self.graph:
            atom.make_lone_pairs()

    def promote_pi_bonds(self):
        for atom in self.graph:
            atom.promote_pi_bonds_to_d_orbitals()

    def remove_bond_between_atoms(self, atom_1, atom_2):
        bond = self.bond_lookup[atom_1][atom_2]

        del self.bond_lookup[atom_1][atom_2]
        del self.bond_lookup[atom_2][atom_1]
        del self.bonds[bond.nr]

        if atom_2 in self.graph[atom_1]:
            self.graph[atom_1].remove(atom_2)
        if atom_1 in self.graph[atom_2]:
            self.graph[atom_2].remove(atom_1)

    def break_bond(self, bond):
        """
        Break a bond in the structure

        Input
        -----
        bond: Bond object
        """

        atom_1 = bond.atom_1
        atom_2 = bond.atom_2

        bond.break_bond()

        del self.bond_lookup[atom_1][atom_2]
        del self.bond_lookup[atom_2][atom_1]
        del self.bonds[bond.nr]

        if atom_2 in self.graph[atom_1]:
            self.graph[atom_1].remove(atom_2)
        if atom_1 in self.graph[atom_2]:
            self.graph[atom_2].remove(atom_1)

    def break_bond_by_nr(self, bond):
        """
        Break a bond in the structure based on bond number

        Input
        -----
        bond: Bond object or int, with int the number of a bond in the graph
        """
        if type(bond) == int:
            bond = self.bonds[bond]

        self.break_bond(bond)

    def break_bond_between_atoms(self, atom_1, atom_2):
        """
        Break a bond in the structure based on the numbers of neighbouring
        atoms

        Input
        -----
        atom_1: Atom object or int, with int the number of an atom in the graph
        atom_2: Atom object or int, with int the number of an atom in the graph

        """
        atom_nr_to_atom = self.make_atom_index()

        if atom_1.type == int:
            atom_1 = atom_nr_to_atom[atom_1]
        if atom_2.type == int:
            atom_2 = atom_nr_to_atom[atom_2]

        bond = self.bond_lookup[atom_1][atom_2]
        self.break_bond(bond)

    def remove_atom(self, atom_to_remove):
        """
        Remove an atom from the structure

        Input
        -----
        atom_to_remove: Atom object, atom that is to be removed from the
            structure
        """

        for bond in atom_to_remove.bonds:
            self.break_bond_by_nr(bond.nr)

        for atom in self.graph:
            if atom_to_remove in self.graph[atom]:
                self.graph[atom].remove(atom_to_remove)

        del self.graph[atom_to_remove]

    def set_connectivities(self):
        for atom in self.graph:
            if atom.type != 'H':
                atom.set_connectivity()

    def set_atom_neighbours(self):
        for atom in self.graph:
            atom.set_neighbours(self)

    def get_connectivities(self):
        connectivities = {}
        for atom in self.graph:
            if atom.type != 'H':
                connectivity = atom.connectivity
                if connectivity not in connectivities:
                    connectivities[connectivity] = []
                    connectivities[connectivity].append(atom)

        return connectivities

    def get_chiral_double_bonds(self):
        chiral_bonds = []
        for bond_nr in self.bonds:
            bond = self.bonds[bond_nr]
            if bond.chiral:
                chiral_bonds.append(bond)

        return chiral_bonds

    def check_chiral_double_bonds(self, child, match):
        chirality_matches = True
        chiral_bonds = child.get_chiral_double_bonds()

        for chiral_bond in chiral_bonds:
            neighbour_1, neighbour_2 = chiral_bond.neighbours

            parent_neighbour_1 = match[neighbour_1]
            parent_neighbour_2 = match[neighbour_2]
            parent_bond = self.bond_lookup[parent_neighbour_1][parent_neighbour_2]

            if not parent_bond.chiral:
                chirality_matches = False
                break
            else:
                matching_chirality = chiral_bond.check_same_chirality(parent_bond, match)
                if not matching_chirality:
                    chirality_matches = False
                    break

        return chirality_matches

    @staticmethod
    def check_chiral_centres(child, match):
        chirality_matches = True
        chiral_centres = child.get_steric_atoms()

        for chiral_centre in chiral_centres:
            parent_atom = match[chiral_centre]
            if parent_atom.chiral:

                chirality_matches = check_same_chirality(chiral_centre, parent_atom, match)
                if not chirality_matches:
                    break
            else:
                chirality_matches = False
                break

        return chirality_matches

    def is_substructure_bond_composition(self, substructure):
        bond_summary_to_count_parent = {}
        bond_summary_to_count_child = {}

        for bond_nr, bond in self.bonds:
            bond_summary = bond.bond_summary
            if 'H' not in [atom.type for atom in bond.neighbours]:
                if bond_summary not in bond_summary_to_count_parent:
                    bond_summary_to_count_parent[bond_summary] = 0
                bond_summary_to_count_parent[bond_summary] += 1
        for bond_nr, bond in substructure.bonds:
            if 'H' not in [atom.type for atom in bond.neighbours]:
                bond_summary = bond.bond_summary
                if bond_summary not in bond_summary_to_count_child:
                    bond_summary_to_count_child[bond_summary] = 0
                bond_summary_to_count_child[bond_summary] += 1

        can_be_substructure = True

        for bond_summary, count in bond_summary_to_count_child.items():
            if bond_summary not in bond_summary_to_count_parent:
                can_be_substructure = False
                break
            elif count > bond_summary_to_count_parent[bond_summary]:
                can_be_substructure = False
                break

        return can_be_substructure

    def find_substructures(self, substructure, check_chiral_centres=True, check_chiral_double_bonds=True):
        matches = []
        if self.is_substructure_atom_composition(substructure):
            if self.is_substructure_atom_connectivity(substructure):
                matches = find_substructures(self, substructure)

        if check_chiral_centres:
            final_matches = []

            for match in matches:
                if self.check_chiral_centres(substructure, match.atoms):
                    final_matches.append(match)
            matches = final_matches

        if check_chiral_double_bonds:
            final_matches = []
            for match in matches:
                if self.check_chiral_double_bonds(substructure, match.atoms):
                    final_matches.append(match)
            matches = final_matches

        return matches

    def get_atom(self, atom):
        return self.atoms[atom.nr]

    def get_bond(self, bond):
        return self.bonds[bond.nr]

    def is_substructure_atom_composition(self, child):

        atom_counts_self = self.get_atom_counts()
        atom_counts_child = child.get_atom_counts()

        can_be_substructure = True

        for atom_type in atom_counts_child:
            try:
                atom_nr_self = atom_counts_self[atom_type]
                atom_nr_child = atom_counts_child[atom_type]
                if atom_nr_child > atom_nr_self:
                    can_be_substructure = False
                    break
            except KeyError:
                can_be_substructure = False
                break

        return can_be_substructure

    def get_connectivity_counts(self):
        connectivities = {}
        for atom in self.graph:
            if atom.type != 'H':
                if atom.type not in connectivities:
                    connectivities[atom.type] = {}
                connectivity = atom.connectivity
                if connectivity not in connectivities[atom.type]:
                    connectivities[atom.type][connectivity] = 0
                connectivities[atom.type][connectivity] += 1

        return connectivities

    def get_substructure_connectivity_counts(self, atom_connectivities_child):
        substructure_connectivity_counts = {}
        for atom_type in atom_connectivities_child:
            substructure_connectivity_counts[atom_type] = {}
            for connectivity in atom_connectivities_child[atom_type]:
                substructure_connectivity_counts[atom_type][connectivity] = 0
                for atom in self.graph:
                    if atom.type == atom_type and atom.potential_same_connectivity(connectivity):
                        substructure_connectivity_counts[atom_type][connectivity] += 1

        return substructure_connectivity_counts

    def get_substructure_connectivities(self, atom_connectivities_child):
        substructure_connectivities = {}
        for substructure_connectivity in atom_connectivities_child:
            substructure_connectivities[substructure_connectivity] = []
            for atom in self.graph:
                if atom.type != 'H' and atom.potential_same_connectivity(substructure_connectivity):
                    substructure_connectivities[substructure_connectivity].append(atom)

        return substructure_connectivities

    def is_substructure_atom_connectivity(self, child):

        atom_connectivities_child = child.get_connectivity_counts()
        atom_connectivities_self = self.get_substructure_connectivity_counts(atom_connectivities_child)

        can_be_substructure = True

        for atom_type in atom_connectivities_child:
            for connectivity in atom_connectivities_child[atom_type]:
                connectivity_nr_self = atom_connectivities_self[atom_type][connectivity]
                connectivity_nr_child = atom_connectivities_child[atom_type][connectivity]
                if connectivity_nr_child > connectivity_nr_self:
                    can_be_substructure = False
                    break

        return can_be_substructure

    def make_bond_dict(self):
        bond_dict = {}

        for atom in self.graph:
            if atom.type != 'H':
                bond_dict[atom] = 0
                for neighbour in atom.neighbours:
                    if neighbour.type != 'H':
                        bond_dict[atom] += 1

        return bond_dict

    def drop_electrons(self):
        for atom in self.graph:
            atom.drop_electrons()

    def add_shells_non_hydrogens(self):
        for atom in self.graph:
            if atom.type != 'H':
                atom.add_electron_shells()

    def add_shells(self):
        for atom in self.graph:
            if not atom.shells:
                atom.add_electron_shells()

    def hybridise_atoms(self, atoms=None):

        if not atoms:
            for atom in self.graph:
                atom.hybridise()

        else:
            for atom in atoms:
                atom.hybridise()

    def add_hydrogens(self):
        atom_numbers = [atom.nr for atom in self.graph]

        if len(atom_numbers) > 0:
            max_atom_nr = max(atom_numbers)
        else:
            max_atom_nr = -1

        bond_numbers = list(self.bonds.keys())

        if len(bond_numbers) > 0:
            max_bond_nr = max(bond_numbers)

        else:
            max_bond_nr = -1

        for atom in list(self.graph.keys()):
            hydrogens_to_add = atom.calc_hydrogens()
            for i in range(hydrogens_to_add):
                max_atom_nr += 1
                max_bond_nr += 1
                hydrogen = Atom('H', max_atom_nr, None, 0, False)
                self.add_bond(atom, hydrogen, 'single', max_bond_nr)

    def form_pi_bonds(self):
        for bond_nr in self.bonds:
            bond = self.bonds[bond_nr]
            if bond.type != 'single':
                bond.combine_p_orbitals()

    def form_sigma_bonds(self):
        for bond_nr, bond in self.bonds.items():
            bond.combine_hybrid_orbitals()

    def get_atom_counts(self):
        atom_counts = {}
        for atom in self.graph:
            if atom.type != 'H':
                if atom.type not in atom_counts:
                    atom_counts[atom.type] = 0
                atom_counts[atom.type] += 1

        return atom_counts

    def add_disconnected_atom(self, atom):
        self.graph[atom] = []
        self.atoms[atom.nr] = atom

    def sort_by_nr(self):
        for atom in self.graph:
            self.graph[atom].sort(key=lambda x: x.nr)

    def make_dummy_bond(self, atom_1, atom_2, bond_nr, dummy=False):
        if dummy:
            bond_type = 'dummy'
        else:
            bond_type = 'single'

        if atom_1 in self.graph:
            self.graph[atom_1].append(atom_2)
        else:
            self.graph[atom_1] = [atom_2]

        if atom_2 in self.graph:
            self.graph[atom_2].append(atom_1)
        else:
            self.graph[atom_2] = [atom_1]

        bond = Bond(atom_1, atom_2, bond_type, bond_nr)

        atom_1.add_bond(bond)
        atom_2.add_bond(bond)

        self.bonds[bond_nr] = bond

        if atom_1 not in self.bond_lookup:
            self.bond_lookup[atom_1] = {}
        if atom_2 not in self.bond_lookup:
            self.bond_lookup[atom_2] = {}

        self.bond_lookup[atom_1][atom_2] = bond
        self.bond_lookup[atom_2][atom_1] = bond

    def make_bond(self, atom_1, atom_2, bond_nr):

        bond = Bond(atom_1, atom_2, 'single', bond_nr)

        electron_1 = None
        electron_2 = None

        orbital_1 = None
        orbital_2 = None

        for orbital in atom_1.valence_shell.orbitals:
            if orbital.electron_nr == 1 and not orbital.electrons[0].aromatic:
                orbital_1 = orbital
                electron_1 = orbital_1.electrons[0]
                break

        for orbital in atom_2.valence_shell.orbitals:
            if orbital.electron_nr == 1 and not orbital.electrons[0].aromatic:
                orbital_2 = orbital
                electron_2 = orbital_2.electrons[0]
                break

        orbital_1.add_electron(electron_2)
        orbital_2.add_electron(electron_1)

        orbital_1.set_bond(bond, 'sigma')
        orbital_2.set_bond(bond, 'sigma')

        atom_1.add_bond(bond)
        atom_2.add_bond(bond)

        self.bonds[bond_nr] = bond

        if atom_1 not in self.bond_lookup:
            self.bond_lookup[atom_1] = {}
        if atom_2 not in self.bond_lookup:
            self.bond_lookup[atom_2] = {}

        if atom_1 in self.graph:
            self.graph[atom_1].append(atom_2)
        else:
            self.graph[atom_1] = [atom_2]

        if atom_2 in self.graph:
            self.graph[atom_2].append(atom_1)
        else:
            self.graph[atom_2] = [atom_1]

        self.bond_lookup[atom_1][atom_2] = bond
        self.bond_lookup[atom_2][atom_1] = bond

        bond.electrons.append(electron_1)
        bond.electrons.append(electron_2)

        atom_1.set_neighbours(self)
        atom_2.set_neighbours(self)

    def add_bond(self, atom_1, atom_2, bond_type, bond_nr, chiral_symbol=None):
        if atom_1 in self.graph:
            self.graph[atom_1].append(atom_2)
        else:
            self.graph[atom_1] = [atom_2]

        if atom_2 in self.graph:
            self.graph[atom_2].append(atom_1)
        else:
            self.graph[atom_2] = [atom_1]

        bond = Bond(atom_1, atom_2, bond_type, bond_nr)

        bond.chiral_symbol = chiral_symbol

        atom_1.add_bond(bond)
        atom_2.add_bond(bond)

        self.bonds[bond_nr] = bond

        if atom_1 not in self.bond_lookup:
            self.bond_lookup[atom_1] = {}
        if atom_2 not in self.bond_lookup:
            self.bond_lookup[atom_2] = {}

        self.bond_lookup[atom_1][atom_2] = bond
        self.bond_lookup[atom_2][atom_1] = bond

    def reset_attribute(self, annotation, default=None):
        for atom in self.graph:
            atom.annotations.set_annotation(annotation, default)

    def add_attribute(self, annotation, default=None):
        for atom in self.graph:
            if annotation not in atom.annotations.annotations:
                atom.annotations.add_annotation(annotation, default)

        self.annotations.add((annotation, default))

    def reset_attributes(self, annotations, defaults=None, boolean=False):

        for atom in self.graph:
            for i, annotation in enumerate(annotations):
                if defaults:
                    default = defaults[i]
                elif boolean:
                    default = False
                else:
                    default = None

                atom.annotations.set_annotation(annotation, default)
        
    def add_attributes(self, annotations, defaults=None, boolean=False):
        if defaults:
            assert len(defaults) == len(annotations)

            for i, annotation in enumerate(annotations):
                self.annotations.add((annotation, defaults[i]))

        elif boolean:
            for i, annotation in enumerate(annotations):
                self.annotations.add((annotation, False))

        else:
            for i, annotation in enumerate(annotations):
                self.annotations.add((annotation, None))

        for atom in self.graph:
            for i, annotation in enumerate(annotations):
                if defaults:
                    default = defaults[i]
                elif boolean:
                    default = False
                else:
                    default = None

                if annotation not in atom.annotations.annotations:
                    atom.annotations.add_annotation(annotation, default)

    @staticmethod
    def set_attribute(atoms, annotation, value):
        for atom in atoms:
            atom.annotations.set_annotation(annotation, value)

    def get_atoms_of_type(self, atom_type):
        atoms = []
        for atom in self.atoms.values():
            if atom.type == atom_type:
                atoms.append(atom)

        return atoms

    def print_graph(self):
        pprint(self.graph)

    def print_atoms(self):
        pprint(self.atoms)

    def print_bonds(self):
        pprint(self.bonds)

    def find_pi_subgraph(self, prune=True):
        pi_subgraph = {}

        for bond in self.bonds.values():
            if bond.type == 'aromatic':

                # prune the subgraph as kekulisation can only occur in atoms
                # that have an unpaired electron

                unpaired_electrons_1 = 0
                unpaired_electrons_2 = 0

                if len(bond.aromatic_system.get_contributed_electrons(bond.atom_1)) == 1:
                    unpaired_electrons_1 += 1

                if len(bond.aromatic_system.get_contributed_electrons(bond.atom_2)) == 1:
                    unpaired_electrons_2 += 1

                if unpaired_electrons_1 and unpaired_electrons_2:

                    if bond.atom_1 not in pi_subgraph:
                        pi_subgraph[bond.atom_1] = []
                    if bond.atom_2 not in pi_subgraph:
                        pi_subgraph[bond.atom_2] = []

                    pi_subgraph[bond.atom_1].append(bond.atom_2)
                    pi_subgraph[bond.atom_2].append(bond.atom_1)

                elif not prune:

                    if bond.atom_1 not in pi_subgraph:
                        pi_subgraph[bond.atom_1] = []
                    if bond.atom_2 not in pi_subgraph:
                        pi_subgraph[bond.atom_2] = []

                    pi_subgraph[bond.atom_1].append(bond.atom_2)
                    pi_subgraph[bond.atom_2].append(bond.atom_1)

        return pi_subgraph

    def kekulise(self):

        kekule_structure = self.deepcopy()

        pruned = kekule_structure.find_pi_subgraph(prune=True)
        unpruned = kekule_structure.find_pi_subgraph(prune=False)

        aromatic_unmatched = set(unpruned.keys()) - set(pruned.keys())

        matching = Match.from_structure(kekule_structure)
        unmatched_nodes = matching.unmatched_nodes()
        if unmatched_nodes != 0:
            raise Exception("This structure cannot be kekulised!")
        else:
            double_bond_pairs = set()
            single_bond_pairs = set()

            for node in matching.nodes:
                double_bond_pair = tuple(sorted([node.atom, node.mate.atom], key=lambda x: x.nr))
                if double_bond_pair not in double_bond_pairs:
                    double_bond_pairs.add(double_bond_pair)

                for neighbour in node.neighbors:
                    if neighbour.index != node.mate.index:
                        single_bond_pair = tuple(sorted([node.atom, neighbour.atom], key=lambda x: x.nr))
                        if single_bond_pair not in single_bond_pairs:
                            single_bond_pairs.add(single_bond_pair)

            # heteroatoms containing lone pairs, sp2-hybridised carbons

            for atom in aromatic_unmatched:
                for neighbour in atom.neighbours:
                    if neighbour in atom.aromatic_system.atoms:
                        single_bond_pair = tuple(sorted([atom, neighbour], key=lambda x: x.nr))
                        if single_bond_pair not in single_bond_pairs:
                            single_bond_pairs.add(single_bond_pair)

        for aromatic_system in kekule_structure.aromatic_systems:
            aromatic_system.relocalise_electrons()

        for pair in double_bond_pairs:

            new_atom_1 = kekule_structure.atoms[pair[0].nr]
            new_atom_2 = kekule_structure.atoms[pair[1].nr]
            
            bond = kekule_structure.bond_lookup[new_atom_1][new_atom_2]
            bond.type = 'double'
            bond.aromatic = False
            
            bond.atom_1.aromatic = False
            bond.atom_2.aromatic = False

            bond.set_bond_summary()
            
            orbitals_1 = new_atom_1.get_orbitals('p')
            orbitals_2 = new_atom_2.get_orbitals('p')
            
            if orbitals_1 and orbitals_2:
                orbital_1 = orbitals_1[0]
                orbital_2 = orbitals_2[0]
                
                if not len(orbital_1.electrons) == 1 or not len(orbital_2.electrons) == 1:
                    raise KekulisationError(bond.aromatic_system.__repr__())

                orbital_1.add_electron(orbital_2.electrons[0])
                orbital_2.add_electron(orbital_1.electrons[0])

                orbital_1.set_bond(bond, 'pi')
                orbital_2.set_bond(bond, 'pi')

                bond.electrons.append(orbital_1.electrons[0])
                bond.electrons.append(orbital_2.electrons[0])

                bond.set_bond_summary()
                bond.aromatic_system = None

        for pair in single_bond_pairs:
            new_atom_1 = kekule_structure.atoms[pair[0].nr]
            new_atom_2 = kekule_structure.atoms[pair[1].nr]

            bond = kekule_structure.bond_lookup[new_atom_1][new_atom_2]
            bond.type = 'single'

            bond.aromatic = False
            bond.atom_1.aromatic = False
            bond.atom_2.aromatic = False
            bond.atom_1.pyrrole = False
            bond.atom_2.pyrrole = False
            bond.atom_1.furan = False
            bond.atom_2.furan = False
            bond.atom_1.thiophene = False
            bond.atom_2.thiophene = False

            bond.aromatic_system = None

            bond.set_bond_summary()

        kekule_structure.aromatic_systems = []
        kekule_structure.aromatic_cycles = []

        return kekule_structure

    def break_any_bond(self, atom_1, atom_2):
        """Remove edges from structure to break any bond

        atom_1: tuple of (str, int), with str atom type and int atom number
        atom_2: tuple of (str, int), with str atom type and int atom number

        """
        bonds_removed = 0
        while atom_2 in self.structure[atom_1]:
            self.structure[atom_1].remove(atom_2)
            bonds_removed += 1
        while atom_1 in self.structure[atom_2]:
            self.structure[atom_2].remove(atom_1)

        self.bonds[atom_1] -= bonds_removed
        self.bonds[atom_2] -= bonds_removed

    def split_disconnected_structures(self):
        """Return list of unconnected structures from structure


        Output:
        new_graphs: list of dicts of {node: [node, ->], ->}, with each dict a
            bidirectional graph that is not connected to the other graphs in
            the list
        """
        working_graph = self.deepcopy()

        new_graphs = []
        working_graph.make_bond_nr_dict()

        working_graph.remove_connectors()

        start_node = list(working_graph.graph.keys())[0]

        paths_collection = []
        paths = []

        while start_node:

            path = working_graph.find_a_path(start_node)
            paths.append(path)

            potential_start_nodes = working_graph.find_start_nodes(paths)

            try:

                start_node = potential_start_nodes[0]

            except IndexError:
                paths_collection.append(paths)
                paths = []
                potential_start_nodes = working_graph.find_new_start_node()

                try:
                    start_node = potential_start_nodes[0]

                except IndexError:
                    paths_collection.append(paths)
                    start_node = None

        counter = 0
        for paths in paths_collection:
            if paths:
                counter += 1
                new_graph = working_graph.put_paths_in_graph(paths)
                new_graphs.append(new_graph)

        # add back connectors

        for new_graph in new_graphs:
            for node in new_graph:

                new_graph[node] = self.graph[node]

        # Add lone atoms
        if working_graph.graph:
            for atom in working_graph.graph:
                new_graph = {atom: []}
                if new_graph not in new_graphs:
                    new_graphs.append(new_graph)

        new_structures = []

        for new_graph in new_graphs:
            new_structures.append(Structure(new_graph))

        for new_structure in new_structures:
            new_structure.infer_bonds()
            new_structure.set_atom_neighbours()
            new_structure.make_bond_lookup()
            new_structure.refresh_structure(find_cycles=True)

        return new_structures

    def infer_bonds(self):
        self.bonds = {}
        for atom in self.graph:
            for bond in atom.bonds:
                if bond.nr not in self.bonds:
                    self.bonds[bond.nr] = bond

    # ========================================================================
    # Auxillary functions
    # ========================================================================

    @staticmethod
    def put_paths_in_graph(paths):
        """Return single structure from bond paths

        Input:
        paths: list of [atom, ->], with each atom a tuple of (str, int), with
            str atom type and int atom number

        Output:
        rest_group_graph: dict of {atom: [atom, ->], ->}, with each atom a tuple
            of (str, int), str representing atom type, and int representing atom
            number. Example: ('O', 1). This graph represents the side chain of
            an amino acid
        """
        rest_group_graph = {}
        for path in paths:
            current_atom = path[0]
            if path[1:]:
                for atom in path[1:]:

                    next_atom = atom
                    if current_atom in rest_group_graph:
                        rest_group_graph[current_atom] += [next_atom]
                    else:
                        rest_group_graph[current_atom] = [next_atom]

                    if next_atom in rest_group_graph:
                        rest_group_graph[next_atom] += [current_atom]
                    else:
                        rest_group_graph[next_atom] = [current_atom]

                    current_atom = atom
            else:
                rest_group_graph = {current_atom: []}

        return rest_group_graph

    def get_bond_nr(self, atom_1, atom_2):
        bond_nr = self.structure[atom_1].count(atom_2)
        return bond_nr

    def make_bond_nr_dict(self):
        """
        """
        self.bond_nr_dict = {}
        for atom, neighbours in self.graph.items():
            self.bond_nr_dict[atom] = len(neighbours)

    def find_a_path(self, start_atom):
        """Return a list of linked atoms from a structure

        Input:
        structure: dict of {atom: [atom, ->], ->}, with each atom a tuple of
            (str, int), str representing atom type, and int representing atom
            number. Example: ('O', 1)
        start_atom: tuple of (str, int), with str atom type and int atom number
        bond_dict: dict of {atom: remaining_bonds, ->}, with atom tuple of
            (str, int), with str atom type and int atom number, and remaining bonds
            int

        Output:
        path: list of [atom, ->], where adjacent atoms are connected, and each atom
            is a tuple of (str, int), with str atom type and int atom number
        """

        current_atom = start_atom
        path = [current_atom]

        if len(self.graph[current_atom]) == 0:
            path = [current_atom]
            return path

        # keep trying to extend the path until there are no bonds to traverse
        while True:
            try:

                next_atom = self.graph[current_atom][0]

                path.append(next_atom)

                # remove traversed bond from structure
                self.graph[current_atom].remove(next_atom)
                if not next_atom == current_atom:
                    self.graph[next_atom].remove(current_atom)

                self.bond_nr_dict[current_atom] -= 1
                self.bond_nr_dict[next_atom] -= 1

                # remove atom from structure if no more untraversed bonds come
                # from it
                if not self.graph[current_atom]:
                    del self.graph[current_atom]

                current_atom = next_atom

                if not self.graph[current_atom]:
                    del self.graph[current_atom]

            except KeyError:
                break

        return path

    def find_start_nodes(self, paths):
        """Return atoms that still have outgoing bonds within an existing path

        Input:
        paths: list of [atom, ->], with each atom a tuple of (str, int), with
            str atom type and int atom number
        bond_dict: dict of {atom: remaining_bonds, ->}, with atom tuple of
            (str, int), with str atom type and int atom number, and remaining
            bonds int

        Output:
        start_atoms: list of [atom, ->], with each atom a tuple of (str, int),
            with str atom type and int atom number


        """

        start_atoms = []
        for path in paths:
            for atom in path:
                if self.bond_nr_dict[atom] > 0:
                    start_atoms.append(atom)

        return start_atoms

    def remove_connectors(self):
        """Remove nodes that only have incoming edges from graph

        Input:
        working_graph: dict of {node: [node, ->], ->}, representing a graph

        """

        for node in self.graph:
            for next_node in self.graph[node]:
                if next_node not in self.graph:
                    self.graph[node].remove(next_node)

    def find_new_start_node(self):
        """Return list of nodes that still have outgoing edges

        Input:
        edges_dict: dict of {node: int, ->}, with int representing the number of
            outgoing edges of that node

        Output:
        start_nodes: list of [node, ->], with each node an immutable object
        """
        start_nodes = []
        for atom in self.graph:
            if self.bond_nr_dict[atom] > 0:
                start_nodes.append(atom)

        return start_nodes


class Tree:
    def __init__(self, is_root):
        self.is_root = is_root
        self.parent = None
        self.children = []
