import copy
from pprint import pprint
import sys

from pikachu.chem.bond_properties import BOND_PROPERTIES
from pikachu.errors import SmilesError
from pikachu.chem.atom import Atom
from pikachu.chem.bond import Bond
from pikachu.chem.kekulisation import Match
from pikachu.chem.substructure_search import check_same_chirality, compare_all_matches, SubstructureMatch, find_substructures
from pikachu.chem.rings.ring_identification import check_five_ring, check_aromatic
import pikachu.chem.rings.find_cycles as find_cycles

sys.setrecursionlimit(100000)


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
        if bonds:
            self.bonds = bonds
            self.make_bond_lookup()
        else:
            self.bonds = {}

        if bond_lookup:
            self.bond_lookup = bond_lookup
        else:
            self.bond_lookup = {}

    def get_atoms(self):
        atoms = []
        for atom in self.graph:
            atoms.append(atom)

        return atoms

    def get_subtree(self, atom, parent_atom):
        pass

    def set_atoms(self):
        self.atoms = {}
        for atom in self.graph:
            self.atoms[atom.nr] = atom

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

    def colour_substructure_single(self, substructure, colour="hot pink", check_chiral_centres=True, check_bond_chirality=True):
        matches = self.find_substructures(substructure,
                                          check_chiral_centres=check_chiral_centres,
                                          check_chiral_double_bonds=check_bond_chirality)

        if matches:
            match = matches[0]
            for parent_atom in match.atoms.values():
                for atom in self.graph:
                    if atom == parent_atom:
                        atom.draw.colour = colour

    def colour_substructure_all(self, substructure, colour="hot pink", check_chiral_centres=True, check_bond_chirality=True):
        matches = self.find_substructures(substructure,
                                          check_chiral_centres=check_chiral_centres,
                                          check_chiral_double_bonds=check_bond_chirality)

        print(matches)

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

                # double bonds can only be chiral if they have exactly three bonds

                if len(bond.atom_1.bonds) == 3 and len(bond.atom_2.bonds) == 3:

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

                    for bond_1 in bond.atom_1.bonds:

                        if bond_1.type == 'single':
                            if bond.atom_1 == bond_1.atom_1:
                                if bond_1.chiral_symbol == '/':
                                    direction = 'up'
                                elif bond_1.chiral_symbol == '\\':
                                    direction = 'down'
                                else:
                                    direction = None

                                if not atom_1_1:
                                    atom_1_1 = bond_1.atom_2
                                    chiral_1_1 = direction
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
                            raise SmilesError('chiral double bond')
                        if chiral_2_2 == chiral_2_1:
                            raise SmilesError('chiral double bond')

                        if chiral_1_1:
                            first_atom = atom_1_1
                            first_other_atom = atom_1_2
                            first_chiral_symbol = chiral_1_1

                        else:
                            first_atom = atom_1_2
                            first_other_atom = atom_1_1
                            first_chiral_symbol = chiral_1_2

                        if chiral_2_1:
                            second_atom = atom_2_1
                            second_other_atom = atom_2_2
                            second_chiral_symbol = chiral_2_1
                        else:
                            second_atom = atom_2_2
                            second_other_atom = atom_2_1
                            second_chiral_symbol = chiral_2_2

                        bond.chiral_dict[first_atom] = {}
                        bond.chiral_dict[first_other_atom] = {}
                        bond.chiral_dict[second_atom] = {}
                        bond.chiral_dict[second_other_atom] = {}

                        if first_chiral_symbol == second_chiral_symbol:
                            bond.chiral_dict[first_atom][second_atom] = 'cis'
                            bond.chiral_dict[second_atom][first_atom] = 'cis'

                            bond.chiral_dict[first_other_atom][second_other_atom] = 'cis'
                            bond.chiral_dict[second_other_atom][first_other_atom] = 'cis'

                            bond.chiral_dict[first_atom][second_other_atom] = 'trans'
                            bond.chiral_dict[second_other_atom][first_atom] = 'trans'

                            bond.chiral_dict[first_other_atom][second_atom] = 'trans'
                            bond.chiral_dict[second_atom][first_other_atom] = 'trans'

                        else:
                            bond.chiral_dict[first_atom][second_atom] = 'trans'
                            bond.chiral_dict[second_atom][first_atom] = 'trans'

                            bond.chiral_dict[first_other_atom][second_other_atom] = 'trans'
                            bond.chiral_dict[second_other_atom][first_other_atom] = 'trans'

                            bond.chiral_dict[first_atom][second_other_atom] = 'cis'
                            bond.chiral_dict[second_other_atom][first_atom] = 'cis'

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

        if next_atom_nr == None:
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

        if next_bond_nr == None:
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
            if not atom_1 in self.bond_lookup:
                self.bond_lookup[atom_1] = {}
            if not atom_2 in self.bond_lookup:
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
        atom.add_shell_layout()

        for i, neighbour in enumerate(neighbours):
            next_bond_nr = self.find_next_bond_nr()
            self.make_bond(atom, neighbour, next_bond_nr)

        return atom

    def find_cycles(self):
        """
        Find all cycles in a structure and store them
        """
        self.cycles = find_cycles.Cycles(self)

    def promote_electrons_in_five_rings(self):
        """
        Promote electrons in aromatic five rings to p-orbitals
        """
        five_rings = self.cycles.find_five_membered()
        for five_ring in five_rings:
            aromatic, heteroatom = check_five_ring(five_ring)
            if aromatic:
                heteroatom.promote_lone_pair_to_p_orbital()
                if heteroatom.type == 'N':
                    heteroatom.pyrrole = True

    def find_aromatic_cycles(self):
        """
        Returns cycles that are aromatic

        Output
        ------
        aromatic_cycles: list of [[Atom, ->], ->], with each list of
            atoms the atoms that comprise an aromatic cycle

        """
        aromatic_cycles = []
        for cycle in self.cycles.unique_cycles:
            if check_aromatic(cycle):
                aromatic_cycles.append(cycle)

        return aromatic_cycles

    def find_aromatic_systems(self):
        """
        Return ring systems that are aromatic

        Output
        ------
        aromatic_ring_systems: list of [[Atom, ->], ->], with each list of
            atoms the atoms that comprise an aromatic ring system
        """

        ring_systems = self.cycles.find_cyclic_systems()
        aromatic_ring_systems = []
        for ring_system in ring_systems:
            if check_aromatic(ring_system):
                aromatic_ring_systems.append(ring_system)

        return aromatic_ring_systems

    def set_bonds_to_aromatic(self, aromatic_systems):
        for aromatic_system in aromatic_systems:
            for atom_1 in aromatic_system:
                atom_1.aromatic = True
                for atom_2 in aromatic_system:
                    if atom_1 in self.graph[atom_2]:
                        bond = self.bond_lookup[atom_1][atom_2]
                        if bond.type != 'aromatic':
                            bond.make_aromatic()

    def refine_structure(self):
        """
        """

        self.add_shells()
        self.add_hydrogens()
        self.add_shells()

        self.sort_by_nr()

        self.refine_p_bonds()
        self.hybridise_atoms()
        self.check_d_orbitals()

        self.refine_s_bonds()
        self.drop_electrons()
        self.set_atom_neighbours()

        self.make_lone_pairs()
        self.make_bond_lookup()
        self.find_cycles()
        self.promote_electrons_in_five_rings()
        aromatic_systems = self.find_aromatic_systems()
        aromatic_cycles = self.find_aromatic_cycles()
        self.set_bonds_to_aromatic(aromatic_systems)
        self.set_bonds_to_aromatic(aromatic_cycles)
        self.set_connectivities()
        self.set_atoms()

    def make_lone_pairs(self):
        for atom in self.graph:
            atom.make_lone_pairs()

    def check_d_orbitals(self):
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
            bond_nr = bond
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
                if not connectivity in connectivities:
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

    def check_chiral_centres(self, child, match):
        chirality_matches = True
        chiral_centres = child.get_steric_atoms()

        for chiral_centre in chiral_centres:
            parent_atom = match[chiral_centre]
            if parent_atom.chiral:

                chirality_matches = check_same_chirality(chiral_centre, parent_atom, match)
                if not chirality_matches:
                    print(chiral_centre, parent_atom, match)

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

    def get_unvisited_bonds(self, bond_to_visited, bond_to_matching_attempts, child_bond, bond_match, parent_atom):
        bonds = []
        for parent_bond in parent_atom.bonds:
            if parent_bond not in bond_match.values():
                if not bond_to_visited[parent_bond]:
                    bonds.append(parent_bond)
                elif child_bond not in bond_to_matching_attempts[parent_bond]:
                    bonds.append(parent_bond)

        return bonds

    def make_bond_match_dict(self):
        bond_to_bond = {}
        for bond in self.bonds.values():
            bond_to_bond[bond] = None

        return bond_to_bond

    def find_substructures_bonds(self, substructure, check_chiral_centres=True, check_chiral_double_bonds=True):
        matches = []

        if self.is_substructure_atom_composition(substructure):
            if self.is_substructure_bond_composition(substructure):
                if self.is_substructure_atom_connectivity(substructure):
                    pass

    def find_substructures(self, substructure, check_chiral_centres=True, check_chiral_double_bonds=True):
        matches = []
        if self.is_substructure_atom_composition(substructure):
            if self.is_substructure_atom_connectivity(substructure):
                matches = find_substructures(self, substructure)
                # matches = self.find_substructure_in_structure(substructure)

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
                if not atom.type in connectivities:
                    connectivities[atom.type] = {}
                connectivity = atom.connectivity
                if not connectivity in connectivities[atom.type]:
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

    def make_match_dict(self):
        match_dict = {}
        for atom in self.graph:
            if atom.type != 'H':
                match_dict[atom] = None

        return match_dict

    def compress_graph(self):
        compressed_graph = {}
        for atom in self.graph:
            if atom.type != 'H':
                new_atom = copy.deepcopy(atom)
                compressed_graph[new_atom] = copy.deepcopy(self.graph[atom])

    def reset_visited(self):
        for atom in self.graph:
            atom.visited = False
        for bond_nr in self.bonds:
            bond = self.bonds[bond_nr]
            bond.visited = False
            bond.visited_child_bonds = set()

    def remove_visited(self):
        for atom in self.graph:
            delattr(atom, 'visited')
        for bond_nr in self.bonds:
            bond = self.bonds[bond_nr]
            delattr(bond, 'visited')

    def traceback(self, placed_atoms_parent, match_dict, reverse_match_dict, parent_options, atoms_to_place,
                  child_bond_dict, bonds_visited_through_ring_closure):

        new_parent_candidate = None
        new_child_candidate = None
        previous_child_candidate = None

        for i in range(len(placed_atoms_parent) - 1, -1, -1):

            current_atom = placed_atoms_parent[i]

            try:
                for option in parent_options[current_atom]:

                    if not option.visited:
                        new_parent_candidate = current_atom
                        new_child_candidate = reverse_match_dict[current_atom]
                        break

            except KeyError:
                pass

            if new_parent_candidate:
                break
            else:
                del placed_atoms_parent[i]
                child_current = reverse_match_dict[current_atom]
                atoms_to_place.add(child_current)

                match_dict[child_current] = None

                child_bond_dict[child_current] += 1

                for bond in child_current.bonds:
                    for atom in bond.neighbours:
                        if atom != child_current:
                            if atom in match_dict:
                                child_bond_dict[atom] += 1

                del reverse_match_dict[current_atom]

        return new_child_candidate, new_parent_candidate

    def find_substructure_in_structure(self, child):

        atom_connectivities_child = child.get_connectivities()
        atom_connectivities_parent = self.get_substructure_connectivities(atom_connectivities_child)
        # Sort based on the complexity of the connectivity

        connectivities = sorted(list(atom_connectivities_child.keys()),
                                key=lambda x: len(set(x)), reverse=True)

        matches = []

        starting_connectivity = connectivities[0]

        starting_atom = atom_connectivities_child[starting_connectivity][0]

        seeds = []

        for atom in atom_connectivities_parent[starting_connectivity]:
            if atom.type == starting_atom.type:
                seeds.append(atom)

        for seed in seeds:

            self.reset_visited()

            child_bond_dict = child.make_bond_dict()

            match_active = True
            match = child.make_match_dict()
            reverse_match = {}
            match[starting_atom] = seed
            reverse_match[seed] = starting_atom
            visited_through_ring_closure = []

            atoms_to_place = set(copy.copy(list(child.graph.keys())))
            for atom in copy.deepcopy(atoms_to_place):
                if atom.type == 'H':
                    atoms_to_place.remove(atom)

            atoms_to_place.remove(starting_atom)

            placed_atoms_child = {starting_atom}
            placed_atoms_parent = [seed]
            parent_options = {}

            current_atom_child = starting_atom
            current_atom_parent = seed

            while atoms_to_place and match_active:

                child_neighbours = current_atom_child.neighbours
                next_atom_child = None

                for neighbour in child_neighbours:
                    if neighbour in atoms_to_place and neighbour.type != 'H':
                        next_atom_child = neighbour

                        break

                if next_atom_child:

                    child_bond = child.bond_lookup[current_atom_child][next_atom_child]
                    next_atom_parent_options = []
                    parent_neighbours = current_atom_parent.neighbours
                    for neighbour in parent_neighbours:
                        parent_bond = self.bond_lookup[current_atom_parent][neighbour]

                        if neighbour.type == next_atom_child.type and child_bond.type == parent_bond.type:
                            bond = self.bond_lookup[current_atom_parent][neighbour]
                            if child_bond not in bond.visited_child_bonds and neighbour not in match.values():
                                if neighbour.potential_same_connectivity(next_atom_child.connectivity):
                                    next_atom_parent_options.append(neighbour)

                    if not current_atom_parent in parent_options:
                        parent_options[current_atom_parent] = []

                    if next_atom_parent_options:

                        next_atom_parent_options = sorted(next_atom_parent_options,
                                                          key=lambda x: len(set(x.connectivity)), reverse=True)
                        for option in next_atom_parent_options:
                            bond = self.bond_lookup[current_atom_parent][option]
                            if bond not in parent_options[current_atom_parent]:
                                parent_options[current_atom_parent].append(bond)

                        next_atom_parent = next_atom_parent_options[0]
                        bond = self.bond_lookup[current_atom_parent][next_atom_parent]
                        bond.visited = True
                        bond.visited_child_bonds.add(child_bond)

                        match[next_atom_child] = next_atom_parent
                        reverse_match[next_atom_parent] = next_atom_child
                        placed_atoms_parent.append(next_atom_parent)

                        child_bond_dict[current_atom_child] -= 1
                        child_bond_dict[next_atom_child] -= 1

                        for neighbour in next_atom_child.neighbours:
                            if neighbour.type != 'H' and neighbour != current_atom_child \
                                    and neighbour not in atoms_to_place:
                                bond = self.bond_lookup[match[neighbour]][next_atom_parent]
                                if not bond.visited:
                                    visited_through_ring_closure.append(bond)
                                    child_bond_dict[neighbour] -= 1
                                    child_bond_dict[next_atom_child] -= 1
                                    bond.visited = True
                                    bond.visited_child_bonds.add(bond)


                        current_atom_child = next_atom_child
                        current_atom_parent = next_atom_parent
                        atoms_to_place.remove(current_atom_child)

                    else:
                        new_child_candidate, new_parent_candidate = self.traceback(placed_atoms_parent, match,
                                                                                   reverse_match, parent_options,
                                                                                   atoms_to_place, child_bond_dict,
                                                                                   visited_through_ring_closure)
                        if not new_child_candidate:
                            match_active = False
                        else:
                            current_atom_child = new_child_candidate
                            current_atom_parent = new_parent_candidate

                else:
                    if atoms_to_place:

                        for atom in match.keys():
                            if atom != current_atom_child:
                                if child_bond_dict[atom] > 0:
                                    current_atom_child = atom
                                    current_atom_parent = match[atom]
                                    break
                    else:
                        pass

            if match_active:
                matches.append(match)

        compare_all_matches(matches)

        self.remove_visited()

        return matches

    def drop_electrons(self):
        for atom in self.graph:
            atom.drop_electrons()

    def add_shells(self):
        for atom in self.graph:
            if not atom.shells:
                atom.add_shell_layout()

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

    def refine_p_bonds(self):
        for bond_nr in self.bonds:
            bond = self.bonds[bond_nr]
            if bond.type != 'single':
                bond.combine_p_orbitals()

    def refine_s_bonds(self):
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

        if not atom_1 in self.bond_lookup:
            self.bond_lookup[atom_1] = {}
        if not atom_2 in self.bond_lookup:
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
            if atom_1.valence_shell.orbitals[orbital].electron_nr == 1 and not \
            atom_1.valence_shell.orbitals[orbital].electrons[0].aromatic:
                orbital_1 = atom_1.valence_shell.orbitals[orbital]
                electron_1 = orbital_1.electrons[0]
                break

        for orbital in atom_2.valence_shell.orbitals:
            if atom_2.valence_shell.orbitals[orbital].electron_nr == 1 and not \
            atom_2.valence_shell.orbitals[orbital].electrons[0].aromatic:
                orbital_2 = atom_2.valence_shell.orbitals[orbital]
                electron_2 = orbital_2.electrons[0]
                break

        orbital_1.add_electron(electron_2)
        orbital_2.add_electron(electron_1)

        orbital_1.set_bond(bond, 'sigma')
        orbital_2.set_bond(bond, 'sigma')

        atom_1.add_bond(bond)
        atom_2.add_bond(bond)

        self.bonds[bond_nr] = bond

        if not atom_1 in self.bond_lookup:
            self.bond_lookup[atom_1] = {}
        if not atom_2 in self.bond_lookup:
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

        if not atom_1 in self.bond_lookup:
            self.bond_lookup[atom_1] = {}
        if not atom_2 in self.bond_lookup:
            self.bond_lookup[atom_2] = {}

        self.bond_lookup[atom_1][atom_2] = bond
        self.bond_lookup[atom_2][atom_1] = bond

    def print_graph(self):
        pprint(self.graph)

    def print_atoms(self):
        pprint(self.atoms)

    def print_bonds(self):
        pprint(self.bonds)

    def find_aromatic_bonds(self):
        aromatic_bonds = []
        for bond in self.bonds:
            if self.bonds[bond].type == 'aromatic':
                aromatic_bonds.append(self.bonds[bond])

        return aromatic_bonds

    def find_aromatic_atoms(self):
        aromatic_bonds = self.find_aromatic_bonds()
        aromatic_atoms = []
        for atom in self.graph:
            for bond in atom.bonds:
                if bond in aromatic_bonds:
                    aromatic_atoms.append(atom)

        return aromatic_atoms

    def find_aromatic_graphs(self):
        aromatic_atoms = self.find_aromatic_atoms()
        aromatic_bonds = self.find_aromatic_bonds()
        aromatic_graph = {}

        for atom in copy.deepcopy(self.graph):
            if atom in aromatic_atoms:
                aromatic_graph[copy.deepcopy(atom)] = copy.deepcopy(self.graph[atom])

        for atom in aromatic_graph:
            for bond in atom.bonds:
                if bond.type != 'aromatic':
                    atom.bonds.remove(bond)

        for atom_1 in aromatic_graph:

            for atom_2 in aromatic_graph[atom_1]:
                aromatic = False
                for bond in atom_1.bonds:
                    if bond.atom_1 == atom_2 or bond.atom_2 == atom_2:
                        aromatic = True

                if not aromatic:
                    aromatic_graph[atom_1].remove(atom_2)

        aromatic_structure = Structure(aromatic_graph, aromatic_bonds)

        aromatic_structures = aromatic_structure.split_disconnected_structures()

        return aromatic_structures

    def find_pi_subgraph(self, prune=True):
        pi_subgraph = {}

        for bond_nr, bond in self.bonds.items():
            if bond.type == 'aromatic':

                # prune the subgraph as kekulisation can only occur in atoms
                # #that have an unpaired electron

                unpaired_electrons_1 = 0
                unpaired_electrons_2 = 0

                for orbital_name, orbital in bond.atom_1.valence_shell.orbitals.items():
                    if len(orbital.electrons) == 1:
                        unpaired_electrons_1 += 1

                for orbital_name, orbital in bond.atom_2.valence_shell.orbitals.items():
                    if len(orbital.electrons) == 1:
                        unpaired_electrons_2 += 1

                if unpaired_electrons_1 and unpaired_electrons_2:

                    if not bond.atom_1 in pi_subgraph:
                        pi_subgraph[bond.atom_1] = []
                    if not bond.atom_2 in pi_subgraph:
                        pi_subgraph[bond.atom_2] = []

                    pi_subgraph[bond.atom_1].append(bond.atom_2)
                    pi_subgraph[bond.atom_2].append(bond.atom_1)

                elif not prune:

                    if not bond.atom_1 in pi_subgraph:
                        pi_subgraph[bond.atom_1] = []
                    if not bond.atom_2 in pi_subgraph:
                        pi_subgraph[bond.atom_2] = []

                    pi_subgraph[bond.atom_1].append(bond.atom_2)
                    pi_subgraph[bond.atom_2].append(bond.atom_1)

        return pi_subgraph

    def kekulise(self):

        pruned = self.find_pi_subgraph(prune=True)
        unpruned = self.find_pi_subgraph(prune=False)

        aromatic_unmatched = set(unpruned.keys()) - set(pruned.keys())

        matching = Match.from_structure(self)
        unmatched_nodes = matching.unmatched_nodes()
        if unmatched_nodes != 0:
            raise Exception("This structure cannot be kekulised!")
        else:
            double_bond_pairs = set()
            single_bond_pairs = set()

            for node in matching.nodes:
                double_bond_pair = tuple(sorted([node.atom, node.mate.atom], key=lambda x: x.nr))
                if not double_bond_pair in double_bond_pairs:
                    double_bond_pairs.add(double_bond_pair)

                for neighbour in node.neighbors:
                    if neighbour.index != node.mate.index:
                        single_bond_pair = tuple(sorted([node.atom, neighbour.atom], key=lambda x: x.nr))
                        if not single_bond_pair in single_bond_pairs:
                            single_bond_pairs.add(single_bond_pair)

            for atom in aromatic_unmatched:
                for neighbour in atom.neighbours:
                    if neighbour in pruned:
                        single_bond_pair = tuple(sorted([atom, neighbour], key=lambda x: x.nr))
                        if not single_bond_pair in single_bond_pairs:
                            single_bond_pairs.add(single_bond_pair)

        kekule_structure = copy.deepcopy(self)

        for pair in double_bond_pairs:
            bond = kekule_structure.bond_lookup[pair[0]][pair[1]]
            bond.type = 'double'
            bond.aromatic = False
            bond.atom_1.aromatic = False
            bond.atom_2.aromatic = False

        for pair in single_bond_pairs:
            bond = kekule_structure.bond_lookup[pair[0]][pair[1]]
            bond.type = 'single'
            bond.aromatic = False
            bond.atom_1.aromatic = False
            bond.atom_2.aromatic = False
            bond.atom_1.pyrrole = False
            bond.atom_2.pyrrole = False

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

    def remove_an_atom(self, atom):
        """Remove a node from structure to remove an atom

        atom_1: tuple of (str, int), with str atom type and int atom number
        """
        del self.strucrure[atom]

    def remove_isolated_atoms(self):
        """Remove nodes without outgoing edges from structure

        atom_1: tuple of (str, int), with str atom type and int atom number
        """
        for atom in list(self.graph.keys()):
            if not self.graph[atom]:
                del self.graph[atom]

    def split_disconnected_structures(self):
        """Return list of unconnected structures from structure


        Output:
        new_graphs: list of dicts of {node: [node, ->], ->}, with each dict a
            bidirectional graph that is not connected to the other graphs in
            the list
        """
        working_graph = copy.deepcopy(self)
        working_graph.refresh_structure()
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
            new_structure.refresh_structure()

        return new_structures

    def infer_bonds(self):
        self.bonds = {}
        for atom in self.graph:
            for bond in atom.bonds:
                if not bond.nr in self.bonds:
                    self.bonds[bond.nr] = bond

    def find_lowest_atom_nr(self):
        lowest_atom_nr = 4242424242424242

        for atom in self.graph:
            if atom.nr < lowest_atom_nr:
                lowest_atom_nr = atom.nr

        return lowest_atom_nr

    def sort_substructures_by_nr(self):
        disconnected_structures = self.split_disconnected_structures()

        index_dict = {}

        for i, structure in enumerate(disconnected_structures):
            lowest_atom_nr = structure.find_lowest_atom_nr()
            index_dict[lowest_atom_nr] = i

        new_disconnected_structures = []

        sorted_indices = list(index_dict.keys())
        sorted_indices.sort()

        for index in sorted_indices:
            new_disconnected_structures += [disconnected_structures[index_dict[index]]]

        return new_disconnected_structures

    def make_sorted_copy(self):
        copy = Structure(self.structure)
        for atom in copy.structure:
            copy.structure[atom].sort()

        return copy

    def compare_structures(self, structure_2):
        sorted_copy_group_1 = self.make_sorted_copy()
        sorted_copy_group_2 = structure_2.make_sorted_copy()

        if sorted_copy_group_1.structure == sorted_copy_group_2.structure:
            return True
        else:
            return False

    def label_smiles_nrs(self, K_to_S_dict):
        for atom in self.structure:
            atom.set_smiles_nr(K_to_S_dict)

    def get_true_neighbour(self, atom_1, atom_2):
        bond = self.bond_lookup[atom_1][atom_2]
        true_neighbour = None
        for atom in self.graph[atom_1]:
            if atom in self.graph[atom_2] and self.bond_lookup[atom_1][atom].type != 'dummy' and \
                    self.bond_lookup[atom_2][atom].type != 'dummy':
                true_neighbour = atom
                break

        return true_neighbour

    # ========================================================================
    # Auxillary functions
    # ========================================================================

    def make_H_string(self, atom, valence_dict):
        H_nr = valence_dict[atom]['max'] - valence_dict[atom]['current']
        if H_nr > 1:
            H_string = "[%sH%d]" % (atom.type, H_nr)
        elif H_nr == 1:
            H_string = "[%sH]" % (atom.type)
        else:
            H_string = "[%s]" % (atom.type)

        return H_string

    def make_nr_string(self, bond):
        bond_nr, bond_type = bond
        if bond_type == 1:
            string = ''
        elif bond_type == 2:
            string = '='
        elif bond_type == 3:
            string = '#'

        if bond_nr > 9:
            bond_string = '%%%d' % bond_nr
        else:
            bond_string = str(bond_nr)

        return ''.join([string, bond_string])

    def make_valence_dict(self):
        valence_dict = {}
        for atom in self.structure:
            valence_dict[atom] = {}
            valence_dict[atom]['current'] = len(self.structure[atom])
            valence_dict[atom]['max'] = atom.get_valence(len(self.structure[atom]))

        return valence_dict

    def record_bonds(self):
        bond_record = {}

        counter = 1

        self.bond_record = {}
        for atom_1 in self.structure:
            self.bond_record[atom_1] = []

            for atom_2 in self.structure[atom_1]:
                bond_type = self.structure[atom_1].count(atom_2)

                if atom_1.nr > atom_2.nr:
                    pair = (atom_2, atom_1)
                else:
                    pair = (atom_1, atom_2)

                if not pair in bond_record:
                    bond_nr = counter
                    bond_record[pair] = bond_nr
                    counter += 1
                else:
                    bond_nr = bond_record[pair]

                if not (bond_nr, bond_type) in self.bond_record[atom_1]:
                    self.bond_record[atom_1] += [(bond_nr, bond_type)]
            self.bond_record[atom_1].sort()

    def put_paths_in_graph(self, paths):
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
        for atom in self.graph:
            self.bond_nr_dict[atom] = len(atom.bonds)

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

        nodes = list(self.graph.keys())

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

