#!/usr/bin/env python
import pikachu
import find_cycles
import copy
from pprint import pprint

def get_cyclic_label(cycle_nr):
    if cycle_nr > 9:
        return '%' + str(cycle_nr)
    else:
        return str(cycle_nr)

class GraphToSmiles:
    def __init__(self, structure):
        self.structure = copy.deepcopy(structure)

        self.remove_hydrogens()

        self.add_representations()
        pprint(self.representations)

        self.find_branch_points()
        self.find_terminal_nodes()
        self.find_cycles()

        self.collapse_graph()
        self.make_simplified_graph_from_collapsed()

        pprint(self.collapsed_graph)
        self.map_node_pair_to_edges()
        pprint(self.node_pair_to_edge)
        pprint(self.node_pair_to_atoms)

        print('branch points', self.branch_points)
        print('terminal nodes', self.terminal_nodes)
        print('cycles', self.cycles)

        self.make_smiles_components()
        print(self.components)

    def make_smiles_components(self):
        self.components = []
        self.atom_to_index = {}
        is_branched = {}
        cycle_nr = 0

        for atom in self.structure:
            self.atom_to_index[atom] = None
            is_branched[atom] = False

        if self.terminal_nodes:
            first_atom = list(self.terminal_nodes)[0]
        elif self.branch_points:
            first_atom = list(self.branch_points)[0]
        else: #this happens when the entire graph is cyclic
            first_atom = list(self.structure.graph.keys())[0]

        self.atom_to_index[first_atom] = 0

        current_atom = first_atom
        self.components.append(self.representations[current_atom])

        working_graph = copy.deepcopy(self.structure)

        while working_graph:
            cyclic = False
            if working_graph[current_atom]:
                next_atom = working_graph[current_atom][0]
                if self.atom_to_index[next_atom] != None:
                    cyclic = True
                    cycle_nr += 1
                    cyclic_label = get_cyclic_label(cycle_nr)

                if cyclic:
                    cyclic_label_idx_1 = self.atom_to_index[next_atom]
                    cyclic_label_idx_2 = self.atom_to_index[current_atom]


                if not cyclic:
                    is_branched[current_atom] = True





            else:
                del working_graph[current_atom]



    def make_smiles_components(self):
        self.components = []
        self.atom_order = []

        index = 0
        used_nodes = set()

        self.node_to_branch_nr = self.make_branch_lookup()
        self.node_to_index = self.make_index_dict()

        cycle_nr = 0

        while self.node_pair_to_edge:

            last_node = None

            for node in self.node_pair_to_edge:
                if node in used_nodes:
                    last_node = node
                    break

            if not last_node:
                last_node = list(self.node_pair_to_edge.keys())[0]

            if not self.components:
                self.components = [self.representations[last_node]]
                self.node_to_index[last_node] = 0

            possible_next_step = True

            while possible_next_step:
                next_node = None

                self.node_to_branch_nr[last_node] += 1
                used_nodes.add(last_node)

                if last_node in self.node_pair_to_edge:
                    next_node = list(self.node_pair_to_edge[last_node].keys())[0]
                    edges = self.node_pair_to_edge[last_node][next_node]

                    if last_node == next_node:
                        cyclic_single = True
                        cyclic_multiple = False
                        cyclic_nr += 1

                    elif len(edges) > 1:
                        cyclic_single = False
                        cyclic_multiple = True
                        cyclic_nr += 1

                    else:
                        cyclic_single = False
                        cyclic_multiple = False

                    break_point = self.node_to_index[last_node]

                    if self.node_to_branch_nr[last_node] == 1:
                        # This is the first branch we are making for this node

                        if cyclic_single:


                        if edge:
                            insert = [edge, self.representations[next_node]]
                        else:
                            insert = [self.representations[next_node]]

                    else:

                        if edge:
                            insert = ['(', edge, self.representations[next_node], ')']
                        else:
                            insert = ['(', self.representations[next_node], ')']

                    self.add_insert(insert, break_point)
                    self.node_to_index[next_node] = self.node_to_index[last_node] + len(insert)

                    used_nodes.add(next_node)

                else:
                    possible_next_step = False
                    continue

                del self.node_pair_to_edge[last_node][next_node]
                del self.node_pair_to_edge[next_node][last_node]

                if not self.node_pair_to_edge[last_node]:
                    del self.node_pair_to_edge[last_node]

                if not self.node_pair_to_edge[next_node]:
                    del self.node_pair_to_edge[next_node]

    def collapse_graph(self):
        self.collapsed_graph = {}

        if self.terminal_nodes:
            first_node = list(self.terminal_nodes)[0]
        elif self.branch_points:
            first_node = list(self.branch_points)[0]
        else: #this happens when the entire graph is cyclic
            first_node = list(self.structure.graph.keys())[0]
            self.collapsed_graph[first_node] = []

        for node in self.branch_points:
            self.collapsed_graph[node] = []

        for node in self.terminal_nodes:
            self.collapsed_graph[node] = []



        last_collapsed_node = first_node
        current_node = first_node
        previous_node = None

        edge = []
        edge_atoms = []

        while current_node:
            if current_node in self.structure.graph:

                next_node = self.structure.graph[current_node][0]
                bond_symbol = self.get_bond_symbol(current_node, next_node)

                edge.append(bond_symbol)

                if next_node in self.collapsed_graph:
                    edge = ''.join(edge)

                    self.collapsed_graph[last_collapsed_node].append((next_node, edge, edge_atoms))

                    if next_node != last_collapsed_node: #only do this when molecule is not cyclic
                        self.collapsed_graph[next_node].append((last_collapsed_node, edge[::-1], edge_atoms[::-1]))

                    edge = []
                    edge_atoms = []

                    self.remove_bonds_and_nodes(current_node, next_node)

                    last_collapsed_node = next_node


                else:

                    edge.append(self.representations[next_node])
                    edge_atoms.append(next_node)

                    self.remove_bonds_and_nodes(current_node, next_node)

                current_node = next_node
            else:

                candidate_node = None

                for node in self.collapsed_graph:
                    if node in self.structure.graph:
                        candidate_node = node
                        break

                if not candidate_node:
                    current_node = None

                else:
                    current_node = candidate_node
                    last_collapsed_node = current_node


    def get_neighbour_nr(self, node):
        neighbours = self.structure.graph[node]
        neighbour_nr = 0
        for neighbour in neighbours:
            if neighbour.type != 'H':
                neighbour_nr += 1

        return neighbour_nr

    def remove_hydrogens(self):
        for atom in list(self.structure.graph.keys()):
            if atom.type == 'H' and atom.charge == 0:
                neighbour = atom.neighbours[0]
                if neighbour.type in ['B', 'C', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I'] \
                   and neighbour.charge == 0 and not neighbour.chiral:
                    bond = atom.bonds[0]
                    del self.structure.bonds[bond.nr]
                    del self.structure.bond_lookup[atom][neighbour]
                    del self.structure.bond_lookup[neighbour][atom]
                    
                    self.structure.graph[neighbour].remove(atom)
                    self.structure.graph[atom].remove(neighbour)

                    if not self.structure.graph[neighbour]:
                        del self.structure.graph[neighbour]

                    if not self.structure.graph[atom]:
                        del self.structure.graph[atom]

    def get_bond_symbol(self, last_node, next_node):
        bond_type = self.structure.bond_lookup[last_node][next_node].type
        symbol = pikachu.BOND_PROPERTIES.bond_type_to_symbol[bond_type]
        return symbol

    def remove_explicit_hydrogen(self, hydrogen):
        neighbour = hydrogen.neighbours[0]

        bond = hydrogen.bonds[0]

        del self.structure.bonds[bond.nr]
        del self.structure.bond_lookup[hydrogen][neighbour]
        del self.structure.bond_lookup[neighbour][hydrogen]

        self.structure.graph[neighbour].remove(hydrogen)

        del self.structure.graph[hydrogen]


    def remove_bonds_and_nodes(self, last_node, next_node):
        bond = self.structure.bond_lookup[last_node][next_node]
        
        del self.structure.bonds[bond.nr]
        del self.structure.bond_lookup[last_node][next_node]
        del self.structure.bond_lookup[next_node][last_node]

        if next_node in self.structure.graph[last_node]:
            self.structure.graph[last_node].remove(next_node)
        if last_node in self.structure.graph[next_node]:
            self.structure.graph[next_node].remove(last_node)

        if not self.structure.graph[last_node]:
            del self.structure.graph[last_node]

        if not self.structure.graph[next_node]:
            del self.structure.graph[next_node]

    def make_simplified_graph_from_collapsed(self):
        self.simplified_graph = {}
        for atom, next_atom_and_edges in self.collapsed_graph.items():
            self.simplified_graph[atom] = []
            for next_atom, next_edge, edge_atoms in next_atom_and_edges:
                self.simplified_graph[atom].append(next_atom)

    def make_branch_lookup(self):
        branch_lookup = {}
        for node in self.simplified_graph:
            branch_lookup[node] = 0

        return branch_lookup

    def make_index_dict(self):
        index_dict = {}
        
        for node in self.simplified_graph:
            index_dict[node] = None

        return index_dict

    def map_node_pair_to_edges(self):
        self.node_pair_to_edge = {}
        self.node_pair_to_atoms = {}

        for atom, neighbours_and_edges in self.collapsed_graph.items():
            if atom.type != 'H':
                self.node_pair_to_edge[atom] = {}
                self.node_pair_to_atoms[atom] = {}
                for neighbour, edge, atoms in neighbours_and_edges:
                    if neighbour.type != 'H':
                        if not neighbour in self.node_pair_to_edge[atom]:
                            self.node_pair_to_edge[atom][neighbour] = []
                            self.node_pair_to_atoms[atom][neighbour] = []

                        self.node_pair_to_edge[atom][neighbour].append(edge)
                        self.node_pair_to_atoms[atom][neighbour].append(atoms)

        print('nodes to edges')
        pprint(self.node_pair_to_edge)

        for atom, neighbour_to_edge in self.node_pair_to_edge.items():
            for neighbour, edges in neighbour_to_edge.items():
                zipped_edges = zip(self.node_pair_to_edge[atom][neighbour],
                                   self.node_pair_to_atoms[atom][neighbour])
                sorted_edges_and_atoms = sorted(zipped_edges, key = lambda x: x[0], reverse = True)
                sorted_edges, sorted_atoms = [list(tuple) for tuple in zip(*sorted_edges_and_atoms)]


                self.node_pair_to_edge[atom][neighbour] = sorted_edges
                self.node_pair_to_atoms[atom][neighbour] = sorted_atoms

    def adjust_atom_indices(self, insert, break_point):
        index_jump = len(insert)

        for atom, index in self.atom_to_index.items():
            if index != None and index > break_point:
                self.atom_to_index[atom] += index_jump

    def add_insert(self, insert, break_point):
        half_1 = self.components[:break_point + 1]
        half_2 = self.components[break_point + 1:]

        self.adjust_node_indices(insert, break_point)

        self.components = half_1 + insert + half_2               



    def add_representations(self):
        self.representations = {}
        for atom in self.structure.graph:
            if atom.type != 'H':
                if atom.aromatic:
                    self.representations[atom] = atom.type.lower()
                else:
                    self.representations[atom] = atom.type

        self.add_chiral_placeholders()
        self.add_hydrogen_representations()
        self.add_charge_representations()

    def add_chiral_placeholders(self):
        for atom in self.structure.graph:
            if atom.chiral:
                self.representations[atom] = f'[{self.representations[atom]}X]'

    def add_charge_representations(self):
        for atom in self.structure.graph:
            if atom.charge != 0:
                if charge == 1:
                    charge_string = '+'
                elif charge == -1:
                    charge_string = '-'
                elif charge_string > 1:
                    charge_string = '+%d' % charge
                elif charge_string < -1:
                    charge_string = '-%d' % charge
                    
                representation = self.representations[atom]
                if representation[-1] == ']':
                    self.representations[atom] = representation[:-1] + charge_string + ']'
                else:
                    self.representations[atom] = '[' + representation + charge_string + ']'

    def count_hydrogens(self):
        hydrogen_counts = {}
        atom_to_hydrogens = {}

        for atom, neighbours in self.structure.graph.items():
            hydrogen_counts[atom] = 0
            atom_to_hydrogens[atom] = []
            for neighbour in neighbours:
                if neighbour.type == 'H':
                    hydrogen_counts[atom] += 1
                    atom_to_hydrogens[atom].append(neighbour)

        return hydrogen_counts, atom_to_hydrogens

    def add_hydrogen_representations(self):
        hydrogen_counts, atom_to_hydrogens = self.count_hydrogens()
        
        for atom, count in hydrogen_counts.items():
            if count > 0:
                if count > 1:
                    hydrogen_string = f'H{count}'
                elif count == 1:
                    hydrogen_string = 'H'



                representation = self.representations[atom]
                if representation[-1] == ']':
                    self.representations[atom] = representation[:-1] + hydrogen_string + ']'
                else:
                    self.representations[atom] = '[' + representation + hydrogen_string + ']'

                explicit_hydrogens = atom_to_hydrogens[atom]
                for hydrogen in explicit_hydrogens:
                    self.remove_explicit_hydrogen(hydrogen)

    def find_cycles(self):
        cycles = find_cycles.Cycles(self.structure)

        self.cycles = cycles.find_minimal_cycles()

    def find_branch_points(self):
        self.branch_points = set()
        for node in self.structure.graph:
            if self.is_branch_point(node):
                self.branch_points.add(node)

    def find_terminal_nodes(self):
        self.terminal_nodes = set()
        for node in self.structure.graph:
            if self.is_terminal_node(node):
                self.terminal_nodes.add(node)

    def find_cyclic_pairs(self):
        self.cyclic_pairs = set()
        self.cyclic_nodes = set()
        for cycle in self.cycles:
            for i, atom in enumerate(cycle):
                previous_atom = cycle[i - 1]
                self.cyclic_nodes.add(previous_atom)
                self.cyclic_nodes.add(atom)
                pair = tuple(sorted((previous_atom, atom), key = lambda x: x.nr))
                self.cyclic_pairs.add(pair)

    def is_terminal_node(self, node):
        if self.get_neighbour_nr(node) == 1:
            return True
        else:
            return False

    def is_branch_point(self, node):
        if self.get_neighbour_nr(node) >= 3:
            return True
        else:
            return False




if __name__ == "__main__":
    smiles = 'c1ccc2c(c1)c(c[nH]2)C[C@@H](C(=O)O)N'
 #   smiles = 'C1=CC=C2C(=C1)C(=CN2)C[C@@H](C(=O)O)N'
 #   smiles = 'c1(C)c(C)c(C)c(C)c(C)c1(C)'
 #   smiles = 'c1(c(c(c(c(c1C)C)C)C)C)C'
    smiles = 'c12ccccc1cccc2'
 #   smiles = 'c1ccccc1'
 #   smiles = 'C1=CC=C(C=C1)C[C@@H](C(=O)O)N'
    structure = pikachu.Smiles(smiles).smiles_to_structure()
    collapsed_structure = GraphToSmiles(structure)
 #   pprint(collapsed_structure.collapsed_graph)
 #   pprint(collapsed_structure.simplified_graph)
  #  pprint(collapsed_structure.representations)
 #   pprint(collapsed_structure.cycles)
#    pprint(collapsed_structure.simplified_graph)
#    pprint(collapsed_structure.representations)
#    pprint(collapsed_structure.node_to_edge_dict)
    
    

        
