#!/usr/bin/env python
import pikachu
import find_cycles
import copy
from pprint import pprint

class CollapsedGraph:
    def __init__(self, structure):
        self.structure = copy.deepcopy(structure)
        self.remove_hydrogens()

        
        self.collapsed_graph = self.collapse_graph()
        pprint(self.collapsed_graph)
        self.simplified_graph = self.make_simplified_graph()
        self.find_simplified_cycles()
        self.find_cyclic_pairs()
        pprint(self.simplified_graph)
        print('simplified cycles')
        print(self.cycles)
        print('cyclic pairs')
        print(self.cyclic_pairs)
        print('cyclic nodes')
        print(self.cyclic_nodes)
        self.make_node_to_edge_dict()
        
        self.add_representations()
        self.make_smiles_components()
        print(self.components)

    def find_simplified_cycles(self):
        cycles = list(find_cycles.simple_cycles(self.simplified_graph))
        self.cycles = []
        for cycle in cycles:
            if len(cycle) > 2:
                new_cycle = True
                for cycle_2 in self.cycles:
                    if find_cycles.is_reverse_cycle(cycle, cycle_2):
                        new_cycle = False
                        break
                
                if new_cycle:
                    self.cycles.append(cycle)

    def find_cyclic_pairs(self):
        self.cyclic_pairs = set([])
        self.cyclic_nodes = set([])
        for cycle in self.cycles:
            for i, atom in enumerate(cycle):
                previous_atom = cycle[i - 1]
                self.cyclic_nodes.add(previous_atom)
                self.cyclic_nodes.add(atom)
                pair = tuple(sorted((previous_atom, atom), key = lambda x: x.nr))
                self.cyclic_pairs.add(pair)
                

    def is_branch_point(self, node):
        if self.get_neighbour_nr(node) >= 3:
            return True
        else:
            return False

    def get_neighbour_nr(self, node):
        neighbours = self.structure.graph[node]
        neighbour_nr = 0
        for neighbour in neighbours:
            if neighbour.type != 'H':
                neighbour_nr += 1

        return neighbour_nr

    def find_collapsed_graph_nodes(self):
        internal_graph_nodes = []
        end_nodes = []
        
        for node in self.structure.graph:
            if self.is_branch_point(node):
                internal_graph_nodes += [node]
            else:
                neighbour_nr = self.get_neighbour_nr(node)
                if neighbour_nr == 1:
                    end_nodes += [node]

        return internal_graph_nodes, end_nodes

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

    def make_simplified_graph(self):
        simplified_graph = {}
        for atom, next_atom_and_edges in self.collapsed_graph.items():
            simplified_graph[atom] = []
            for next_atom, next_edge, edge_atoms in next_atom_and_edges:
                simplified_graph[atom].append(next_atom)
        return simplified_graph

    def make_branch_dict(self):
        branch_dict = {}
        for node in self.simplified_graph:
            branch_dict[node] = 0

        return branch_dict

    def make_index_dict(self):
        index_dict = {}
        
        for node in self.simplified_graph:
            index_dict[node] = None

        return index_dict

    def make_node_to_edge_dict(self):
        self.node_to_edge_dict = {}
        self.node_to_atoms_dict = {}
        for atom, neighbours_and_edges in self.collapsed_graph.items():
            if atom.type != 'H':
                self.node_to_edge_dict[atom] = {}
                self.node_to_atoms_dict[atom] = {}
                for neighbour, edge, atoms in neighbours_and_edges:
                    if neighbour.type != 'H':
                        if not neighbour in self.node_to_edge_dict[atom]:
                            self.node_to_edge_dict[atom][neighbour] = []
                            self.node_to_atoms_dict[atom][neighbour] = []

                        self.node_to_edge_dict[atom][neighbour].append(edge)
                        self.node_to_atoms_dict[atom][neighbour].append(atoms)

    def adjust_node_indices(self, insert, break_point):
        index_jump = len(insert)

        for node, index in self.node_to_index.items():
            if index != None and index > break_point:
                self.node_to_index[node] += index_jump

    def add_insert(self, insert, break_point):
        half_1 = self.components[:break_point + 1]
        half_2 = self.components[break_point + 1:]

        self.adjust_node_indices(insert, break_point)

        self.components = half_1 + insert + half_2               

    def make_smiles_components(self):
        self.components = []
        self.atom_order = []
        working_graph = self.node_to_edge_dict

        index = 0
        used_nodes = set([])
        self.node_to_branch_nr = self.make_branch_dict()
        self.node_to_index = self.make_index_dict()

        while working_graph:
            last_node = None
            
            for node in working_graph:
                if node in used_nodes:
                    last_node = node
                    break

            if not last_node:
                last_node = list(working_graph.keys())[0]

            if not self.components:
                self.components = [self.representations[last_node]]
                self.node_to_index[last_node] = 0

            possible_next_step = True

            while possible_next_step:
                next_node = None
                
                self.node_to_branch_nr[last_node] += 1
                used_nodes.add(last_node)

                if last_node in working_graph:
                    next_node = list(working_graph[last_node].keys())[0]
                    sequence = working_graph[last_node][next_node]

                    break_point = self.node_to_index[last_node]

                    if self.node_to_branch_nr[last_node] == 1:
                        #This is the first branch we are making for this node
                        
                        if sequence:
                            insert = [sequence, self.representations[next_node]]
                        else:
                            insert = [self.representations[next_node]]

                    else:

                        if sequence:
                            insert = ['(', sequence, self.representations[next_node], ')']
                        else:
                            insert = ['(', self.representations[next_node], ')']

                    print('insert', insert)

                    self.add_insert(insert, break_point)
                    self.node_to_index[next_node] = self.node_to_index[last_node] + len(insert)
                        

                    used_nodes.add(next_node)
                    
                else:
                    possible_next_step = False
                    continue

                del working_graph[last_node][next_node]
                del working_graph[next_node][last_node]

                if not working_graph[last_node]:
                    del working_graph[last_node]
                    
                if not working_graph[next_node]:
                    del working_graph[next_node]

        print(self.components)

    def add_representations(self):
        self.representations = {}
        for atom in self.simplified_graph:
            if atom.type != 'H':
                if atom.aromatic:
                    self.representations[atom] = atom.type.lower()
                else:
                    self.representations[atom] = atom.type

        self.add_chiral_placeholders()
        self.add_hydrogen_representations()
        
        self.add_charge_representations()

    def add_chiral_placeholders(self):
        for atom in self.simplified_graph:
            if atom.chiral:
                self.representations[atom] = f'[{self.representations[atom]}X]'

    def add_charge_representations(self):
        for atom in self.simplified_graph:
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
        for atom, neighbours in self.simplified_graph.items():
            hydrogen_counts[atom] = 0
            for neighbour in neighbours:
                if neighbour.type == 'H':
                    hydrogen_counts[atom] += 1

        return hydrogen_counts

    def add_hydrogen_representations(self):
        hydrogen_counts = self.count_hydrogens()
        
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



    def collapse_graph(self):
        working_graph = self.structure.graph
        internal_nodes, end_nodes = self.find_collapsed_graph_nodes()
        
        try:
            first_node = end_nodes[0]
        except IndexError:
            try:
                first_node = internal_nodes[0]

            except IndexError:
                #this happens when the entire graph is cyclic
                first_node = list(self.structure.graph.keys())[0]
                internal_nodes = [first_node, first_node.neighbours[0]]

        collapsed_graph = {}

        for internal_node in internal_nodes:
            collapsed_graph[internal_node] = []
            
        for end_node in end_nodes:
            collapsed_graph[end_node] = []
            
        last_collapsed_node = first_node
        last_node = first_node
        previous_node = None
        edge_string = ''
        edge_atoms = []

        while last_node:
            if last_node in self.structure.graph:
                        
                next_node = self.structure.graph[last_node][0]
                
                node_type = next_node.type
                node_nr = next_node.nr

                node_type_last = last_collapsed_node.type
                node_nr_last = last_collapsed_node.nr
                
                edge_symbol = self.get_bond_symbol(last_node, next_node)

                if next_node in collapsed_graph:
                    edge_string += edge_symbol
                    collapsed_graph[last_collapsed_node].append((next_node, edge_string, edge_atoms))
                    collapsed_graph[next_node].append((last_collapsed_node, edge_string[::-1], edge_atoms))
                    
                    edge_string = ''
                    edge_atoms = []
                    self.remove_bonds_and_nodes(last_node, next_node)
                    
                    previous_collapsed_node = last_collapsed_node
                    last_collapsed_node = next_node
                    
                    
                else:
                    edge_string += edge_symbol
                    if not next_node.aromatic:
                        edge_string += node_type
                    else:
                        edge_string += node_type.lower()
                    edge_atoms.append(next_node)
                    
                    self.remove_bonds_and_nodes(last_node, next_node)
                
                last_node = next_node
            else:
                node_type = last_collapsed_node.type
                node_nr = last_collapsed_node.nr
                node_type_prev = previous_collapsed_node.type
                node_nr_prev = previous_collapsed_node.nr
                
                new_node = None

                for node in collapsed_graph:
                    if node in working_graph:
                        new_node = node
                        break
                    
                if not new_node:
                    last_node = None
                else:
                    last_node = new_node
                    last_collapsed_node = new_node

        return collapsed_graph

if __name__ == "__main__":
    smiles = 'c1ccc2c(c1)c(c[nH]2)C[C@@H](C(=O)O)N'
    smiles = 'C1=CC=C2C(=C1)C(=CN2)C[C@@H](C(=O)O)N'
 #   smiles = 'c12ccccc1cccc2'
 #   smiles = 'c1ccccc1'
    structure = pikachu.Smiles(smiles).smiles_to_structure()
    collapsed_structure = CollapsedGraph(structure)
 #   pprint(collapsed_structure.collapsed_graph)
 #   pprint(collapsed_structure.simplified_graph)
  #  pprint(collapsed_structure.representations)
 #   pprint(collapsed_structure.cycles)
#    pprint(collapsed_structure.simplified_graph)
#    pprint(collapsed_structure.representations)
#    pprint(collapsed_structure.node_to_edge_dict)
    
    

        
