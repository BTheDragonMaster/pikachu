#!/usr/bin/env python

class CollapsedGraph:
    def __init__(self, structure):
        self.structure = copy.deepcopy(structure)

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
            if is_branch_point(graph, node):
                internal_graph_nodes += [node]
            else:
                neighbour_nr = get_neighbour_nr(graph, node)
                if neighbour_nr == 1:
                    end_nodes += [node]

        return internal_graph_nodes, end_nodes

    def get_bond_symbol(self, last_node, next_node):
        bond_type = self.structure.graph[last_node][next_node].type

    def collapse_graph(self):
        working_graph = self.structure.graph
        internal_nodes, end_nodes = self.find_collapsed_graph_nodes()
        
        try:
            first_node = end_nodes[0]
        except IndexError:
            try:
                first_node = internal_nodes[0]

            except IndexError:
                print("What kind of wonky group do wo have this time..")
                pprint.pprint(self.structure.graph)

        collapsed_graph = {}

        for internal_node in internal_nodes:
            collapsed_graph[internal_node] = []
            
        for end_node in end_nodes:
            collapsed_graph[end_node] = []
            
        last_collapsed_node = first_node
        last_node = first_node
        previous_node = None
        edge_string = ''

        while last_node:
            try:
                next_node = working_graph[last_node][0]
                
                node_type = next_node[0]
                node_nr = next_node[1]

                node_type_last = last_collapsed_node[0]
                node_nr_last = last_collapsed_node[1]
                
                edge_symbol = find_edge_valence_symbol(working_graph,
                                                       last_node, next_node)

                if next_node in collapsed_graph:
                    edge_string += edge_symbol
                    collapsed_graph[last_collapsed_node] += [(node_type, node_nr,
                                                              edge_string)]
                    if not last_collapsed_node == next_node:
                        collapsed_graph[next_node] += [(node_type_last,
                                                        node_nr_last,
                                                        edge_string[::-1])]
                    
                    edge_string = ''
                    remove_bonds_and_nodes(working_graph, last_node, next_node)
                    
                    previous_collapsed_node = last_collapsed_node
                    last_collapsed_node = next_node
                    
                    
                else:
                    edge_string += edge_symbol
                    edge_string += node_type
                    remove_bonds_and_nodes(working_graph, last_node, next_node)
                    
                
                
                last_node = next_node
            except KeyError:
                node_type = last_collapsed_node[0]
                node_nr = last_collapsed_node[1]
                node_type_prev = previous_collapsed_node[0]
                node_nr_prev = previous_collapsed_node[1]
                
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

        
