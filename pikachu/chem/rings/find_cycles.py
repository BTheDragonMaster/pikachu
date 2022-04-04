# A dependency-free version of networkx's implementation of Johnson's cycle finding algorithm
# Original implementation: https://github.com/networkx/networkx/blob/master/networkx/algorithms/cycles.py#L109
# Original paper: Donald B Johnson. "Finding all the elementary circuits of a directed graph." SIAM Journal on Computing. 1975.

from collections import OrderedDict, defaultdict
import copy


def in_cycle(cycles, atom):
    for cycle in cycles:
        if atom in cycle:
            return True
    return False


def find_cycles(graph):
    all_cycles = simple_cycles(graph)
    cycles = []
    sorted_cycles = set()
    for cycle in all_cycles:
        if len(cycle) > 2:
            sorted_cycle = tuple(sorted(cycle, key=lambda x: x.nr))
            if sorted_cycle not in sorted_cycles:
                sorted_cycles.add(sorted_cycle)
                cycles.append(cycle)

    return cycles


def is_reverse_cycle(cycle_1, cycle_2):
    reversed_2 = list(reversed(cycle_2))

    starting_index = None

    for i, atom in enumerate(reversed_2):
        if atom == cycle_1[0]:
            starting_index = i

    if starting_index == None:
        return False

    else:
        new_reversed = reversed_2[i:] + reversed_2[:i]
        if new_reversed == cycle_1:
            return True
        else:
            return False


def simple_cycles(G):
    # Yield every elementary cycle in python graph G exactly once
    # Expects a dictionary mapping from vertices to iterables of vertices
    def _unblock(thisnode, blocked, B):
        stack = set([thisnode])
        while stack:
            node = stack.pop()
            if node in blocked:
                blocked.remove(node)
                stack.update(B[node])
                B[node].clear()
    G = {v: set(nbrs) for (v, nbrs) in G.items()} # make a copy of the graph
    sccs = strongly_connected_components(G)
    while sccs:
        scc = sccs.pop()
        startnode = scc.pop()
        path = [startnode]
        blocked = set()
        closed = set()
        blocked.add(startnode)
        B = defaultdict(set)
        stack = [ (startnode, list(G[startnode])) ]
        while stack:
            thisnode, nbrs = stack[-1]
            if nbrs:
                nextnode = nbrs.pop()
                if nextnode == startnode:
                    yield path[:]
                    closed.update(path)
                elif nextnode not in blocked:
                    path.append(nextnode)
                    stack.append( (nextnode,list(G[nextnode])) )
                    closed.discard(nextnode)
                    blocked.add(nextnode)
                    continue
            if not nbrs:
                if thisnode in closed:
                    _unblock(thisnode,blocked,B)
                else:
                    for nbr in G[thisnode]:
                        if thisnode not in B[nbr]:
                            B[nbr].add(thisnode)
                stack.pop()
                path.pop()
        remove_node(G, startnode)
        H = subgraph(G, set(scc))
        sccs.extend(strongly_connected_components(H))


def strongly_connected_components(graph):
    # Tarjan's algorithm for finding SCC's
    # Robert Tarjan. "Depth-first search and linear graph algorithms." SIAM journal on computing. 1972.
    # Code by Dries Verdegem, November 2012
    # Downloaded from http://www.logarithmic.net/pfh/blog/01208083168

    index_counter = [0]
    stack = []
    lowlink = {}
    index = {}
    result = []
    
    def _strong_connect(node):
        index[node] = index_counter[0]
        lowlink[node] = index_counter[0]
        index_counter[0] += 1
        stack.append(node)
    
        successors = graph[node]
        for successor in successors:
            if successor not in index:
                _strong_connect(successor)
                lowlink[node] = min(lowlink[node],lowlink[successor])
            elif successor in stack:
                lowlink[node] = min(lowlink[node],index[successor])

        if lowlink[node] == index[node]:
            connected_component = []

            while True:
                successor = stack.pop()
                connected_component.append(successor)
                if successor == node: break
            result.append(connected_component[:])
    
    for node in graph:
        if node not in index:
            _strong_connect(node)
    
    return result


def remove_node(G, target):
    # Completely remove a node from the graph
    # Expects values of G to be sets
    del G[target]
    for nbrs in G.values():
        nbrs.discard(target)


def subgraph(G, vertices):
    # Get the subgraph of G induced by set vertices
    # Expects values of G to be sets
    return {v: G[v] & vertices for v in vertices}


class Cycles:
    def __init__(self, structure):
        self.find_unique_cycles(structure)
        self.make_microcycle_graph()

    def find_x_membered(self, x):
        five_membered = []
        for cycle in self.unique_cycles:
            if len(cycle) == x:
                five_membered.append(cycle)

        return five_membered

    def find_minimal_cycles(self):
        nodes = set()
        cycles = sorted(self.unique_cycles, key=lambda x: len(x))

        minimal_cycles = []

        for cycle in cycles:
            already_covered = True
            for node in cycle:
                if node not in nodes:
                    already_covered = False
                nodes.add(node)

            if not already_covered:
                minimal_cycles.append(cycle)

        return minimal_cycles

    def find_sssr(self):
        sorted_cycles = sorted(self.all_cycles, key=lambda x: len(x))
        sssr = []
        atoms = set()
        for cycle in sorted_cycles:
            add_cycle = False
            for atom in cycle:

                if atom not in atoms:
                    add_cycle = True
                    break
            if add_cycle:
                sssr.append(cycle)
                
            for atom in cycle:
                atoms.add(atom)

        return sssr

    def find_unique_cycles(self, structure):
        all_cycles = simple_cycles(structure.graph)
        self.all_cycles = []

        unique_cycles = set()
        for cycle in all_cycles:

            if len(cycle) > 2:
                self.all_cycles.append(cycle)
                cycle_components = sorted(cycle, key=lambda x: x.nr)
                cycle_components = tuple(cycle_components)
                if len(cycle_components) < 10:
                    unique_cycles.add(cycle_components)
        
        self.unique_cycles = unique_cycles

    def make_microcycle_graph(self):
        self.graph = {}

        for cycle_1 in self.unique_cycles:
            if len(cycle_1) < 10:
                self.graph[cycle_1] = []
            for cycle_2 in self.unique_cycles:
                if cycle_1 != cycle_2:
                    if len(cycle_2) < 10:
                        if len(set(cycle_1).intersection(set(cycle_2))) > 1:
                            self.graph[cycle_1].append(cycle_2)

    def make_bond_nr_dict(self):
        self.bond_nr_dict = {}
        for cycle in self.graph:
            self.bond_nr_dict[cycle] = len(self.graph[cycle])

    def find_cyclic_systems(self):
        working_graph = copy.deepcopy(self)

        if working_graph.graph:
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

            for paths in paths_collection:
                if paths:
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

            new_cycles = []

            for new_graph in new_graphs:
                new_cycle = set([])
                for cycle in new_graph.keys():
                    for atom in cycle:
                        new_cycle.add(atom)
                new_cycles.append(tuple(new_cycle))

        else:
            new_cycles = []

        return new_cycles

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

    def remove_connectors(self):
        """Remove nodes that only have incoming edges from graph

        Input:
        working_graph: dict of {node: [node, ->], ->}, representing a graph

        """

        for node in self.graph:
            for next_node in self.graph[node]:
                if next_node not in self.graph:
                    self.graph[node].remove(next_node)

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
        
