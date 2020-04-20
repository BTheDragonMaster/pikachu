# A dependency-free version of networkx's implementation of Johnson's cycle finding algorithm
# Original implementation: https://github.com/networkx/networkx/blob/master/networkx/algorithms/cycles.py#L109
# Original paper: Donald B Johnson. "Finding all the elementary circuits of a directed graph." SIAM Journal on Computing. 1975.

from collections import defaultdict

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
    G = {v: set(nbrs) for (v,nbrs) in G.items()} # make a copy of the graph
    sccs = strongly_connected_components(G)
    while sccs:
        scc = sccs.pop()
        startnode = scc.pop()
        path=[startnode]
        blocked = set()
        closed = set()
        blocked.add(startnode)
        B = defaultdict(set)
        stack = [ (startnode,list(G[startnode])) ]
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

if __name__ == "__main__":

    smiles = "CCCCCCCCCC(=O)N[C@@H](CC1=CNC2=CC=CC=C21)C(=O)N[C@H](CC(=O)N)C(=O)N[C@@H](CC(=O)O)C(=O)N[C@H]3[C@H](OC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@H](NC(=O)CNC(=O)[C@@H](NC(=O)[C@H](NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)CNC3=O)CCCN)CC(=O)O)C)CC(=O)O)CO)[C@H](C)CC(=O)O)CC(=O)C4=CC=CC=C4N)C"
 #   smiles = "CCCCCCCCCC(=O)N[C@@H]1[C@H]([C@@H]([C@H](O[C@H]1OC2=C3C=C4C=C2OC5=C(C=C(C=C5)[C@H]([C@H]6C(=O)N[C@@H](C7=C(C(=CC(=C7)O)OC8[C@H]([C@H]([C@@H]([C@H](O8)CO)O)O)O)C9=C(C=CC(=C9)[C@H](C(=O)N6)NC(=O)[C@@H]4NC(=O)[C@@H]1C2=CC(=CC(=C2)OC2=C(C=CC(=C2)[C@H](C(=O)N[C@H](CC2=CC(=C(O3)C=C2)Cl)C(=O)N1)N)O)O)O)C(=O)O)O[C@H]1[C@@H]([C@H]([C@@H]([C@H](O1)CO)O)O)NC(=O)C)Cl)CO)O)O"
    structure = Smiles(smiles).smiles_to_structure()
    peptide_bonds = find_bonds(structure, PEPTIDEBOND)
    coc_ester_bonds = find_bonds(structure, ESTERCOCBOND)
    bonds = coc_ester_bonds + peptide_bonds

    product = copy.deepcopy(structure)
 #   for bond in peptide_bonds:
 #       hydrolyse(bond, product)
        

    structures = product.split_disconnected_structures()
    for structure in structures:
        graph = structure.graph
        

    

        unique_cycles = set([])

        for cycle in (tuple(simple_cycles(graph))):
            if len(cycle) > 2:
                cycle_components = sorted(cycle, key = lambda x: x.nr)
                cycle_components = tuple(cycle_components)
                if len(cycle_components) < 10:
                    unique_cycles.add(cycle_components)

        for unique_cycle in unique_cycles:
            if len(unique_cycle) == 5:
                aromatic, heteroatom = check_five_ring(unique_cycle)
                if aromatic:
                    heteroatom.promote_lone_pair_to_p_orbital()
            print(unique_cycle, check_aromatic(unique_cycle))
            
            

        print(unique_cycles)
            
