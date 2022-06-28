import logging


class Node:

    index = 0

    def __init__(self):
        self.atom = None
        self.neighbors = []
        self.is_visited = False
        self.parent = None
        self.mate = None
        self.index = Node.index
        Node.index += 1

    def __repr__(self):
        return str(self.atom)

    def set_atom(self, atom):
        self.atom = atom


class SuperNode(Node):

    def __init__(self):
        Node.__init__(self)
        self.subnodes = []
        self.original_edges = []

    def circle(self, node):
        for i, v in enumerate(self.subnodes):
            if v == node:
                break
        assert i < len(self.subnodes)

        if i > 0 and self.subnodes[i].mate == self.subnodes[i - 1] or i == 0 and self.subnodes[i].mate == self.subnodes[
            -1]:
            return self.subnodes[i::-1] + self.subnodes[:i:-1]
        else:
            return self.subnodes[i::] + self.subnodes[:i]


class Path:

    def __init__(self):
        self.nodes = []

    def head(self):
        return self.nodes[0]

    def tail(self):
        return self.nodes[-1]

    def replace_head(self):
        assert isinstance(self.nodes[0], SuperNode)
        snode = self.nodes.pop(0)
        for node in snode.subnodes:
            if self.nodes[0] in node.neighbors:
                if node.mate is None:
                    self.nodes.insert(0, node)
                else:
                    for v in snode.circle(node):
                        self.nodes.insert(0, v)
                return
        logging.error("cannot replace head node.")

    def replace_tail(self):
        assert isinstance(self.nodes[-1], SuperNode)
        snode = self.nodes.pop()
        for node in snode.subnodes:
            if self.nodes[-1] in node.neighbors:
                if node.mate is None:
                    self.nodes.append(node)
                else:
                    for v in snode.circle(node):
                        self.nodes.append(v)
                return
        logging.error("cannot replace tail node.")

    def __repr__(self):
        return str(self.nodes)


class Match:

    def __init__(self, nodes):
        self.nodes = nodes
        self.freenodes = []
        for node in nodes:
            self.freenodes.append(node)
        self.supernodes = []

    @staticmethod
    def from_structure(structure):
        pi_subgraph = structure.find_pi_subgraph()
        node_nr = len(pi_subgraph.keys())

        atom_to_node = {}
        nodes = [Node() for i in range(node_nr)]

        for i, atom in enumerate(pi_subgraph):
            atom_to_node[atom] = i

        for atom, neighbours in pi_subgraph.items():
            index_1 = atom_to_node[atom]
            nodes[index_1].set_atom(atom)

            for neighbour in neighbours:
                index_2 = atom_to_node[neighbour]
                nodes[index_1].neighbors.append(nodes[index_2])

        return Match(nodes)

    @staticmethod
    def from_edges(N, edges):
        nodes = [Node() for i in range(N)]
        for i, j in edges:
            nodes[i].neighbors.append(nodes[j])
        return Match(nodes)

    def clear_nodes(self):
        for node in self.nodes:
            node.is_visited = False
            node.parent = None

    def find_augmenting_path(self, root):
        self.clear_nodes()
        queue = [root]
        while len(queue) > 0:
            cur_node = queue.pop(0)
            cur_node.is_visited = True
            for node in cur_node.neighbors:
                if node == cur_node.parent:
                    continue

                elif node.is_visited:
                    cycle = self.find_cycles(node, cur_node)
                    if len(cycle) % 2 == 1:
                        logging.debug('blossom: {}'.format(cycle))
                        snode = self.shrink_blossom(cycle)
                        self.supernodes.append(snode)
                        for v in cycle:
                            if v in queue:
                                queue.remove(v)
                            if v.is_visited:
                                snode.is_visited = True
                                snode.parent = v.parent
                        queue.append(snode)
                        break

                elif node.mate is None:
                    node.parent = cur_node
                    return self.construct_augmenting_path(node)

                else:
                    node.is_visited = True
                    node.mate.is_visited = True
                    node.parent = cur_node
                    node.mate.parent = node
                    queue.append(node.mate)
        raise Exception('cannot find an augmenting path')

    def unmatched_nodes(self):
        self.maximum_matching()

        count = 0
        for node in self.nodes:
            if node.mate is not None:
                count += 1

        return len(self.nodes) - count

    def maximum_matching(self):
        while len(self.freenodes) > 0:
            logging.debug('freenodes: {}'.format(self.freenodes))

            for node in self.freenodes:
                try:
                    path = self.find_augmenting_path(node)
                    logging.debug('augmenting path: {}'.format(path.nodes))
                    self.invert_path(path)
                    self.freenodes.remove(path.nodes[0])
                    self.freenodes.remove(path.nodes[-1])
                    break
                except Exception as e:
                    logging.info(e)
            else:
                logging.info('Tried all free nodes, no more augmenting path.')

                break

    def invert_path(self, path):
        assert len(path.nodes) % 2 == 0
        for i in range(0, len(path.nodes), 2):
            path.nodes[i].mate = path.nodes[i + 1]
            path.nodes[i + 1].mate = path.nodes[i]

    def construct_augmenting_path(self, node):
        path = Path()
        path.nodes.append(node)
        node = node.parent
        path.nodes.append(node)
        while node.mate is not None:
            node = node.parent
            path.nodes.append(node)

        while len(self.supernodes) > 0:
            snode = self.supernodes.pop()
            self.expand_supernode(snode)
            if snode == path.head():
                path.replace_head()
            elif snode == path.tail():
                path.replace_tail()

        while path.nodes[0].mate is not None:
            path.nodes.insert(path.nodes[0].parent, 0)

        while path.nodes[-1].mate is not None:
            path.nodes.append(path.nodes[-1].parent)

        return path

    def find_cycles(self, node1, node2):
        def find_ancestors(node):
            ancestors = [node]
            while node.parent is not None:
                node = node.parent
                ancestors.append(node)
            return ancestors

        ancestors1 = find_ancestors(node1)
        ancestors2 = find_ancestors(node2)
        i = len(ancestors1) - 1
        j = len(ancestors2) - 1
        while ancestors1[i] == ancestors2[j]:
            i -= 1
            j -= 1

        cycle = ancestors1[:i + 1] + ancestors2[j + 1::-1]
        return cycle

    def shrink_blossom(self, blossom):
        snode = SuperNode()
        for node in blossom:
            snode.subnodes.append(node)
            for adj_node in node.neighbors:
                if adj_node not in blossom:
                    snode.original_edges.append((node, adj_node))

        for node1, node2 in snode.original_edges:
            node1.neighbors.remove(node2)
            node2.neighbors.remove(node1)
            node2.neighbors.append(snode)
            snode.neighbors.append(node2)

        return snode

    def expand_supernode(self, snode):
        assert isinstance(snode, SuperNode)
        for node1, node2 in snode.original_edges:
            node1.neighbors.append(node2)
            node2.neighbors.append(node1)
            node2.neighbors.remove(snode)
            snode.neighbors.remove(node2)
