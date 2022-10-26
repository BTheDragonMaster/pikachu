#!/usr/bin/env python

from pprint import pprint
from copy import copy
from pikachu.chem import structure


class SSSR(structure.Structure):
    def __init__(self, molecule):
        super().__init__(molecule.graph, molecule.bonds, molecule.bond_lookup)
        self.time = 0

    def get_rings(self):
        adjacency_matrix = self.get_component_adjacency_matrix()
        if not adjacency_matrix:
            return None

        connected_components = self.get_graph_components(adjacency_matrix)
        rings = []

        for component in connected_components:
            cc_adjacency_matrix = self.get_subgraph_adjacency_matrix(component)
            bond_counts = {}
            ring_counts = {}

            for atom_1 in cc_adjacency_matrix:
                bond_counts[atom_1] = 0
                ring_counts[atom_1] = 0
                for atom_2 in cc_adjacency_matrix[atom_1]:
                    bond_counts[atom_1] += cc_adjacency_matrix[atom_1][atom_2]

            edge_nr = 0
            for atom_1, atoms in cc_adjacency_matrix.items():
                for atom_2 in atoms:
                    edge_nr += cc_adjacency_matrix[atom_1][atom_2]

            assert edge_nr % 2 == 0
            edge_nr = edge_nr / 2

            sssr_nr = edge_nr - len(cc_adjacency_matrix) + 1

            all_three = True
            for atom, bond_count in bond_counts.items():
                if bond_count != 3:
                    all_three = False
                    break

            if all_three:
                sssr_nr = 2.0 + edge_nr - len(cc_adjacency_matrix)

            if sssr_nr == 1:
                rings.append(component)
                continue

            d, pe, pe_prime = self.get_path_included_distance_matrices(cc_adjacency_matrix)

            ring_candidates = self.get_ring_candidates(d, pe, pe_prime)
            c_sssr = self.get_sssr(ring_candidates, cc_adjacency_matrix, bond_counts, ring_counts, sssr_nr)

            for ring in c_sssr:
                original_ring = self.get_original_ring_order(list(ring))
                rings.append(original_ring)

        return rings

    def get_original_ring_order(self, ring):
        current_atom = ring[0]
        atoms = set(ring[1:])
        ordered_ring = [current_atom]

        while atoms:
            neighbours = self.graph[current_atom]
            for neighbour in neighbours:
                if neighbour in atoms:
                    atoms.remove(neighbour)
                    ordered_ring.append(neighbour)
                    current_atom = neighbour
                    break

        return ordered_ring

    def get_graph_components(self, adjacency_matrix):
        visited = {}
        components = []
        count = 0

        for atom in self.graph:
            visited[atom] = False

        for atom in self.graph:
            if not visited[atom]:
                component = []
                visited[atom] = True
                component.append(atom)
                count += 1
                self.dfs_components(atom, visited, adjacency_matrix, component)
                if len(component) > 1:
                    components.append(component)

        return components

    def dfs_components(self, atom, visited, adjacency_matrix, component):
        for neighbour in adjacency_matrix[atom]:
            is_adjacent = adjacency_matrix[atom][neighbour]

            if not is_adjacent or visited[neighbour] or atom == neighbour:
                continue

            visited[neighbour] = True
            component.append(neighbour)
            self.dfs_components(neighbour, visited, adjacency_matrix, component)

    def get_component_adjacency_matrix(self):
        adjacency_matrix = {}

        for atom_1 in self.graph:
            adjacency_matrix[atom_1] = {}
            for atom_2 in self.graph:
                adjacency_matrix[atom_1][atom_2] = 0

        for bond_name, bond in self.bonds.items():
            adjacency_matrix[bond.atom_1][bond.atom_2] = 1
            adjacency_matrix[bond.atom_2][bond.atom_1] = 1


        bridges = self.get_bridges()

        for bond in bridges:
            adjacency_matrix[bond.atom_1][bond.atom_2] = 0
            adjacency_matrix[bond.atom_2][bond.atom_1] = 0

        return adjacency_matrix

    def get_subgraph_adjacency_matrix(self, atoms):
        adjacency_matrix = {}

        for atom_1 in atoms:
            adjacency_matrix[atom_1] = {}
            for atom_2 in atoms:
                adjacency_matrix[atom_1][atom_2] = 0

        for atom_1 in atoms:
            for atom_2 in atoms:
                if atom_1 != atom_2:
                    if atom_1 in self.bond_lookup[atom_2]:
                        adjacency_matrix[atom_1][atom_2] = 1

        return adjacency_matrix

    def get_bridges(self):
        visited = {}
        disc = {}
        low = {}
        parent = {}
        bridges = []
        self.time = 0

        for atom in self.graph:
            visited[atom] = False
          #  disc[atom] = 0
            parent[atom] = None
          #  low[atom] = 0

        for atom in self.graph:
            if not visited[atom]:
                self.dfs_bridges(atom, visited, disc, low, parent, bridges)

        return bridges

    def dfs_bridges(self, atom, visited, disc, low, parent, bridges):
        visited[atom] = True
        disc[atom] = self.time
        low[atom] = self.time
        self.time += 1

        for neighbour in self.graph[atom]:

            if not visited[neighbour]:
                parent[neighbour] = atom
                self.dfs_bridges(neighbour, visited, disc, low, parent, bridges)

                low[atom] = min(low[atom], low[neighbour])

                if low[neighbour] > disc[atom]:
                    bridges.append(self.bond_lookup[atom][neighbour])

            elif not parent[atom]:
                low[atom] = min(low[atom], disc[neighbour])
            elif neighbour != parent[atom]:
                low[atom] = min(low[atom], disc[neighbour])

    def get_path_included_distance_matrices(self, adjacency_matrix):
        """
            Use Floyd-Warshall algorithm to compute the shortest paths between all vertex pairs in a graph

        """
        atoms = list(adjacency_matrix.keys())
        length = len(atoms)
        d = {}
        pe = {}
        pe_prime = {}

        i = length

        for atom_1 in atoms:

            d[atom_1] = {}
            pe[atom_1] = {}
            pe_prime[atom_1] = {}

            for atom_2 in atoms:
                # Initialise shortest paths to infinity, 0 (same atom) or 1 (neighbours)
                if atom_1 == atom_2 or adjacency_matrix[atom_1][atom_2] == 1:
                    d[atom_1][atom_2] = adjacency_matrix[atom_1][atom_2]
                else:
                    d[atom_1][atom_2] = float('inf')

                # For neighbours: set the pe

                if d[atom_1][atom_2] == 1:
                    pe[atom_1][atom_2] = [[self.bond_lookup[atom_1][atom_2]]]
                else:
                    pe[atom_1][atom_2] = []

                pe_prime[atom_1][atom_2] = set()

        # For each atom triplet, check if the sum of the distances between two of the three
        # atom pairs is smaller than the current distance recorded between the remaining atom
        # pair

        for atom_k in atoms:
            for atom_i in atoms:
                for atom_j in atoms:

                    previous_path_length = d[atom_i][atom_j]
                    new_path_length = d[atom_i][atom_k] + d[atom_k][atom_j]

                    # The path from atom_i to atom_j through atom_k is shorter than the current path
                    # from atom_i to atom_j

                    if previous_path_length > new_path_length:
                        # pe_prime[atom_i][atom_j] = []
                        # If it is exactly one shorter, we need to store the old path in the pe_prime matrix
                        # if previous_path_length == new_path_length + 1:
                        #     for path in pe[atom_i][atom_j]:
                        #         pe_prime[atom_i][atom_j].append(path[:])

                        d[atom_i][atom_j] = new_path_length
                        pe[atom_i][atom_j] = []

                        for path_1 in pe[atom_i][atom_k]:
                            for path_2 in pe[atom_k][atom_j]:
                                new_path = path_1[:] + path_2[:]
                                pe[atom_i][atom_j].append(new_path)

                    elif previous_path_length == new_path_length:
                        if len(pe[atom_i][atom_k]) and len(pe[atom_k][atom_j]):
                            for path_1 in pe[atom_i][atom_k]:
                                for path_2 in pe[atom_k][atom_j]:
                                    new_path = path_1[:] + path_2[:]
                                    pe[atom_i][atom_j].append(new_path)

                    # elif previous_path_length == new_path_length - 1:
                    #     for path_1 in pe[atom_i][atom_k]:
                    #         for path_2 in pe[atom_k][atom_j]:
                    #             new_path = path_1[:] + path_2[:]
                    #             pe_prime[atom_i][atom_j].append(new_path)

        for atom_k in atoms:
            for atom_i in atoms:
                for atom_j in atoms:
                    shortest_path_length = d[atom_i][atom_j]
                    new_path_length = d[atom_i][atom_k] + d[atom_k][atom_j]
                    if new_path_length - 1 == shortest_path_length:
                        for path_1 in pe[atom_i][atom_k]:
                            for path_2 in pe[atom_k][atom_j]:
                                new_path = tuple(path_1[:] + path_2[:])
                                pe_prime[atom_i][atom_j].add(new_path)

        return d, pe, pe_prime

    def get_ring_candidates(self, d, pe, pe_prime):
        candidates = []

        vertices_in_cycle = 0

        for atom_1 in d:
            for atom_2 in d[atom_1]:
                # If the atom distance is 0, or if there's only one shortest path and there is no 
                # shortest path one longer than the shortest path
                if d[atom_1][atom_2] == 0 or (len(pe[atom_1][atom_2]) == 1 and len(pe_prime[atom_1][atom_2]) == 0):
                    continue
                else:
                    if len(pe[atom_1][atom_2]) > 1:
                        vertices_in_cycle = 2 * (d[atom_1][atom_2])
                    elif len(pe_prime[atom_1][atom_2]) != 0:
                        vertices_in_cycle = 2 * (d[atom_1][atom_2]) + 1
                    # else:
                    #     vertices_in_cycle = 2 * (d[atom_1][atom_2])

                    if vertices_in_cycle != float('inf'):
                        candidates.append([vertices_in_cycle, pe[atom_1][atom_2], pe_prime[atom_1][atom_2]])

        candidates = sorted(candidates, key=lambda x: x[0])

        return candidates

    def get_sssr(self, ring_candidates, cc_adjacency_matrix, bond_counts, ring_counts, sssr_nr):
        c_sssr = []
        all_bonds = set()

        for candidate in ring_candidates:

            ring_size, paths, extended_paths = candidate
            if ring_size % 2 != 0:
                for extended_path in extended_paths:
                    bonds = paths[0] + list(extended_path)

                    atoms = self.bonds_to_atoms(bonds)
                    bond_count = self.get_bond_count(atoms, cc_adjacency_matrix)

                    if bond_count == len(atoms) and not self.path_sets_contain(c_sssr, atoms, bonds, all_bonds,
                                                                               bond_counts, ring_counts):
                        c_sssr.append(atoms)
                        for bond in bonds:
                            all_bonds.add(bond)

                    if len(c_sssr) > sssr_nr:
                        return c_sssr

            else:
                for i in range(len(paths) - 1):
                    bonds = paths[i] + paths[i + 1]

                    atoms = self.bonds_to_atoms(bonds)
                    bond_count = self.get_bond_count(atoms, cc_adjacency_matrix)

                    if bond_count == len(atoms) and not self.path_sets_contain(c_sssr, atoms, bonds, all_bonds,
                                                                          bond_counts, ring_counts):
                        c_sssr.append(atoms)
                        for bond in bonds:
                            all_bonds.add(bond)

                    if len(c_sssr) > sssr_nr:
                        return c_sssr

        return c_sssr

    def bonds_to_atoms(self, bonds):
        atoms = set()
        for bond in bonds:
            atoms.add(bond.atom_1)
            atoms.add(bond.atom_2)

        return atoms

    def get_bond_count(self, atoms, adjacency_matrix):
        count = 0

        for atom_1 in atoms:
            for atom_2 in atoms:
                if not atom_1 == atom_2:
                    count += adjacency_matrix[atom_1][atom_2]

        return count / 2

    def is_superset(self, set_1, set_2):
        for element in set_2:
            if element not in set_1:
                return False
        return True

    def sets_equal(self, set_1, set_2):
        if len(set_1) != len(set_2):
            return False

        for atom in set_1:
            if not atom in set_2:
                return False

        return True

    def path_sets_contain(self, c_sssr, atoms, bonds, all_bonds, bond_counts, ring_counts):
        for candidate_ring in c_sssr:
            if self.is_superset(atoms, candidate_ring):
                return True

            if len(candidate_ring) != len(atoms):
                continue

            if self.sets_equal(atoms, candidate_ring):
                return True

        all_contained = False

        if self.is_superset(all_bonds, bonds):
            all_contained = True

        #special_case - see smiles drawer code
        special_case = False

        if all_contained:
            for atom in atoms:
                if ring_counts[atom] < bond_counts[atom]:
                    special_case = True
                    break

        if all_contained and not special_case:
            return True

        for atom in atoms:
            ring_counts[atom] += 1

        return False
