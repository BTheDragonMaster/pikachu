#!/usr/bin/env python

from pprint import pprint
from copy import copy
import pikachu

class SSSR(pikachu.Structure):
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
            c_sssr = self.get_sssr(ring_candidates, d, cc_adjacency_matrix, pe, pe_prime, bond_counts, ring_counts, sssr_nr)

            for ring in c_sssr:
                rings.append(list(ring))

        return rings




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

            elif neighbour != parent[atom]:
                low[atom] = min(low[atom], disc[neighbour])

    def get_path_included_distance_matrices(self, adjacency_matrix):
        """
            Use Floyd-Warshall algorithm to compute the shortest paths between all vertice pairs in a graph

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
                if atom_1 == atom_2 or adjacency_matrix[atom_1][atom_2] == 1:
                    d[atom_1][atom_2] = adjacency_matrix[atom_1][atom_2]
                else:
                    d[atom_1][atom_2] = float('inf')

                if d[atom_1][atom_2] == 1:
                    pe[atom_1][atom_2] = [[[atom_1, atom_2]]]
                else:
                    pe[atom_1][atom_2] = []

                pe_prime[atom_1][atom_2] = []

        for atom_k in atoms:
            for atom_i in atoms:
                for atom_j in atoms:
                    previous_path_length = d[atom_i][atom_j]
                    new_path_length = d[atom_i][atom_k] + d[atom_k][atom_j]

                    if previous_path_length > new_path_length:
                        if previous_path_length == new_path_length + 1:
                            length_l = len(pe[atom_i][atom_j])
                            pe_prime[atom_i][atom_j] = [None] * length_l
                            for l in range(length_l - 1, -1, -1):
                                length_m = len(pe[atom_i][atom_j][l])
                                pe_prime[atom_i][atom_j][l] = [None] * length_m
                                for m in range(length_m - 1, -1, -1):
                                    length_n = len(pe[atom_i][atom_j][l][m])
                                    pe_prime[atom_i][atom_j][l][m] = [None] * length_n
                                    for n in range(length_n - 1, -1, -1):
                                        pe_prime[atom_i][atom_j][l][m][n] = [pe[atom_i][atom_j][l][m][0],
                                                                             pe[atom_i][atom_j][l][m][1]]
                        else:
                            pe_prime[atom_i][atom_j] = []

                        d[atom_i][atom_j] = new_path_length
                        pe[atom_i][atom_j] = [[]]


                        pe[atom_i][atom_j][0] = list(reversed(pe[atom_i][atom_k][0])) + list(reversed(pe[atom_k][atom_j][0]))

                    elif previous_path_length == new_path_length:
                        if len(pe[atom_i][atom_k]) and len(pe[atom_k][atom_j]):
                            tmp = list(reversed(pe[atom_i][atom_k][0])) + list(reversed(pe[atom_k][atom_j][0]))
                            pe[atom_i][atom_j].append(tmp)

                    elif previous_path_length == new_path_length - 1:
                        tmp = list(reversed(pe[atom_i][atom_k][0])) + list(reversed(pe[atom_k][atom_j][0]))
                        pe_prime[atom_i][atom_j].append(tmp)

        return d, pe, pe_prime

    def get_ring_candidates(self, d, pe, pe_prime):
        candidates = []

        vertices_in_cycle = 0

        for atom_1 in d:
            for atom_2 in d[atom_1]:
                if d[atom_1][atom_2] == 0 or (len(pe[atom_1][atom_2]) == 1 and pe_prime[atom_1][atom_2] == 0):
                    continue
                else:
                    if len(pe_prime[atom_1][atom_2]) != 0:
                        vertices_in_cycle = 2 * (d[atom_1][atom_2]) + 1
                    else:
                        vertices_in_cycle = 2 * (d[atom_1][atom_2])

                    if vertices_in_cycle != float('inf'):
                        candidates.append([vertices_in_cycle, pe[atom_1][atom_2], pe_prime[atom_1][atom_2]])

        candidates = sorted(candidates, key=lambda x: x[0])

        return candidates

    def get_sssr(self, ring_candidates, d, cc_adjacency_matrix, pe, pe_prime, bond_counts, ring_counts, sssr_nr):
        c_sssr = []
        all_bonds = []

        for candidate in ring_candidates:

            ring_size, paths, path_and_vertices = candidate
            if ring_size % 2 != 0:
                for path_and_vertex in path_and_vertices:
                    bonds = paths[0] + path_and_vertex
                    for i in range(len(bonds)):
                        if type(bonds[i][0]) == list:
                            bonds[i] = bonds[i][0]

                    atoms = self.bonds_to_atoms(bonds)
                    bond_count = self.get_bond_count(atoms, cc_adjacency_matrix)

                    if bond_count == len(atoms) and not self.path_sets_contain(c_sssr, atoms, bonds, all_bonds, bond_counts, ring_counts):
                        c_sssr.append(atoms)
                        print(c_sssr)
                        all_bonds += bonds

                    if len(c_sssr) > sssr_nr:
                        return c_sssr

            else:
                for i in range(len(paths) - 1):
                    bonds = paths[i] + paths[i + 1]
                    for i in range(len(bonds)):
                        if type(bonds[i][0]) == list:
                            bonds[i] = bonds[i][0]

                    atoms = self.bonds_to_atoms(bonds)
                    bond_count = self.get_bond_count(atoms, cc_adjacency_matrix)

                    if bond_count == len(atoms) and not self.path_sets_contain(c_sssr, atoms, bonds, all_bonds,
                                                                          bond_counts, ring_counts):
                        c_sssr.append(atoms)
                        all_bonds += bonds

                    if len(c_sssr) > sssr_nr:
                        return c_sssr

        return c_sssr


    def bonds_to_atoms(self, bonds):
        atoms = set()
        for atom_1, atom_2 in bonds:
            atoms.add(atom_1)
            atoms.add(atom_2)

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
        for i in range(len(c_sssr) - 1, -1, -1):
            if self.is_superset(atoms, c_sssr[i]):
                return True

            if len(c_sssr[i]) != len(atoms):
                continue

            if self.sets_equal(atoms, c_sssr[i]):
                return True

        count = 0
        all_contained = False

        for i in range(len(bonds) - 1, -1, -1):
            for j in range(len(all_bonds) - 1, -1, -1):
                if (bonds[i][0] == all_bonds[j][0] and bonds[i][1] == all_bonds[j][1]) or\
                        (bonds[i][1] == all_bonds[j][0] and bonds[i][0] == all_bonds[j][1]):
                    count += 1

                if count == len(bonds):
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


if __name__ == "__main__":
    daptomycin_smiles = pikachu.Smiles('CCCCCCCCCC(=O)N[C@@H](CC1=CNC2=CC=CC=C21)C(=O)N[C@@H](CC(=O)N)C(=O)N[C@@H](CC(=O)O)C(=O)N[C@H]3[C@H](OC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@H](NC(=O)CNC(=O)[C@@H](NC(=O)[C@H](NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)CNC3=O)CCCN)CC(=O)O)C)CC(=O)O)CO)[C@H](C)CC(=O)O)CC(=O)C4=CC=CC=C4N)C')
    daptomycin_structure = daptomycin_smiles.smiles_to_structure()
    teicoplanin_smiles = pikachu.Smiles('CCCCCCCCCC(=O)N[C@@H]1[C@H]([C@@H]([C@H](O[C@H]1OC2=C3C=C4C=C2OC5=C(C=C(C=C5)[C@H]([C@H]6C(=O)N[C@@H](C7=C(C(=CC(=C7)O)OC8[C@H]([C@H]([C@@H]([C@H](O8)CO)O)O)O)C9=C(C=CC(=C9)[C@H](C(=O)N6)NC(=O)[C@@H]4NC(=O)[C@@H]1C2=CC(=CC(=C2)OC2=C(C=CC(=C2)[C@H](C(=O)N[C@H](CC2=CC(=C(O3)C=C2)Cl)C(=O)N1)N)O)O)O)C(=O)O)O[C@H]1[C@@H]([C@H]([C@@H]([C@H](O1)CO)O)O)NC(=O)C)Cl)CO)O)O')
    teicoplanin_structure = teicoplanin_smiles.smiles_to_structure()
    graph = SSSR(teicoplanin_structure)
    rings = graph.get_rings()
    print(rings)



