#!/usr/bin/env python

import pikachu

def get_subgraph_size(atom, graph, seen_atoms):
    seen_atoms.add(atom)
    for neighbour in graph[atom]:
        if neighbour not in seen_atoms:
            get_subgraph_size(neighbour, graph, seen_atoms)

    return len(seen_atoms) - 1

if __name__ == "__main__":
    test_graph = {1: [2],
                  2: [1, 3],
                  3: [2, 4, 5],
                  4: [3, 6],
                  5: [3, 6],
                  6: [4, 5]}

    test_structure = pikachu.read_smiles('CCCCCCCCCC(=O)N[C@@H](CC1=CNC2=CC=CC=C21)C(=O)N[C@@H](CC(=O)N)C(=O)N[C@@H](CC(=O)O)C(=O)N[C@H]3[C@H](OC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@H](NC(=O)CNC(=O)[C@@H](NC(=O)[C@H](NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)CNC3=O)CCCN)CC(=O)O)C)CC(=O)O)CO)[C@H](C)CC(=O)O)CC(=O)C4=CC=CC=C4N)C')
    test_structure = pikachu.read_smiles('CCCCCCCCCC(=O)N[C@@H]1[C@H]([C@@H]([C@H](O[C@H]1OC2=C3C=C4C=C2OC5=C(C=C(C=C5)[C@H]([C@H]6C(=O)N[C@@H](C7=C(C(=CC(=C7)O)OC8[C@H]([C@H]([C@@H]([C@H](O8)CO)O)O)O)C9=C(C=CC(=C9)[C@H](C(=O)N6)NC(=O)[C@@H]4NC(=O)[C@@H]1C2=CC(=CC(=C2)OC2=C(C=CC(=C2)[C@H](C(=O)N[C@H](CC2=CC(=C(O3)C=C2)Cl)C(=O)N1)N)O)O)O)C(=O)O)O[C@H]1[C@@H]([C@H]([C@@H]([C@H](O1)CO)O)O)NC(=O)C)Cl)CO)O)O')

    for atom in test_structure.graph:
        if atom.type == 'C' and atom.nr == 15:
            atom_1 = atom
            break

    print(atom.neighbours)
    print(atom)
    print(atom.type)
    print(atom.nr)

    for atom in test_structure.graph:
        if atom.type == 'C' and atom.nr == 14:
            atom_2 = atom
            break

    for atom in test_structure.graph:
        if atom.type == 'C' and atom.nr == 23:
            atom_3 = atom
            break

    for atom in test_structure.graph:
        if atom.type == 'C' and atom.nr == 16:
            atom_4 = atom
            break

    print(len(test_structure.graph))

    print(get_subgraph_size(atom_2, test_structure.graph, set([atom_1])))
    print(get_subgraph_size(atom_3, test_structure.graph, set([atom_1])))
    print(get_subgraph_size(atom_4, test_structure.graph, set([atom_1])))
