#!/usr/bin/env python

import pikachu
import drawing
import math_functions

def get_subgraph_size(atom, graph, seen_atoms):
    seen_atoms.add(atom)
    for neighbour in graph[atom]:
        if neighbour not in seen_atoms:
            get_subgraph_size(neighbour, graph, seen_atoms)

    return len(seen_atoms) - 1

if __name__ == "__main__":
    test_structure = pikachu.read_smiles("C([C@@H](C(=O)O)N)O")
    drawing.drawer(test_structure)

