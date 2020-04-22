
def calculate_energy_E(x_coords, y_coords, distances):
    pass


def calculate_length_lij(L, dij):
    return L * dij

def calculate_L(display_size, dijs):
    max_dij = max(dijs)
    return display_size / max_dij

def calculate_strength_kij(K, dij):
    return K / dij

def make_adjacency_list(structure, use_true_lengths = False):
    adjacency_list = {}
    for atom_1, neighbours in structure.graph.items():
        adjacency_list[atom_1.nr] = {}
        for atom_2 in structure.graph:
            if atom_1 == atom_2:
                adjacency_list[atom_1.nr][atom_2.nr] = 0
            
            elif atom_2 in neighbours:
                if true_lengths:
                    adjacency_list[atom_1.nr][atom_2.nr] = structure.graph.bond_lookup[atom_1][atom_2].bond_length
                else:
                    adjacency_list[atom_1.nr][atom_2.nr] = 1.0
            else:
                adjacency_list[atom_1.nr][atom_2.nr] = float('inf')

    return adjacency_list
    

def calculate_shortest_paths_floyd(structure, use_true_lengths = False):
    adjacency_list = make_adjacency_list(structure, use_true_lengths)

    for nr_1 in adjacency_list:
        for nr_2 in adjacency_list:
            for nr_3 in adjacency_list:
                adjacency_list[nr_2][nr_3] = min([adjacency_list[nr_2][nr_1] + adjacency_list[nr_1][nr_3],
                                                  adjacency_list[nr_2][nr_3]])

    return adjacency_list

if __name__ == "__main__":
    smiles = 'CCCCC'
    structure = Smiles(smiles).smiles_to_structure()
    pprint(calculate_shortest_path_floyd(structure))
    
                
                
    

    
    

    
                    
        
        
        
    
