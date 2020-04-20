#import matplotlib
#matplotlib.use('TkAgg')
#from matplotlib import pyplot as plt
from pprint import pprint
import copy
from collections import OrderedDict, defaultdict
from math import sqrt


def squared_distance(vector_1, vector_2):

    assert len(vector_1) == len(vector_2)

    total_sum = 0

    for i, coord_1 in enumerate(vector_1):
        coord_2 = vector_2[i]
        diff = coord_1 - coord_2
        total_sum += (diff**2)
        
    return total_sum

def force_to_vector(force, atom_1, atom_2):

    eu_diff = euclidean(atom_2.get_coords(), atom_1.get_coords())
    scale = force / eu_diff
    
    x_diff = atom_2.x - atom_1.x
    y_diff = atom_2.y - atom_1.y
    z_diff = atom_2.z - atom_1.z

    x_diff_scaled = x_diff * scale
    y_diff_scaled = y_diff * scale
    z_diff_scaled = z_diff * scale

    vector = [x_diff_scaled, y_diff_scaled, z_diff_scaled]
    return vector
    

def euclidean(vector_1, vector_2):
    assert len(vector_1) == len(vector_2)

    total_sum = 0

    for i, coord_1 in enumerate(vector_1):
        coord_2 = vector_2[i]
        diff = coord_1 - coord_2
        total_sum += (diff**2)

    euclidean = sqrt(total_sum)
    return euclidean

def compare_matches(match_1, match_2):
    matching = True
    for key in match_1:
        if not key in match_2:
            matching = False
            break

        if match_1[key] != match_2[key]:
            matching = False
            break

    return matching

def compare_all_matches(matches):
    matching_pairs = set([])
    
    for i, match_1 in enumerate(matches):
        for j, match_2 in enumerate(matches):
            if i != j:
                if compare_matches(match_1, match_2):
                    matching_pairs.add(tuple(sorted([i, j])))

    matches_to_remove = set([])

    for matching_pair in matching_pairs:
        matches_to_remove.add(matching_pair[1])

    matches_to_remove = sorted(list(matches_to_remove), reverse = True)
    
    for match_to_remove in matches_to_remove:
        del matches[match_to_remove]

def check_aromatic(atom_set):
    aromatic = True
    for atom in atom_set:
        if atom.hybridisation == 'sp2':
            pass
        else:
            aromatic = False
            break

    if aromatic:
        pi_electron_nr = 0
        for atom in atom_set:
            for orbital_name in atom.valence_shell.orbitals:
                orbital = atom.valence_shell.orbitals[orbital_name]
                if orbital.orbital_type == 'p':
                    for electron in orbital.electrons:
                        if electron.atom == atom:
                            pi_electron_nr += 1

        if not pi_electron_nr % 4 == 2:
            aromatic = False

    return aromatic

def check_five_ring(atom_set):
    assert len(atom_set) == 5

    sp2_hybridised = []
    sp3_hybridised_lone_pair = []

    aromatic = False
    heteroatom = None

    for atom in atom_set:
        if atom.hybridisation == 'sp2':
            sp2_hybridised.append(atom)
        elif atom.hybridisation == 'sp3':
            if atom.calc_electron_pair_nr() > 0:
                sp3_hybridised_lone_pair.append(atom)

    if len(sp2_hybridised) == 4 and len(sp3_hybridised_lone_pair) == 1:
        
        pi_electron_nr = 0
        for atom in sp2_hybridised:
            for orbital_name in atom.valence_shell.orbitals:
                orbital = atom.valence_shell.orbitals[orbital_name]
                if orbital.orbital_type == 'p':
                    for electron in orbital.electrons:
                        if electron.atom == atom:
                            pi_electron_nr += 1
            if pi_electron_nr % 4 == 0:
                aromatic = True
                heteroatom = sp3_hybridised_lone_pair[0]

    return aromatic, heteroatom
        
    

def get_chiral_permutations(order):
    permutations = [tuple(order)]
    permutations.append((order[0], order[3], order[1], order[2]))
    permutations.append((order[0], order[2], order[3], order[1]))
    permutations.append((order[1], order[0], order[3], order[2]))
    permutations.append((order[1], order[2], order[0], order[3]))
    permutations.append((order[1], order[3], order[2], order[0]))
    permutations.append((order[2], order[0], order[1], order[3]))
    permutations.append((order[2], order[3], order[0], order[1]))
    permutations.append((order[2], order[1], order[3], order[0]))
    permutations.append((order[3], order[0], order[2], order[1]))
    permutations.append((order[3], order[1], order[0], order[2]))
    permutations.append((order[3], order[2], order[1], order[0]))

    return permutations


def check_same_chirality(neighbours_1, chirality_1, neighbours_2, chirality_2, match):
    equivalent_atom_list = []
    for atom in neighbours_1:
        if atom.type == 'H':
            for atom_2 in neighbours_2:
                if atom_2.type == 'H':
                    equivalent_atom_list.append(atom_2)
                    break
        else:
            equivalent_atom_list.append(match[atom])
    chiral_permutations_1 = get_chiral_permutations(equivalent_atom_list)

    if chirality_1 == chirality_2:
        if tuple(neighbours_2) in chiral_permutations_1:
            return True
        else:
            return False
    else:
        if tuple(neighbours_2) in chiral_permutations_1:
            return False
        else:
            return True
    


     



class Structure():

    """
    structure: dict of {atom: [atom, ->], ->}, with each atom a tuple of
        (str, int), str representing atom type, and int representing atom
        number. Example: ('O', 1)
    """
   
    def __init__(self, graph = None, bonds = None, bond_lookup = None):
        if graph:
            self.graph = graph
        else:
            self.graph = {}
        if bonds:
            self.bonds = bonds
        else:
            self.bonds = {}

        if bond_lookup:
            self.bond_lookup = bond_lookup
        else:
            self.bond_lookup = {}

    

    def find_next_atom_nr(self):
        atom_nrs = []

        for atom in self.graph:
            atom_nrs.append(atom.nr)

        atom_nrs.sort()

        next_atom_nr = None

        for i, atom_nr in enumerate(atom_nrs):
            if i != atom_nr:
                next_atom_nr = i
                break
        
        if next_atom_nr == None:
            next_atom_nr = atom_nrs[-1] + 1

        return next_atom_nr
            

    def get_steric_atoms(self):
        steric_atoms = []
        for atom in self.graph:
            if atom.chiral:
                steric_atoms.append(atom)

        return steric_atoms

        

    def find_next_bond_nr(self):
        bond_nrs = list(self.bonds.keys())
        bond_nrs.sort()

        next_bond_nr = None

        for i, bond_nr in enumerate(bond_nrs):
            if i != bond_nr:
                next_bond_nr = i
                break
        
        if next_bond_nr == None:
            next_bond_nr = bond_nrs[-1] + 1

        return next_bond_nr

    def make_bond_lookup(self):
        self.bond_lookup = {}
        for bond in self.bonds:
            atom_1 = self.bonds[bond].atom_1
            atom_2 = self.bonds[bond].atom_2
            if not atom_1 in self.bond_lookup:
                self.bond_lookup[atom_1] = {}
            if not atom_2 in self.bond_lookup:
                self.bond_lookup[atom_2] = {}
            self.bond_lookup[atom_1][atom_2] = self.bonds[bond]
            self.bond_lookup[atom_2][atom_1] = self.bonds[bond]

    def make_atom_dict(self):
        self.atom_dict = {}
        for atom in self.graph:
            self.atom_dict[atom.nr] = atom

    def make_atom(self, atom_type, neighbours, chiral = None, charge = 0):

        
        next_atom_nr = self.find_next_atom_nr()
        atom = Atom(atom_type, next_atom_nr, chiral, charge)
        atom.add_shell_layout()

        
        for i, neighbour in enumerate(neighbours):
            next_bond_nr = self.find_next_bond_nr()
            self.make_bond(atom, neighbour, next_bond_nr)

        return atom

    def find_cycles(self):
        self.cycles = find_cycles.Cycles(self)

    def fix_electrons_five_rings(self):
        five_rings = self.cycles.find_five_membered()
        for five_ring in five_rings:
            aromatic, heteroatom = check_five_ring(five_ring)
            if aromatic:
                heteroatom.promote_lone_pair_to_p_orbital()
                
    def find_aromatic_systems(self):
        ring_systems = self.cycles.find_cyclic_systems()
        for ring_system in ring_systems:
            print(ring_system, check_aromatic(ring_system))


    def refine_structure(self):
        self.make_atoms_upper_case()
        
        self.add_shells()
        self.add_hydrogens()
        self.add_shells()
        
        self.sort_by_nr()

        self.refine_p_bonds()
        self.hybridise_atoms()
        self.refine_s_bonds()
        self.drop_electrons()
        self.set_atom_neighbours()
        self.set_connectivities()
        self.make_bond_lookup()
        self.make_atom_dict()
        self.find_cycles()
        self.fix_electrons_five_rings()

    def break_bond_nr(self, bond):
        if type(bond) == int:
            bond_nr = bond
            bond = self.bonds[bond]
            

        else:
            bond_nr = bond.nr
            
        bond.break_bond()
        
        del self.bonds[bond_nr]


    def remove_atom(self, atom_to_remove):

        for bond in atom_to_remove.bonds:
            self.break_bond_nr(bond.nr)
            
        for atom in self.graph:
            if atom_to_remove in self.graph[atom]:
                self.graph[atom].remove(atom_to_remove)

        del self.graph[atom_to_remove]

    def break_bond_atoms(self, atom_1, atom_2):
        if atom_1.type == int:
            atom_1 = self.atom_dict[atom_1]
        if atom_2.type == int:
            atom_2 = self.atom_dict[atom_2]

        bond = self.bond_lookup[atom_1][atom_2]

        bond_nr = bond.nr

        bond.break_bond()

        del self.bonds[bond_nr]
        
        self.graph[bond.atom_1].remove(bond.atom_2)
        self.graph[bond.atom_2].remove(bond.atom_1)

    def set_connectivities(self):
        for atom in self.graph:
            if atom.type != 'H':
                atom.set_connectivity()

    def set_atom_neighbours(self):
        for atom in self.graph:
            atom.set_neighbours(self)

    def get_connectivities(self):
        connectivities = {}
        for atom in self.graph:
            if atom.type != 'H':
                connectivity = atom.connectivity
                if not connectivity in connectivities:
                    connectivities[connectivity] = []
                
                    connectivities[connectivity].append(atom)

        
        return connectivities

    def get_chiral_double_bonds(self):
        chiral_bonds = []
        for bond_nr in self.bonds:
            bond = self.bonds[bond_nr]
            if bond.chiral:
                chiral_bonds.append(bond)

        return chiral_bonds
                

    def check_chiral_double_bonds(self, child, match):
        chirality_matches = True
        chiral_bonds = child.get_chiral_double_bonds()

        for chiral_bond in chiral_bonds:
            neighbour_1, neighbour_2 = chiral_bond.neighbours

            parent_neighbour_1 = match[neighbour_1]
            parent_neighbour_2 = match[neighbour_2]
            parent_bond = self.bond_lookup[parent_neighbour_1][parent_neighbour_2]

            if not parent_bond.chiral:
                chirality_matches = False
                break
            else:
                matching_chirality = chiral_bond.check_same_chirality(parent_bond, match)
                if not matching_chirality:
                    chirality_matches = False
                    break

        return chirality_matches

            
        

    def check_chiral_centres(self, child, match):
        chirality_matches = True
        chiral_centres = child.get_steric_atoms()
        
        for chiral_centre in chiral_centres:
            parent_atom = match[chiral_centre]
            if parent_atom.chiral:
                child_neighbours = chiral_centre.neighbours
                parent_neighbours = parent_atom.neighbours
                
                chirality_matches = check_same_chirality(child_neighbours, chiral_centre.chiral,
                                                         parent_neighbours, parent_atom.chiral,
                                                         match)
                if not chirality_matches:
                    break
            else:
                chirality_matches = False
                break
            
        return chirality_matches

    def substructure_matching(self, child, check_chiral_centres = True, check_chiral_double_bonds = True):
        matches = []
        if self.is_substructure_atom_composition(child):
            if self.is_substructure_atom_connectivity(child):
                matches = self.is_substructure_connected_atoms(child)

        if check_chiral_centres:
            final_matches = []

            for match in matches:
                if self.check_chiral_centres(child, match):
                    final_matches.append(match)
            matches = final_matches

        if check_chiral_double_bonds:
            final_matches = []
            for match in matches:
                if self.check_chiral_double_bonds(child, match):
                    final_matches.append(match)
            matches = final_matches

        return matches

    def is_substructure_atom_composition(self, child):
        
        atom_counts_self = self.get_atom_counts()
        atom_counts_child = child.get_atom_counts()

        can_be_substructure = True

        for atom_type in atom_counts_child:
            try:
                atom_nr_self = atom_counts_self[atom_type]
                atom_nr_child = atom_counts_child[atom_type]
                if atom_nr_child > atom_nr_self:
                    can_be_substructure = False
                    break
            except KeyError:
                can_be_substructure = False
                break

        return can_be_substructure

    def get_connectivity_counts(self):
        connectivities = {}
        for atom in self.graph:
            if atom.type != 'H':
                if not atom.type in connectivities:
                    connectivities[atom.type] = {}
                connectivity = atom.connectivity
                if not connectivity in connectivities[atom.type]:
                    connectivities[atom.type][connectivity] = 0
                connectivities[atom.type][connectivity] += 1

        return connectivities

    def get_substructure_connectivity_counts(self, atom_connectivities_child):
        substructure_connectivity_counts = {}
        for atom_type in atom_connectivities_child:
            substructure_connectivity_counts[atom_type] = {}
            for connectivity in atom_connectivities_child[atom_type]:
                substructure_connectivity_counts[atom_type][connectivity] = 0
                for atom in self.graph:
                    if atom.type == atom_type and atom.potential_same_connectivity(connectivity):
                            substructure_connectivity_counts[atom_type][connectivity] += 1

        return substructure_connectivity_counts

    def get_substructure_connectivities(self, atom_connectivities_child):
        substructure_connectivities = {}
        for substructure_connectivity in atom_connectivities_child:
            substructure_connectivities[substructure_connectivity] = []
            for atom in self.graph:
                if atom.type != 'H' and atom.potential_same_connectivity(substructure_connectivity):
                        substructure_connectivities[substructure_connectivity].append(atom)

        
        return substructure_connectivities

    def is_substructure_atom_connectivity(self, child):
        
        atom_connectivities_child = child.get_connectivity_counts()
        atom_connectivities_self = self.get_substructure_connectivity_counts(atom_connectivities_child)

        can_be_substructure = True

        for atom_type in atom_connectivities_child:
            for connectivity in atom_connectivities_child[atom_type]:
                connectivity_nr_self = atom_connectivities_self[atom_type][connectivity]
                connectivity_nr_child = atom_connectivities_child[atom_type][connectivity]
                if connectivity_nr_child > connectivity_nr_self:
                    can_be_substructure = False
                    break

        return can_be_substructure

    def make_bond_dict(self):
        bond_dict = {}

        for atom in self.graph:
            if atom.type != 'H':
                bond_dict[atom] = 0
                for neighbour in atom.neighbours:
                    if neighbour.type != 'H':
                        bond_dict[atom] += 1
                

        return bond_dict

    def make_atom_dict(self):
        atom_dict = {}
        for atom in self.graph:
            atom_dict[atom.nr] = atom
        return atom_dict

    def make_match_dict(self):
        match_dict = {}
        for atom in self.graph:
            if atom.type != 'H':
                match_dict[atom] = None

        return match_dict

    def compress_graph(self):
        compressed_graph = {}
        for atom in self.graph:
            if atom.type != 'H':
                new_atom = copy.deepcopy(atom)
                compressed_graph[new_atom] = copy.deepcopy(self.graph[atom])

    def reset_visited(self):
        for atom in self.graph:
            atom.visited = False
        for bond_nr in self.bonds:
            bond = self.bonds[bond_nr]
            bond.visited = False

    def remove_visited(self):
        for atom in self.graph:
            delattr(atom, 'visited')
        for bond_nr in self.bonds:
            bond = self.bonds[bond_nr]
            delattr(bond, 'visited')

    def traceback(self, placed_atoms_parent, match_dict, reverse_match_dict, parent_options, atoms_to_place, child_bond_dict):

        new_parent_candidate = None
        new_child_candidate = None
        previous_child_candidate = None


        for i in range(len(placed_atoms_parent) - 1, -1, -1):
            current_atom = placed_atoms_parent[i]

            try:
                for option in parent_options[current_atom]:
                    bond = self.bond_lookup[current_atom][option]
                        
                    if not bond.visited:
                        new_parent_candidate = current_atom
                        new_child_candidate = reverse_match_dict[current_atom]
                        break
            except KeyError:
                pass

            if new_parent_candidate:
                break
            else:
                del placed_atoms_parent[i]
                child_current = reverse_match_dict[current_atom]
                atoms_to_place.add(child_current)
                    
                match_dict[child_current] = None

                child_bond_dict[child_current] += 1
                
                for bond in child_current.bonds:
                    for atom in bond.neighbours:
                        if atom != child_current:
                            if atom in match_dict:
                                child_bond_dict[atom] += 1
                    
                    
                del reverse_match_dict[current_atom]

        return new_child_candidate, new_parent_candidate

    def is_substructure_connected_atoms(self, child):

 #       parent_dfs_structure, parent_edges, parent_edge_lookup = self.make_dfs_graph()
        
        atom_connectivities_child = child.get_connectivities()
 #       atom_connectivities_parent = Structure(parent_dfs_structure, parent_edges).get_substructure_connectivities(atom_connectivities_child)
        atom_connectivities_parent = self.get_substructure_connectivities(atom_connectivities_child)
        #Sort based on the complexity of the connectivity
        
        connectivities = sorted(list(atom_connectivities_child.keys()),
                                key = lambda x: len(set(x)), reverse = True)

        starting_connectivity = connectivities[0]
        starting_atom = atom_connectivities_child[starting_connectivity][0]
        
        seeds = atom_connectivities_parent[starting_connectivity]
        matches = []

        for seed in seeds:

            self.reset_visited()

 #           parent_dfs_structure, parent_edges, parent_edge_lookup = self.make_dfs_graph()

            child_bond_dict = child.make_bond_dict()
            
            
            match_active = True
            match = child.make_match_dict()
            reverse_match = {}
            match[starting_atom] = seed
            reverse_match[seed] = starting_atom
            
            atoms_to_place = set(copy.copy(list(child.graph.keys())))
            for atom in copy.deepcopy(atoms_to_place):
                if atom.type == 'H':
                    atoms_to_place.remove(atom)
                    
            atoms_to_place.remove(starting_atom)
            
            placed_atoms_child = set([starting_atom])
            placed_atoms_parent = [seed]
            parent_options = {}
            
            current_atom_child = starting_atom
            current_atom_parent = seed

            while atoms_to_place and match_active:

                child_neighbours = current_atom_child.neighbours
                next_atom_child = None
                
                for neighbour in child_neighbours:
                    if neighbour in atoms_to_place and neighbour.type != 'H':
                        next_atom_child = neighbour
                        
                        break

                if next_atom_child:
                    child_bond = child.bond_lookup[current_atom_child][next_atom_child]
                    next_atom_parent_options = []
                    parent_neighbours = current_atom_parent.neighbours
                    for neighbour in parent_neighbours:
                        parent_bond = self.bond_lookup[current_atom_parent][neighbour]
                        
                        
                        if neighbour.type == next_atom_child.type and child_bond.type == parent_bond.type:
                            if not self.bond_lookup[current_atom_parent][neighbour].visited:
                                if neighbour.potential_same_connectivity(next_atom_child.connectivity):
                                    next_atom_parent_options.append(neighbour)

                    parent_options[current_atom_parent] = []

                    if next_atom_parent_options:
                        
                        next_atom_parent_options = sorted(next_atom_parent_options, key = lambda x: len(set(x.connectivity)), reverse = True)
                        for option in next_atom_parent_options:
                            bond = self.bond_lookup[current_atom_parent][option]
                            parent_options[current_atom_parent].append(bond)


                        next_atom_parent = next_atom_parent_options[0]
                        self.bond_lookup[current_atom_parent][next_atom_parent].visited = True
                        
                        match[next_atom_child] = next_atom_parent
                        reverse_match[next_atom_parent] = next_atom_child
                        placed_atoms_parent.append(next_atom_parent)
                        
                        child_bond_dict[current_atom_child] -= 1
                        child_bond_dict[next_atom_child] -= 1
                        
                        current_atom_child = next_atom_child
                        current_atom_parent = next_atom_parent
                        atoms_to_place.remove(current_atom_child)


                    else:
                        new_child_candidate, new_parent_candidate = self.traceback(placed_atoms_parent, match, reverse_match, parent_options, atoms_to_place, child_bond_dict)
                        if not new_child_candidate:
                            match_active = False
                        else:
                            current_atom_child = new_child_candidate
                            current_atom_parent = new_parent_candidate


                    

                else:
                    if atoms_to_place:

                        for atom in match.keys():
                            if atom != current_atom_child:
                                if child_bond_dict[atom] > 0:
                                    current_atom_child = atom
                                    current_atom_parent = match[atom]
                                    break
                    else:
                        pass

            if match_active:
                matches.append(match)

        compare_all_matches(matches)

        self.remove_visited()

        return matches
                        

    def make_atoms_upper_case(self):
        for atom in self.graph:
            if len(atom.type) == 1:
                atom.type = atom.type.upper()

    def drop_electrons(self):
        for atom in self.graph:
            atom.drop_electrons()

    def add_shells(self):
        for atom in self.graph:
            if not atom.shells:
                atom.add_shell_layout()

    def hybridise_atoms(self, atoms = None):

        if not atoms:
            for atom in self.graph:
                atom.hybridise()

        else:
            for atom in atoms:
                atom.hybridise()

    def add_hydrogens(self):
        max_atom_nr = max([atom.nr for atom in self.graph])
        max_bond_nr = max(list(self.bonds.keys()))
        
        for atom in list(self.graph.keys()):
            hydrogens_to_add = atom.calc_hydrogens()
            for i in range(hydrogens_to_add):
                max_atom_nr += 1
                max_bond_nr += 1
                hydrogen = Atom('H', max_atom_nr, None, 0)
                self.add_bond(atom, hydrogen, 'single', max_bond_nr)
                

    def refine_p_bonds(self):
        for bond_nr in self.bonds:
            bond = self.bonds[bond_nr]
            if bond.type != 'single':

                bond.combine_p_orbitals()


    def refine_s_bonds(self):
        for bond_nr in self.bonds:
            bond = self.bonds[bond_nr]
            bond.combine_hybrid_orbitals()
        
    def get_atom_counts(self):
        atom_counts = {}
        for atom in self.graph:
            if atom.type != 'H':
                if atom.type not in atom_counts:
                    atom_counts[atom.type] = 0
                atom_counts[atom.type] += 1

        return atom_counts
            

    def add_atom(self, atom):
        self.graph[atom] = []

    def sort_by_nr(self):
        for atom in self.graph:
            self.graph[atom].sort(key = lambda x: x.nr)

    def make_dummy_bond(self, atom_1, atom_2, bond_nr, dummy = False):
        if dummy:
            bond_type = 'dummy'
        else:
            bond_type = 'single'
            
        if atom_1 in self.graph:
            self.graph[atom_1].append(atom_2)
        else:
            self.graph[atom_1] = [atom_2]

        if atom_2 in self.graph:
            self.graph[atom_2].append(atom_1)
        else:
            self.graph[atom_2] = [atom_1]

        bond = Bond(atom_1, atom_2, bond_type, bond_nr)

        
        atom_1.add_bond(bond)
        atom_2.add_bond(bond)
        
        self.bonds[bond_nr] = bond

        if not atom_1 in self.bond_lookup:
            self.bond_lookup[atom_1] = {}
        if not atom_2 in self.bond_lookup:
            self.bond_lookup[atom_2] = {}

        self.bond_lookup[atom_1][atom_2] = bond
        self.bond_lookup[atom_2][atom_1] = bond



    def make_bond(self, atom_1, atom_2, bond_nr):
        if atom_1 in self.graph:
            self.graph[atom_1].append(atom_2)
        else:
            self.graph[atom_1] = [atom_2]

        if atom_2 in self.graph:
            self.graph[atom_2].append(atom_1)
        else:
            self.graph[atom_2] = [atom_1]

        bond = Bond(atom_1, atom_2, 'single', bond_nr)
        
        atom_1.add_bond(bond)
        atom_2.add_bond(bond)
        
        self.bonds[bond_nr] = bond

        electron_1 = None
        electron_2 = None

        orbital_1 = None
        orbital_2 = None

        for orbital in atom_1.valence_shell.orbitals:
            if atom_1.valence_shell.orbitals[orbital].electron_nr == 1 and not atom_1.valence_shell.orbitals[orbital].electrons[0].aromatic:
                orbital_1 = atom_1.valence_shell.orbitals[orbital]
                electron_1 = orbital_1.electrons[0]
                break
                

        for orbital in atom_2.valence_shell.orbitals:
            if atom_2.valence_shell.orbitals[orbital].electron_nr == 1 and not atom_2.valence_shell.orbitals[orbital].electrons[0].aromatic:
                orbital_2 = atom_2.valence_shell.orbitals[orbital]
                electron_2 = orbital_2.electrons[0]
                break

        orbital_1.add_electron(electron_2)
        orbital_2.add_electron(electron_1)

        orbital_1.set_bond(bond, 'sigma')
        orbital_2.set_bond(bond, 'sigma')

        if not atom_1 in self.bond_lookup:
            self.bond_lookup[atom_1] = {}
        if not atom_2 in self.bond_lookup:
            self.bond_lookup[atom_2] = {}

        self.bond_lookup[atom_1][atom_2] = bond
        self.bond_lookup[atom_2][atom_1] = bond

    def add_bond(self, atom_1, atom_2, bond_type, bond_nr):
        if atom_1 in self.graph:
            self.graph[atom_1].append(atom_2)
        else:
            self.graph[atom_1] = [atom_2]

        if atom_2 in self.graph:
            self.graph[atom_2].append(atom_1)
        else:
            self.graph[atom_2] = [atom_1]

        bond = Bond(atom_1, atom_2, bond_type, bond_nr)
        
        atom_1.add_bond(bond)
        atom_2.add_bond(bond)
        
        self.bonds[bond_nr] = bond

        if not atom_1 in self.bond_lookup:
            self.bond_lookup[atom_1] = {}
        if not atom_2 in self.bond_lookup:
            self.bond_lookup[atom_2] = {}

        self.bond_lookup[atom_1][atom_2] = bond
        self.bond_lookup[atom_2][atom_1] = bond


    def print_graph(self):
        pprint(self.graph)

    def print_bonds(self):
        pprint(self.bonds)

    def find_aromatic_bonds(self):
        aromatic_bonds = []
        for bond in self.bonds:
            if self.bonds[bond].type == 'aromatic':
                aromatic_bonds.append(self.bonds[bond])

        return aromatic_bonds

    def find_aromatic_atoms(self):
        aromatic_bonds = self.find_aromatic_bonds()
        aromatic_atoms = []
        for atom in self.graph:
            for bond in atom.bonds:
                if bond in aromatic_bonds:
                    aromatic_atoms.append(atom)

        return aromatic_atoms

    def find_aromatic_graphs(self):
        aromatic_atoms = self.find_aromatic_atoms()
        aromatic_bonds = self.find_aromatic_bonds()
        aromatic_graph = {}
        
        for atom in copy.deepcopy(self.graph):
            if atom in aromatic_atoms:
                aromatic_graph[copy.deepcopy(atom)] = copy.deepcopy(self.graph[atom])

        for atom in aromatic_graph:
            for bond in atom.bonds:
                if bond.type != 'aromatic':
                    atom.bonds.remove(bond)

        for atom_1 in aromatic_graph:
            
            for atom_2 in aromatic_graph[atom_1]:
                aromatic = False
                for bond in atom_1.bonds:
                    if bond.atom_1 == atom_2 or bond.atom_2 == atom_2:
                        aromatic = True

                if not aromatic:
                    aromatic_graph[atom_1].remove(atom_2)

        aromatic_structure = Structure(aromatic_graph, aromatic_bonds)
        
                
        aromatic_structures = aromatic_structure.split_disconnected_structures()

        return aromatic_structures
        
        

    def kekulise(self):
        aromatic_structures = self.find_aromatic_graphs()
        for aromatic_structure in aromatic_structures:
            pass
        

    def break_any_bond(self, atom_1, atom_2):
        """Remove edges from structure to break any bond

        atom_1: tuple of (str, int), with str atom type and int atom number
        atom_2: tuple of (str, int), with str atom type and int atom number

        """
        bonds_removed = 0
        while atom_2 in self.structure[atom_1]:
            self.structure[atom_1].remove(atom_2)
            bonds_removed += 1
        while atom_1 in self.structure[atom_2]:
            self.structure[atom_2].remove(atom_1)
            
        self.bonds[atom_1] -= bonds_removed
        self.bonds[atom_2] -= bonds_removed

    def remove_an_atom(self, atom):
        """Remove a node from structure to remove an atom

        atom_1: tuple of (str, int), with str atom type and int atom number
        """
        del self.strucrure[atom]

    def remove_isolated_atoms(self):
        """Remove nodes without outgoing edges from structure

        atom_1: tuple of (str, int), with str atom type and int atom number
        """
        for atom in list(self.graph.keys()):
            if not self.graph[atom]:
                del self.graph[atom]

    def split_disconnected_structures(self):
        """Return list of unconnected structures from structure


        Output:
        new_graphs: list of dicts of {node: [node, ->], ->}, with each dict a
            bidirectional graph that is not connected to the other graphs in
            the list
        """
        working_graph = copy.deepcopy(self)
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

        #add back connectors

        for new_graph in new_graphs:
            for node in new_graph:
                new_graph[node] = self.graph[node]

        #Add lone atoms
        if working_graph.graph:
            for atom in working_graph.graph:
                new_graph = {atom: []}
                if new_graph not in new_graphs:
                    new_graphs.append(new_graph)

        new_structures = []

        for new_graph in new_graphs:
            new_structures.append(Structure(new_graph))

        for new_structure in new_structures:
            new_structure.infer_bonds()
            new_structure.make_bond_lookup()

        return new_structures

    def infer_bonds(self):
        self.bonds = {}
        for atom in self.graph:
            for bond in atom.bonds:
                if not bond.nr in self.bonds:
                    self.bonds[bond.nr] = bond

    def find_lowest_atom_nr(self):
        lowest_atom_nr = 4242424242424242

        for atom in self.structure:
            if atom.nr < lowest_atom_nr:
                lowest_atom_nr = atom.nr

        return lowest_atom_nr

    def sort_substructures_by_nr(self):
        disconnected_structures = self.split_disconnected_structures()

        index_dict = {}

        for i, structure in enumerate(disconnected_structures):
            lowest_atom_nr = structure.find_lowest_atom_nr()
            index_dict[lowest_atom_nr] = i
            
        new_disconnected_structures = []

        sorted_indices = list(index_dict.keys())
        sorted_indices.sort()
        
        for index in sorted_indices:
            new_disconnected_structures += [disconnected_structures[index_dict[index]]]

        return new_disconnected_structures

    def make_sorted_copy(self):
        copy = Structure(self.structure)
        for atom in copy.structure:
            copy.structure[atom].sort()

        return copy

    def compare_structures(self, structure_2):
        sorted_copy_group_1 = self.make_sorted_copy()
        sorted_copy_group_2 = structure_2.make_sorted_copy()
        
        if sorted_copy_group_1.structure == sorted_copy_group_2.structure:
            return True
        else:
            return False

    def label_smiles_nrs(self, K_to_S_dict):
        for atom in self.structure:
            atom.set_smiles_nr(K_to_S_dict)

    def get_true_neighbour(self, atom_1, atom_2):
        bond = self.bond_lookup[atom_1][atom_2]
        true_neighbour = None
        for atom in self.graph[atom_1]:
            if atom in self.graph[atom_2] and self.bond_lookup[atom_1][atom].type != 'dummy' and  self.bond_lookup[atom_2][atom].type != 'dummy':
                true_neighbour = atom
                break

        return true_neighbour

    def structure_to_smiles(self):
        valence_dict = self.make_valence_dict()
        strings = []
        for atom in self.structure:
            atom_string = []
            H_string = self.make_H_string(atom, valence_dict)
            atom_string += [H_string]
            
            bonds = self.bond_record[atom]
            for bond in bonds:
                nr_string = self.make_nr_string(bond)
                atom_string += [nr_string]

            atom_string += ['.']
            

            strings += [''.join(atom_string)]

        string = ''.join(strings)[:-1]
        mol = Chem.MolFromSmiles(string)
        for i in range(mol.GetNumAtoms()):
            mol.GetAtomWithIdx(i).SetAtomMapNum(i)
            

        return Smiles(Chem.MolToSmiles(mol))

    def kekulise(self, aromatic_structures):
        for aromatic in aromatic_structures:
            atoms = list(aromatic.keys())
            atoms.sort()
            

            for atom in atoms:
                aromatic_bonds = []
                for bond in atom.bonds:
                    if bond.type == 'aromatic':
                        aromatic_bonds.append(bond)

                if len(aromatic_bonds) == 2:
                    pass
                    
                
        
            
        
        
                

    #========================================================================
    #Auxillary functions
    #========================================================================

    def make_H_string(self, atom, valence_dict):
        H_nr = valence_dict[atom]['max'] - valence_dict[atom]['current']
        if H_nr > 1:
            H_string = "[%sH%d]" % (atom.type, H_nr)
        elif H_nr == 1:
            H_string = "[%sH]" % (atom.type)
        else:
            H_string = "[%s]" % (atom.type)
            
        return H_string

    def make_nr_string(self, bond):
        bond_nr, bond_type = bond
        if bond_type == 1:
            string = ''
        elif bond_type == 2:
            string = '='
        elif bond_type == 3:
            string = '#'

        if bond_nr > 9:
            bond_string = '%%%d' % bond_nr
        else:
            bond_string = str(bond_nr)

        return ''.join([string, bond_string])
        
            
    def make_valence_dict(self):
        valence_dict = {}
        for atom in self.structure:
            valence_dict[atom] = {}
            valence_dict[atom]['current'] = len(self.structure[atom])
            valence_dict[atom]['max'] = atom.get_valence(len(self.structure[atom]))

        return valence_dict

    def record_bonds(self):
        bond_record = {}

        counter = 1

        self.bond_record = {}
        for atom_1 in self.structure:
            self.bond_record[atom_1] = []
            
            for atom_2 in self.structure[atom_1]:
                bond_type = self.structure[atom_1].count(atom_2)
                
                if atom_1.nr > atom_2.nr:
                    pair = (atom_2, atom_1)
                else:
                    pair = (atom_1, atom_2)

                if not pair in bond_record:
                    bond_nr = counter
                    bond_record[pair] = bond_nr
                    counter += 1
                else:
                    bond_nr = bond_record[pair]
                    
                if not (bond_nr, bond_type) in self.bond_record[atom_1]:
                    self.bond_record[atom_1] += [(bond_nr, bond_type)]
            self.bond_record[atom_1].sort()
                

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

    def make_atom_pair(self, atom_1, atom_2):

        atom_pair = [atom_1, atom_2]
        if atom_1.nr > atom_2.nr:
            atom_pair = (atom_2, atom_1)
        else:
            atom_pair = (atom_1, atom_2)

        atom_pair = tuple(atom_pair)

        return atom_pair

    def get_bond_nr(self, atom_1, atom_2):
        bond_nr = self.structure[atom_1].count(atom_2)
        return bond_nr

    def make_bond_nr_dict(self):
        """
        """
        self.bond_nr_dict = {}
        for atom in self.graph:
            self.bond_nr_dict[atom] = len(atom.bonds)


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

    def remove_connectors(self):
        """Remove nodes that only have incoming edges from graph

        Input:
        working_graph: dict of {node: [node, ->], ->}, representing a graph

        """

        for node in self.graph:
            for next_node in self.graph[node]:
                if next_node not in self.graph:
                    self.graph[node].remove(next_node)

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

class WorkingStructure(Structure):
    def __init__(self, graph = {}, bonds = {}):
        Structure.__init__(self, graph)

class RingSystem():
    def __init__(self, atoms, bonds):
        self.atoms = atoms
        self.bonds = bonds

    def get_pi_electrons(self):
        self.pi_electrons = 0
        for bond in self.bonds:
            if bond.type == 'aromatic':
                self.pi_electrons += 1
            if bond.type == 'double':
                self.pi_electrons += 2






            

        

 




        
        

if __name__ == "__main__":
    carbon = Atom('C', 1, 'counterclockwise', 0)
    carbon.add_shell_layout()
    for shell in carbon.shells:
        for orbital_set in carbon.shells[shell].orbital_sets:
            for orbital in carbon.shells[shell].orbital_sets[orbital_set].orbitals:
                print(orbital)
                print(orbital.electrons)

    print('\n')
                
    carbon.excite()
    for orbital_set in carbon.valence_shell.orbital_sets:
        for orbital in carbon.valence_shell.orbital_sets[orbital_set].orbitals:
            print(orbital)
            print(orbital.electrons)

    print('\n')

    carbon = Atom('P', 2, 'counterclockwise', 0)
    carbon.add_shell_layout()
    for shell in carbon.shells:
        for orbital_set in carbon.shells[shell].orbital_sets:
            for orbital in carbon.shells[shell].orbital_sets[orbital_set].orbitals:
                print(orbital)
                print(orbital.electrons)


 #   graph = build_graph()
#    print(graph.nodes)
#    print(graph.edges)
#    print(graph)

#    nx.draw(graph)
#    plt.show()
    
        
