#!/usr/bin/env python
from pikachu.structure_inference.shells import Shell
from pikachu.structure_inference.bond import BOND_PROPERTIES

class AtomProperties:
    valence_dict = {'C': [4],
                    'O': [2],
                    'N': [3],
                    'P': [3, 5],
                    'S': [2, 4, 6],
                    'Cl': [1, 7],
                    'Br': [1, 7],
                    'F': [1],
                    'H': [1],
                    'I': [1, 7],
                    '*': 1}

    atomic_numbers = {'H': 1,
                      'He': 2,
                      'Li': 3,
                      'Be': 4,
                      'B': 5,
                      'C': 6,
                      'N': 7,
                      'O': 8,
                      'F': 9,
                      'Ne': 10,
                      'Na': 11,
                      'Mg': 12,
                      'Al': 13,
                      'Si': 14,
                      'P': 15,
                      'S': 16,
                      'Cl': 17,
                      'Ar': 18,
                      'K': 19,
                      'Ca': 20,
                      'Sc': 21,
                      'Ti': 22,
                      'V': 23,
                      'Cr': 24,
                      'Mn': 25,
                      'Fe': 26,
                      'Co': 27,
                      'Ni': 28,
                      'Cu': 29,
                      'Zn': 30,
                      'Ga': 31,
                      'Ge': 32,
                      'As': 33,
                      'Se': 34,
                      'Br': 35,
                      'Kr': 36,
                      'I': 53}

    atom_radii = {'H': 0.37,
                  'He': 0.32,
                  'Li': 1.34,
                  'Be': 0.90,
                  'B': 0.82,
                  'C': 0.77,
                  'N': 0.75,
                  'O': 0.73,
                  'F': 0.71,
                  'Ne': 0.69,
                  'Na': 1.54,
                  'Mg': 1.30,
                  'Al': 1.18,
                  'Si': 1.11,
                  'P': 1.06,
                  'S': 1.02,
                  'Cl': 0.99,
                  'Ar': 0.97,
                  'K': 1.96,
                  'Ca': 1.74,
                  'Sc': 1.44,
                  'Ti': 1.36,
                  'V': 1.25,
                  'Cr': 1.27,
                  'Mn': 1.39,
                  'Fe': 1.25,
                  'Co': 1.26,
                  'Ni': 1.21,
                  'Cu': 1.38,
                  'Zn': 1.31,
                  'Ga': 1.26,
                  'Ge': 1.22,
                  'As': 1.19,
                  'Se': 1.16,
                  'Br': 1.14,
                  'Kr': 1.10,
                  'I': 1.33}


    valence_electrons = {'H': 1,
                         'C': 4,
                         'O': 6,
                         'N': 5,
                         'P': 5,
                         'S': 6,
                         'Cl': 7,
                         'Br': 7,
                         'F': 7,
                         'I': 7}

    outer_shell = {'H': 1,
                   'C': 2,
                   'N': 2,
                   'O': 2,
                   'F': 2,
                   'P': 3,
                   'S': 3,
                   'Cl': 3,
                   'Br': 4,
                   'I': 5}

    electrons_per_shell = {1: 2,
                           2: 8,
                           3: 18,
                           4: 32,
                           5: 32}
                           

    available_orbitals = {1: ['1s'],
                          2: ['2s', '2p1', '2p2', '2p3'],
                          3: ['3s', '3p1', '3p2', '3p3', '3d1', '3d2', '3d3', '3d4', '3d5'],
                          4: ['4s', '4p1', '4p2', '4p3', '4d1', '4d2', '4d3', '4d4', '4d5', '4f1',
                              '4f2', '4f3', '4f4', '4f5', '4f6', '4f7'],
                          5: ['5s', '5p1', '5p2', '5p3', '5d1', '5d2', '5d3', '5d4', '5d5', '5f1',
                              '5f2', '5f3', '5f4', '5f5', '5f6', '5f7']}

    orbital_order = ('1s', '2s', '2p', '3s', '3p', '4s', '3d', '4p', '5s', '4d', '5p', '6s', '4f', '5d', '6p', '7s', '5f', '6d', '7p', '8s', '5g', '6f', '7d', '8p', '9s')

    reverse_orbital_order = reversed(orbital_order)
    
    orbital_numbers = {'s': 1,
                       'p': 3,
                       'd': 5,
                       'f': 7,
                       'g': 9}
    
ATOM_PROPERTIES = AtomProperties()

class Atom:

    def __new__(cls, atom_type, atom_nr, chiral, charge):
        self = super().__new__(cls)  # Must explicitly create the new object
        # Aside from explicit construction and return, rest of __new__
        # is same as __init__
        self.type = atom_type
        self.nr = atom_nr
        self.chiral = chiral
        self.charge = charge
        self.shells = {}
        
        return self  # __new__ returns the new object

    def __getnewargs__(self):
        # Return the arguments that *must* be passed to __new__
        return (self.type, self.nr, self.chiral, self.charge)

    def __init__(self, atom_type, atom_nr, chiral, charge):
        
        self.type = atom_type
        self.nr = atom_nr
 #       self.atomic_nr = ATOM_PROPERTIES.atomic_numbers[self.type]
        self.bonds = []

        self.chiral = chiral
        # remove?
        self.charge = charge
        self.shells = {}
        
    def __eq__(self, atom):
        return self.nr == atom.nr

    def __hash__(self):
        return self.nr

    def __repr__(self):
        if self.charge == 0:
            charge_string = ''
        elif self.charge > 0:
            if self.charge == 1:
                charge_string = '+'
            else:
                charge_string = str(self.charge) + '+'
        elif self.charge < 0:
            if self.charge == -1:
                charge_string = '-'
            else:
                charge_string = str(abs(self.charge)) + '-'
        
        return f'{self.type}{charge_string}_{self.nr}'

    def set_neighbours(self, structure):
        self.neighbours = structure.graph[self]

    def remove_neighbour(self, neighbour):
        self.neighbours.remove(neighbour)

    def set_connectivity(self):
        self.connectivity = self.get_connectivity()

    def get_connectivity(self):
        connectivity = []

        for bond in self.bonds:
            for atom in bond.neighbours:
                if atom.type != 'H' and atom != self:
                    bond_type = bond.type
                    connectivity.append(f'{atom.type}_{bond_type}')

        connectivity = tuple(sorted(connectivity))
        return connectivity

    def same_connectivity(self, atom):
        if self.type == atom.type:
            if len(self.connectivity) == len(atom.connectivity):
                if self.chiral == atom.chiral == None:
                    if set(self.connectivity) == set(atom.connectivity):
                        return True
                    else:
                        return False
                else:
                    pass
            else:
                return False
                

        else:
            return False

    def potential_same_connectivity(self, substructure_connectivity):
        parent_connectivity_copy = list(copy.copy(self.connectivity))
        substructure_connectivity_copy = list(copy.copy(substructure_connectivity))

        same_connectivity = True

        for atom in substructure_connectivity_copy:
            if atom in parent_connectivity_copy:
                parent_connectivity_copy.remove(atom)
            else:
                same_connectivity = False

        return same_connectivity

    def set_order(self):
        self.order = 0
        for neighbour in self.neighbours:
            if neighbour.type != 'H':
                self.order += 1
            
            

    def add_shell_layout(self):
        self.shell_nr = ATOM_PROPERTIES.outer_shell[self.type]
        self.make_shells()
        self.fill_shells()

        self.excitable = self.valence_shell.is_excitable()

        if self.type == 'C':
            self.excite()
        else:
            bond_weights = [BOND_PROPERTIES.bond_weights[bond.type] for bond in self.bonds]
            nr_of_nonH_bonds = sum(bond_weights)
            if nr_of_nonH_bonds > ATOM_PROPERTIES.valence_dict[self.type][0]:
                self.excite()

    def make_shells(self):
        for i in range(self.shell_nr):
            current_shell = i + 1
            self.shells[current_shell] = Shell(self, current_shell)

        self.valence_shell = self.shells[self.shell_nr]
            

    def fill_shells(self):
        electrons_assigned = 0

        electrons_remaining = ATOM_PROPERTIES.atomic_numbers[self.type] - self.charge
        

        #Iterate over the orbitals in order of them being filled

        
        for orbital in ATOM_PROPERTIES.orbital_order:
            if electrons_remaining > 0:
                shell = int(orbital[0])
                orbital_set = self.shells[shell].orbital_sets[orbital]
                electrons_to_dump = min([electrons_remaining, orbital_set.capacity])
                orbital_set.fill_orbitals(electrons_to_dump)
                electrons_remaining -= electrons_to_dump
            else:
                break

    def excite(self):
        assert self.excitable
        
        self.valence_shell.excite()

    def remove_bond(self, bond):
        self.bonds.remove(bond)

    def calc_electron_pair_nr(self):

        bond_nr = self.calc_bond_nr()
        bonds_accounted_for = 0
        
        electron_nr = 0
        orbital_nr = len(list(self.valence_shell.orbitals.keys()))
        
        for orbital_name in self.valence_shell.orbitals:
            orbital = self.valence_shell.orbitals[orbital_name]
            if orbital.electron_nr == 1:
                electron_nr += 1
            elif orbital.electron_nr == 2:
                if orbital.electrons[0].atom == orbital.electrons[1].atom:
                    electron_nr += 2
                else:
                    bonds_accounted_for += 1

        bonds_to_make = bond_nr - bonds_accounted_for

        unbonded_electrons = electron_nr - bonds_to_make

        

        if unbonded_electrons % 2 != 0:
            print("Warning! Rogue electron.")
            print(self.valence_shell.orbitals)

        electron_pair_nr = int(unbonded_electrons / 2)
        
        return electron_pair_nr

    def drop_electrons(self):
        if self.valence_shell.get_lone_electrons() > 1:
            self.valence_shell.drop_electrons()
        

    def calc_bond_nr(self):

        bond_nr = 0
        aromatic_bond_nr = 0

        for bond in self.bonds:
            if bond.type == 'single':
                bond_nr += 1
            elif bond.type == 'double':
                bond_nr += 2
            elif bond.type == 'triple':
                bond_nr += 3
            elif bond.type == 'quadruple':
                bond_nr += 4
            elif bond.type == 'aromatic':
                aromatic_bond_nr += 1

        if aromatic_bond_nr == 2:
            bond_nr += 3
        elif aromatic_bond_nr == 3:
            bond_nr += 4

        return bond_nr

    def calc_CIP_priority(self):

        atom = self
        cip = CIP(self)
        

        if not cip.ties:
            return cip.priority_list

        else:
            pass


    def is_promotable(self):
        promotable = False
        for orbital_set in self.valence_shell.orbital_sets:
            if 'd' in orbital_set:
                promotable = True

        return promotable

    def promote_lone_pair_to_p_orbital(self):
        assert self.calc_electron_pair_nr() > 0
        assert self.hybridisation == 'sp3'

        self.valence_shell.dehybridise()
        #self.valence_shell.hybridise('sp2')

        p_orbitals = []
        sp2_orbitals = []

        for orbital_name in self.valence_shell.orbitals:
            orbital = self.valence_shell.orbitals[orbital_name]
            if orbital.electron_nr == 2:
                if orbital.electrons[0].atom != orbital.electrons[1].atom:
                    sp2_orbitals.append(orbital)
                elif orbital.electrons[0].atom == orbital.electrons[1].atom == self:
                    p_orbitals.append(orbital)

        if len(p_orbitals) > 1:
            for i in range(1, len(p_orbitals)):
                sp2_orbitals.append(p_orbitals[i])

        p_orbital = p_orbitals[0]

        p_orbital.orbital_type = 'p'

        for orbital in sp2_orbitals:
            orbital.orbital_type = 'sp2'

        self.hybridisation = 'sp2'

        for orbital_name in self.valence_shell.orbitals:
            orbital = self.valence_shell.orbitals[orbital_name]
            for electron in orbital.electrons:
                if electron.atom == self:
                    electron.set_orbital(orbital)
            

    def promote_pi_bond_to_d_orbital(self):
        assert self.is_promotable()
            
        donor_orbitals = []
        receiver_orbitals = []
        for orbital_name in self.valence_shell.orbitals:
            orbital = self.valence_shell.orbitals[orbital_name]
            if orbital.orbital_type == 'p' and orbital.electron_nr == 2:
                if orbital.electrons[0].atom != orbital.electrons[1].atom:
                    donor_orbitals.append(orbital)

            elif orbital.orbital_type == 'd' and orbital.electron_nr == 1:
                receiver_orbitals.append(orbital)

        donor_orbital = donor_orbitals[0]
        receiver_orbital = receiver_orbitals[0]

        moved_electron = None

        for electron in donor_orbital.electrons:
            if electron.atom != self:
                moved_electron = electron


        donor_orbital.remove_electron(moved_electron)
        receiver_orbital.add_electron(moved_electron)

        receiver_orbital.set_bond(donor_orbital.bond, 'pi')
        donor_orbital.remove_bond()

        self.valence_shell.dehybridise()

        self.hybridise()


    def calc_hydrogens(self):
        hydrogens = 0
        if self.type in ['B', 'C', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I']:

            bond_nr = self.calc_bond_nr()
            if bond_nr in ATOM_PROPERTIES.valence_dict[self.type]:
                hydrogens = 0
            else:
                max_bonds = self.valence_shell.get_lone_electrons()
                hydrogens = max_bonds - bond_nr

        return hydrogens

    def add_bond(self, bond):
        self.bonds.append(bond)

    def calc_repulsion_force(self, atom, crep, repulsion_threshold):
        self_coords = self.get_coords()
        atom_coords = atom.get_coords()
        
        distance = squared_distance(self_coords, atom_coords)

        if repulsion_threshold > distance > 0.1:
            vector = []
            
            for i, coord in enumerate(self_coords):
                atom_coord = atom_coords[i]
                diff = atom_coord - coord
                vector.append(crep * diff / distance)

        elif distance <= 0.1:
            vector = force_to_vector(1.0, atom, self)
 #           print('stom_1', self, self.get_coords())
#            print('atom_2', atom, atom.get_coords())
#            print(vector)

        else:
            vector = [0.0, 0.0, 0.0]

        return vector

    def hybridise(self):
        hybridisation = self.get_hybridisation()
        self.valence_shell.hybridise(hybridisation)
        self.set_hybridisation()

    def set_hybridisation(self):
        self.hybridisation = 's'
        for orbital_name in self.valence_shell.orbitals:
            orbital = self.valence_shell.orbitals[orbital_name]
            if orbital.orbital_type in {'sp', 'sp2', 'sp3', 'sp3d', 'sp3d2'}:
                self.hybridisation = orbital.orbital_type
                break

    def get_hybridisation(self):
        steric_number = self.get_steric_number()
        #Make dict

        if steric_number == 1:
            hybridisation = 's'
        elif steric_number == 2:
            hybridisation = 'sp'
        elif steric_number == 3:
            hybridisation = 'sp2'
        elif steric_number == 4:
            hybridisation = 'sp3'
        elif steric_number == 5:
            hybridisation = 'sp3d'
        elif steric_number == 6:
            hybridisation = 'sp3d2'
        

        return hybridisation

    def get_steric_number(self):
        return self.calc_electron_pair_nr() + len(self.bonds) 

    def get_valence(self):
        try:
            valence = ATOM_PROPERTIES.valence_dict[self.type][0]
            return valence

        except KeyError:
            print("Unknown atom: %s" % self.type)
            return "Unknown atom: %s" % self.type

    def get_coords(self):
        return [self.x, self.y, self.z]

    def get_bond_angle(self, structure):
        non_dummy = []
        for atom in structure.graph[self]:
            if structure.bond_lookup[self][atom].type != 'dummy':
                non_dummy.append(atom)

        if len(non_dummy) == 3:
            angle = 120
        elif len(non_dummy) == 4:
            angle = 90
        elif len(non_dummy) == 2:
            angle = 180
        elif len(non_dummy) == 1:
            angle = 0
        else:
            angle = None

        return angle

    def move_atom(self, vector):
        current_position = self.get_coords()
        new_position = [self.x + vector[0],
                        self.y + vector[1],
                        self.z + vector[2]]
        distance = euclidean(current_position, new_position)
        if distance > 1.5:
            x_displacement = vector[0] * 1.5 / distance
            y_displacement = vector[1] * 1.5 / distance
            z_displacement = vector[2] * 1.5 / distance
        else:
            x_displacement = vector[0]
            y_displacement = vector[1]
            z_displacement = vector[2]

        self.x += x_displacement
        self.y += y_displacement
        self.z += z_displacement
            
                

    def calc_distance(self, atom):
        distance = squared_distance(self.get_coords(), atom.get_coords())
        return distance
