class Bond:
    bond_types = {'single', 'double', 'triple', 'quadruple', 'aromatic', 'ionic', 'dummy'}
    
                 
    def __init__(self, atom_1, atom_2, bond_type, bond_nr):
        atoms = [atom_1, atom_2]
        atoms.sort(key = lambda a: a.nr)

        self.atom_1 = atoms[0]
        self.atom_2 = atoms[1]
        self.neighbours = atoms

        assert bond_type in self.bond_types
        
        self.type = bond_type
        self.nr = bond_nr
        self.aromatic = False
        if bond_type == 'aromatic':
            self.aromatic = True
        self.electrons = []

        self.chiral = False

        if self.type == 'dummy':
            self.cbond = 0.3

        else:
            self.cbond = 0.1

    def __eq__(self, bond):
        return self.nr == bond.nr

    def __hash__(self):
        return self.nr

    def __repr__(self):
        return f'{self.type}_{self.nr}:{self.atom_1}_{self.atom_2}'

    def set_double_bond_chirality(self, atom_1, atom_2, orientation):

        self.chiral_dict = {}
        
        self.orientation = orientation
        side_1 = None
        side_2 = None
        for atom in self.neighbours:
            if atom_1 in atom.neighbours:
                side_1 = atom
            elif atom_2 in atom.neighbours:
                side_2 = atom

        atom_1_2 = None
        atom_2_2 = None

        for atom in side_1.neighbours:
            if atom != side_2 and atom != atom_1:
                atom_1_2 = atom

        for atom in side_2.neighbours:
            if atom != side_1 and atom != atom_2:
                atom_2_2 = atom

        self.chiral_dict[atom_1] = {}
        self.chiral_dict[atom_2] = {}
        self.chiral_dict[atom_1_2] = {}
        self.chiral_dict[atom_2_2] = {}

        opposite_orientation = None
        
        if orientation == 'cis':
            opposite_orientation = 'trans'
        else:
            opposite_orientation = 'cis'

        self.chiral_dict[atom_1][atom_2] = orientation
        self.chiral_dict[atom_1][atom_2_2] = opposite_orientation
        
        self.chiral_dict[atom_2][atom_1] = orientation
        self.chiral_dict[atom_2][atom_1_2] = opposite_orientation
        
        self.chiral_dict[atom_1_2][atom_2] = opposite_orientation
        self.chiral_dict[atom_1_2][atom_2_2] = orientation
        
        self.chiral_dict[atom_2_2][atom_1] = opposite_orientation
        self.chiral_dict[atom_2_2][atom_1_2] = opposite_orientation

        self.chiral = True
            
        

        
    def check_same_chirality(self, parent_bond, match):
        same_chirality = True
        for atom in self.chiral_dict:
            if atom.type != 'H':
                parent_atom = match[atom]
                for atom_2 in self.chiral_dict[atom]:
                    if atom_2.type != 'H':
                        parent_atom_2 = match[atom_2]
                        orientation = self.chiral_dict[atom][atom_2]
                        parent_orientation = parent_bond.chiral_dict[parent_atom][parent_atom_2]
                        if orientation != parent_orientation:
                            same_chirality = False
                            break
            if not same_chirality:
                break

        return same_chirality
        
        
            
                

    def break_bond(self):
        assert self.type == 'single'

        electron_1, electron_2 = self.electrons
        orbital_1 = electron_1.orbital
        orbital_2 = electron_2.orbital

        orbital_1.remove_electron(electron_2)
        orbital_2.remove_electron(electron_1)
        
        orbital_1.remove_bond()
        orbital_2.remove_bond()
        
        self.atom_1.remove_neighbour(self.atom_2)
        self.atom_2.remove_neighbour(self.atom_1)
        
        self.atom_1.remove_bond(self)
        self.atom_2.remove_bond(self)

        
    def combine_hybrid_orbitals(self):

        s_bonding_orbital_1 = None
        s_bonding_orbital_2 = None
        
        for orbital_name in self.atom_1.valence_shell.orbitals:
            orbital = self.atom_1.valence_shell.orbitals[orbital_name]
            if 's' in orbital.orbital_type and orbital.electron_nr == 1:
                s_bonding_orbital_1 = orbital

        if not s_bonding_orbital_1:
            promotable = self.atom_1.is_promotable()
            if promotable:
                self.atom_1.promote_pi_bond_to_d_orbital()
                for orbital_name in self.atom_1.valence_shell.orbitals:
                    orbital = self.atom_1.valence_shell.orbitals[orbital_name]
                    if 's' in orbital.orbital_type and orbital.electron_nr == 1:
                        s_bonding_orbital_1 = orbital
                    

        for orbital_name in self.atom_2.valence_shell.orbitals:
            orbital = self.atom_2.valence_shell.orbitals[orbital_name]
            if 's' in orbital.orbital_type and orbital.electron_nr == 1:
                s_bonding_orbital_2 = orbital

        if not s_bonding_orbital_2:
            promotable = self.atom_2.is_promotable()
            if promotable:
                self.atom_2.promote_pi_bond_to_d_orbital()
                for orbital_name in self.atom_2.valence_shell.orbitals:
                    orbital = self.atom_2.valence_shell.orbitals[orbital_name]
                    if 's' in orbital.orbital_type and orbital.electron_nr == 1:
                        s_bonding_orbital_2 = orbital

        electron_1 = s_bonding_orbital_1.electrons[0]
        electron_2 = s_bonding_orbital_2.electrons[0]

        self.electrons.append(electron_1)
        self.electrons.append(electron_2)

        s_bonding_orbital_1.add_electron(electron_2)
        s_bonding_orbital_2.add_electron(electron_1)

        s_bonding_orbital_1.set_bond(self, 'sigma')
        s_bonding_orbital_2.set_bond(self, 'sigma')

    def combine_p_orbitals(self):
        assert self.type != 'single'

        p_bonding_orbitals_1 = []
        electrons_found = 0
        
        for orbital_name in self.atom_1.valence_shell.orbitals:
            orbital = self.atom_1.valence_shell.orbitals[orbital_name]
            if orbital.orbital_type == 'p' and orbital.electron_nr == 1:
                if (orbital.electrons[0].aromatic and self.type == 'aromatic') or not orbital.electrons[0].aromatic:
                    electrons_found += 1
                    p_bonding_orbitals_1.append(orbital)

                    if electrons_found == BOND_PROPERTIES.p_bonding_orbitals[self.type]:
                        break


        p_bonding_orbitals_2 = []
        electrons_found = 0
                
        for orbital_name in self.atom_2.valence_shell.orbitals:
            orbital = self.atom_2.valence_shell.orbitals[orbital_name]
            if orbital.orbital_type == 'p' and orbital.electron_nr == 1:
                if (orbital.electrons[0].aromatic and self.type == 'aromatic') or not orbital.electrons[0].aromatic:
                    electrons_found += 1
                    p_bonding_orbitals_2.append(orbital)
                    
                    if electrons_found == BOND_PROPERTIES.p_bonding_orbitals[self.type]:
                        break


        assert len(p_bonding_orbitals_1) == len(p_bonding_orbitals_2)

        if self.type == 'aromatic':
            assert len(p_bonding_orbitals_1) == len(p_bonding_orbitals_2) == 1
            electron_1 = p_bonding_orbitals_1[0].electrons[0]
            electron_2 = p_bonding_orbitals_2[0].electrons[0]
            
            electron_1.set_aromatic()
            electron_2.set_aromatic()

            self.electrons.append(electron_1)
            self.electrons.append(electron_2)
            
            p_bonding_orbitals_1[0].set_bond(self, 'pi')
            p_bonding_orbitals_2[0].set_bond(self, 'pi')
        

        else:

            for i in range(len(p_bonding_orbitals_1)):
                electron_1 = p_bonding_orbitals_1[i].electrons[0]
                electron_2 = p_bonding_orbitals_2[i].electrons[0]
                

                p_bonding_orbitals_1[i].add_electron(electron_2)
                p_bonding_orbitals_2[i].add_electron(electron_1)

                self.electrons.append(electron_1)
                self.electrons.append(electron_2)
                
                p_bonding_orbitals_1[i].set_bond(self, 'pi')
                p_bonding_orbitals_2[i].set_bond(self, 'pi')



    def calc_bond_spring_force(self):
        #calculate the distance between a bond's two neighbours
        r = self.calc_bond_length_model()
        bond_spring_force = self.cbond * (self.blength - r)
        return bond_spring_force

    def calc_angle_force(self, structure):
        r = self.calc_bond_length_model()
        
        true_neighbour = structure.get_true_neighbour(self.atom_1, self.atom_2)
        bond_angle = true_neighbour.get_bond_angle(structure)
        double_repulsion = False

        if bond_angle == 120:
            angle_force = 0.3 * (sqrt(3) * self.blength - r)
        elif bond_angle == 180:
            angle_force = 0.3 * (0.65 * blength - r)
        elif bond_angle == 90:
            angle_force = 0
            double_repulsion = True
        else:
            angle_force = 0

        return angle_force, double_repulsion
            
        
    def calc_bond_length_model(self):
        r = self.atom_1.calc_distance(self.atom_2)
        return r

    def calc_bond_length(self):
        atom_radii = ATOM_PROPERTIES.atom_radii
        atom_1_radius = atom_radii[atom_1.type]
        atom_2_radius = atom_radii[atom_2.type]
        bond_length = atom_1_radius + atom_2_radius
        return bond_length

class BondProperties():

    bond_weights = {'single': 1,
                    'double': 2,
                    'triple': 3,
                    'quadruple': 4,
                    'aromatic': 1}
    p_bonding_orbitals = {'single': 0,
                          'double': 1,
                          'triple': 2,
                          'quadruple': 3,
                          'aromatic': 1}

BOND_PROPERTIES = BondProperties()   

