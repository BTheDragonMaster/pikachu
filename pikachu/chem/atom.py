import copy

from pikachu.chem.lone_pair import LonePair
from pikachu.chem.atom_properties import ATOM_PROPERTIES
from pikachu.chem.bond_properties import BOND_PROPERTIES
from pikachu.errors import SmilesError
from pikachu.chem.shell import Shell
from pikachu.math_functions import Vector


class Atom:

    def __new__(cls, atom_type, atom_nr, chiral, charge, aromatic):
        self = super().__new__(cls)  # Must explicitly create the new object
        # Aside from explicit construction and return, rest of __new__
        # is same as __init__
        self.type = atom_type
        self.nr = atom_nr
        self.chiral = chiral
        self.charge = charge
        self.aromatic = aromatic
        self.shells = {}

        return self  # __new__ returns the new object

    def __getnewargs__(self):
        # Return the arguments that *must* be passed to __new__
        return (self.type, self.nr, self.chiral, self.charge, self.aromatic)

    def __init__(self, atom_type, atom_nr, chiral, charge, aromatic):

        self.type = atom_type
        self.nr = atom_nr
        self.aromatic = aromatic
        # self.atomic_nr = ATOM_PROPERTIES.element_to_atomic_nr[self.type]
        self.bonds = []

        self.chiral = chiral
        self.charge = charge
        self.pyrrole = False
        self.shells = {}
        self.lone_pairs = []
        self.draw = AtomDrawProperties()

    def __eq__(self, atom):
        if type(atom) == Atom:
            return self.nr == atom.nr
        else:
            return False

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
        else:
            if self.charge == -1:
                charge_string = '-'
            else:
                charge_string = str(abs(self.charge)) + '-'

        return f'{self.type}{charge_string}_{self.nr}'

    def make_lone_pairs(self):

        lone_pair_nr = self.valence_shell.get_lone_pair_nr()
        for i in range(lone_pair_nr):
            self.lone_pairs.append(LonePair(self, i + 10000))

    def get_bonds(self):
        return self.bonds[:]

    def set_neighbours(self, structure):
        self.neighbours = structure.graph[self]

    def set_drawn_neighbours(self):
        self.drawn_neighbours = []
        for neighbour in self.neighbours:
            if neighbour.draw.is_drawn:
                self.drawn_neighbours.append(neighbour)

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
                if set(self.connectivity) == set(atom.connectivity):
                    return True
                else:
                    return False

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

    def get_drawn_neighbours(self):
        drawn_neighbours = []
        for neighbour in self.neighbours:
            if neighbour.draw.is_drawn:
                drawn_neighbours.append(neighbour)

        return drawn_neighbours
    
    def has_neighbour(self, atom_type):
        for neighbour in self.neighbours:
            if neighbour.type == atom_type:
                return True
        
        return False

    def add_shell_layout(self):
        self.shell_nr = ATOM_PROPERTIES.element_to_shell_nr[self.type]
        self.make_shells()
        self.fill_shells()

        self.excitable = self.valence_shell.is_excitable()

        # Can't excite carbon if it has more than 4 electrons

        if (self.type == 'C' or self.type == 'B') and self.charge >= 0:
            self.excite()
        else:

            # Assigns a value to all bonds: 2 for double bond, 1 for single, ~1.5 for aromatic
            bond_weights = [BOND_PROPERTIES.bond_type_to_weight[bond.type] for bond in self.bonds]
            aromatic_count = 0
            for bond in self.bonds:
                if bond.type == 'aromatic':
                    aromatic_count += 1

            # If odd number of aromatic bonds (such as central atoms in Trp), only add 1 'extra' bond for the
            # three outgoing aromatic bonds

            nr_of_nonH_bonds = sum(bond_weights) + int(aromatic_count / 2)

            if self.pyrrole:
                nr_of_nonH_bonds -= 1

            # Does this work for all atoms? Doesn't for carbon. Should this be made general?

            bonding_electrons = self.get_bonding_electrons()

            if nr_of_nonH_bonds > bonding_electrons:
                if self.excitable:
                    self.excite()
                else:
                    raise SmilesError('violated_bonding_laws')

        #    if nr_of_nonH_bonds > ATOM_PROPERTIES.element_to_valences[self.type][0] + self.charge:
         #       if self.excitable:
          #          print("exciting", self)
          #          self.excite()
           #     else:
            #        raise SmilesError('violated_bonding_laws')

    def get_bonding_electrons(self):
        counter = 0
        for orbital_name, orbital in self.valence_shell.orbitals.items():
            if len(orbital.electrons) == 1:
                counter += 1
        return counter

    def make_shells(self):
        for i in range(self.shell_nr):
            current_shell = i + 1
            self.shells[current_shell] = Shell(self, current_shell)

        self.valence_shell = self.shells[self.shell_nr]

    def in_ring(self, structure):
        cycles = structure.cycles.all_cycles

        for cycle in cycles:
            if self in cycle:
                return True

        return False

    def fill_shells(self):
        electrons_assigned = 0

        electrons_remaining = ATOM_PROPERTIES.element_to_atomic_nr[self.type] - self.charge

        # Iterate over the orbitals in order of them being filled

        for orbital in ATOM_PROPERTIES.orbital_order:
            if electrons_remaining > 0:
                shell = int(orbital[0])
                orbital_set = self.shells[shell].orbital_sets[orbital]
                electrons_to_dump = min([electrons_remaining, orbital_set.capacity])
                orbital_set.fill_orbitals(electrons_to_dump)
                electrons_remaining -= electrons_to_dump
            else:
                break
                
    def get_neighbour(self, atom_type):
        for neighbour in self.neighbours:
            if neighbour.type == atom_type:
                return neighbour
        return None

    def excite(self):
        assert self.excitable

        self.valence_shell.excite()

    def get_non_hydrogen_neighbours(self):
        neighbours = []
        for atom in self.neighbours:
            if atom.type != 'H' and atom.type != '*':
                neighbours.append(atom)
        return neighbours

    def get_non_hydrogen_bonds(self):
        bonds = []
        for bond in self.bonds:
            if bond.atom_1.type != 'H' and bond.atom_2.type != 'H':
                bonds.append(bond)
        return bonds

    def remove_bond(self, bond):
        self.bonds.remove(bond)

    def calc_electron_pair_nr(self):

        bond_nr = self.calc_bond_nr()
        bonds_accounted_for = 0

        electron_nr = 0
        orbital_nr = len(list(self.valence_shell.orbitals.keys()))

        for orbital_name, orbital in self.valence_shell.orbitals.items():

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
            print(self)
            print(bond_nr)
            print(self.bonds)
            print(bonds_to_make)
            print(bonds_accounted_for)
            print(electron_nr)

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
            if not self.pyrrole:
                bond_nr += 3
            else:
                bond_nr += 2
        elif aromatic_bond_nr == 3:
            bond_nr += 4

        return bond_nr

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
        # self.valence_shell.hybridise('sp2')

        p_orbitals = []
        sp2_orbitals = []

        for orbital_name, orbital in self.valence_shell.orbitals.items():
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

        for orbital_name, orbital in self.valence_shell.orbitals.items():
            for electron in orbital.electrons:
                if electron.atom == self:
                    electron.set_orbital(orbital)
                    if orbital.orbital_type == 'p':
                        electron.set_aromatic()
                        
    def promote_pi_bonds_to_d_orbitals(self):

        if self.is_promotable() and 'd' in self.hybridisation:
            
            donor_orbitals = []
            receiver_orbitals = []
            for orbital_name in self.valence_shell.orbitals:
                orbital = self.valence_shell.orbitals[orbital_name]
                if 'p' in orbital.orbital_type and orbital.orbital_type != 'p' and orbital.electron_nr == 2:
                    if orbital.electrons[0].atom != orbital.electrons[1].atom:
                        donor_orbitals.append(orbital)
    
                elif orbital.orbital_type == 'd' and orbital.electron_nr == 1:
                    receiver_orbitals.append(orbital)

            if donor_orbitals and receiver_orbitals:
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

    def reset_hybridisation(self):
        self.valence_shell.dehybridise()
        self.hybridise()

    def calc_hydrogens(self):
        hydrogens = 0
        if self.type in ['B', 'C', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I']:

            bond_nr = self.calc_bond_nr()
            if bond_nr in ATOM_PROPERTIES.element_to_valences[self.type]:
                hydrogens = 0
            else:
                max_bonds = self.valence_shell.get_lone_electrons()
                hydrogens = max_bonds - bond_nr

        return hydrogens

    def add_bond(self, bond):
        self.bonds.append(bond)

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
        # Make dict

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
        elif steric_number == 0:
            hybridisation = None
        else:
            hybridisation = None

        return hybridisation

    def get_steric_number(self):
        return self.calc_electron_pair_nr() + len(self.bonds)

    def get_valence(self):
        try:
            valence = ATOM_PROPERTIES.element_to_valences[self.type][0]
            return valence

        except KeyError:
            print("Unknown atom: %s" % self.type)
            return "Unknown atom: %s" % self.type

    def get_coords(self):
        return [self.x, self.y, self.z]

    def get_hydrogen_nr(self, structure):
        hydrogen_count = 0
        for atom in structure.graph[self]:
            if atom.type == 'H':
                hydrogen_count += 1

        return hydrogen_count


class AtomDrawProperties:
    def __init__(self, x=0, y=0):
        self.rings = []
        self.original_rings = []
        self.anchored_rings = []
        self.is_bridge_atom = False
        self.is_bridge = False
        self.bridged_ring = None
        self.is_drawn = True
        self.has_hydrogen = False
        self.positioned = False
        self.previous_position = Vector(0, 0)
        self.position = Vector(x, y)
        self.angle = None
        self.force_positioned = False
        self.connected_to_ring = False
        self.draw_explicit = False
        self.previous_atom = None
        self.colour = 'black'

    def set_position(self, vector):
        self.position = vector

    def set_previous_position(self, previous_atom):
        self.previous_position = previous_atom.draw.position
        self.previous_atom = previous_atom

    def get_angle(self, reference_vector=None):
        if not reference_vector:
            vector = Vector.subtract_vectors(self.position, self.previous_position)
        else:
            vector = Vector.subtract_vectors(self.position, reference_vector)

        return vector.angle()

    def restore_rings(self):
        self.rings = []
        for ring in self.original_rings:
            self.rings.append(ring)