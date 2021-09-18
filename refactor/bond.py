from bond_properties import BOND_PROPERTIES


class Bond:
    bond_types = {'single', 'double', 'triple', 'quadruple', 'aromatic', 'ionic', 'dummy'}

    def __init__(self, atom_1, atom_2, bond_type, bond_nr):
        atoms = [atom_1, atom_2]
        atoms.sort(key=lambda a: a.nr)

        self.atom_1 = atoms[0]
        self.atom_2 = atoms[1]
        self.neighbours = atoms

        self.chiral_symbol = None

        try:

            assert bond_type in self.bond_types

        except AssertionError:
            print(bond_type)

        self.type = bond_type
        self.nr = bond_nr
        self.aromatic = False
        if bond_type == 'aromatic':
            self.aromatic = True
        self.electrons = []

        self.chiral = False
        self.chiral_dict = {}

        if self.type == 'dummy':
            self.cbond = 0.3

        else:
            self.cbond = 0.1

        self.draw = BondDrawProperties()

    def __eq__(self, bond):
        return self.nr == bond.nr

    def __hash__(self):
        return self.nr

    def __repr__(self):
        return f'{self.type}_{self.nr}:{self.atom_1}_{self.atom_2}'

    def make_aromatic(self):
        self.type = 'aromatic'

        self.atom_1.aromatic = True
        self.atom_2.aromatic = True

        for orbital_name, orbital in self.atom_1.valence_shell.orbitals.items():
            if orbital.orbital_type == 'p':
                for electron in orbital.electrons:
                    if electron.atom != self.atom_1:
                        orbital.remove_electron(electron)
                    else:
                        electron.set_aromatic()

        for orbital_name, orbital in self.atom_2.valence_shell.orbitals.items():
            if orbital.orbital_type == 'p':
                for electron in orbital.electrons:
                    if electron.atom != self.atom_2:
                        orbital.remove_electron(electron)
                    else:
                        electron.set_aromatic()

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

        if self.atom_1.pyrrole or self.atom_2.pyrrole:
            pass
        else:
            p_bonding_orbitals_1 = []
            electrons_found = 0

            for orbital_name, orbital in self.atom_1.valence_shell.orbitals.items():
                if orbital.orbital_type == 'p' and orbital.electron_nr == 1:
                    if (orbital.electrons[0].aromatic and self.type == 'aromatic') or not orbital.electrons[0].aromatic:
                        electrons_found += 1
                        p_bonding_orbitals_1.append(orbital)

                        if electrons_found == BOND_PROPERTIES.bond_type_to_p_orbitals[self.type]:
                            break

            p_bonding_orbitals_2 = []
            electrons_found = 0

            for orbital_name, orbital in self.atom_2.valence_shell.orbitals.items():
                if orbital.orbital_type == 'p' and orbital.electron_nr == 1:
                    if (orbital.electrons[0].aromatic and self.type == 'aromatic') or not orbital.electrons[0].aromatic:
                        electrons_found += 1
                        p_bonding_orbitals_2.append(orbital)

                        if electrons_found == BOND_PROPERTIES.bond_type_to_p_orbitals[self.type]:
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


class BondDrawProperties:
    def __init__(self):
        self.center = False
