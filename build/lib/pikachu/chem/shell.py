from pikachu.chem.orbital import OrbitalSet
from pikachu.errors import StructureError
from pikachu.chem.atom_properties import ATOM_PROPERTIES


class Shell:

    def __init__(self, atom, shell_nr):
        self.shell_nr = shell_nr
        self.orbital_sets = {}
        self.orbitals = {}
        self.bonding_orbitals = []
        self.atom = atom

        self.define_orbitals()

        self.find_bonding_orbitals()

    def define_orbitals(self):
        self.orbitals = []
        self.orbital_sets[f'{self.shell_nr}s'] = OrbitalSet(self.atom, self.shell_nr, 's')
        if self.shell_nr >= 2:
            self.orbital_sets[f'{self.shell_nr}p'] = OrbitalSet(self.atom, self.shell_nr, 'p')
        if self.shell_nr >= 3:
            self.orbital_sets[f'{self.shell_nr}d'] = OrbitalSet(self.atom, self.shell_nr, 'd')
        if self.shell_nr >= 4:
            self.orbital_sets[f'{self.shell_nr}f'] = OrbitalSet(self.atom, self.shell_nr, 'f')

        for orbital_set in self.orbital_sets:
            for orbital in self.orbital_sets[orbital_set].orbitals:
                self.orbitals.append(orbital)

    def __hash__(self):
        return f'{self.atom.nr}_{self.shell_nr}'

    def __repr__(self):
        return f'{self.atom.nr}_{self.shell_nr}'

    def hybridise(self, hybridisation):
        if hybridisation == 'sp3':
            self.sp_hybridise(3)
        elif hybridisation == 'sp2':
            self.sp_hybridise(2)
        elif hybridisation == 'sp':
            self.sp_hybridise(1)
        elif hybridisation == 'sp3d':
            self.spd_hybridise(1)
        elif hybridisation == 'sp3d2':
            self.spd_hybridise(2)
        elif hybridisation is None:
            pass

        for orbital in self.orbitals:
            for electron in orbital.electrons:
                if self.atom == electron.atom:
                    electron.set_orbital(orbital)

    def count_p_orbitals(self):
        count = 0
        for orbital in self.orbitals:
            if orbital.orbital_type == 'p':
                count += 1

        return count

    def count_d_orbitals(self):
        count = 0
        for orbital in self.orbitals:
            if orbital.orbital_type == 'd':
                count += 1

        return count

    def dehybridise(self):
        for orbital_set in self.orbital_sets:
            for i, orbital in enumerate(self.orbital_sets[orbital_set].orbitals):
                if orbital.orbital_type not in {'s', 'p', 'd', 'f'}:
                    new_orbital_type = self.orbital_sets[orbital_set].orbital_type
                    if new_orbital_type != 's':
                        new_orbital_nr = i + 1
                    else:
                        new_orbital_nr = None

                    orbital.orbital_type = new_orbital_type
                    orbital.orbital_nr = new_orbital_nr

        for orbital in self.orbitals:
            for electron in orbital.electrons:
                if self.atom == electron.atom:
                    electron.set_orbital(orbital)

    def sp_hybridise(self, p_nr):
        if p_nr == 1:
            orbital_type = 'sp'
        else:
            orbital_type = f'sp{p_nr}'
        hybridised_p = 0

        orbital_nr = 1

        for orbital in self.orbitals:
            if orbital.orbital_type == 's':
                orbital.orbital_nr = orbital_nr
                orbital.orbital_type = orbital_type
                orbital_nr += 1
            elif orbital.orbital_type == 'p':
                if not orbital.bond or orbital.bonding_orbital == 'sigma':
                    if hybridised_p < p_nr:
                        orbital.orbital_type = orbital_type
                        orbital.orbital_nr = orbital_nr
                        hybridised_p += 1
                        orbital_nr += 1

    def spd_hybridise(self, d_nr):
        hybridised_d = 0

        orbital_nr = 1

        if d_nr == 1:
            orbital_type = 'sp3d'
        else:
            orbital_type = f'sp3d{d_nr}'

        for orbital in self.orbitals:
            if orbital.orbital_type == 's':
                orbital.orbital_type = orbital_type
                orbital.orbital_nr = orbital_nr
                orbital_nr += 1
            if orbital.orbital_type == 'p':
                orbital.orbital_type = orbital_type
                orbital.orbital_nr = orbital_nr
                orbital_nr += 1
            elif orbital.orbital_type == 'd':
                if not orbital.bond or orbital.bonding_orbital == 'sigma':
                    if hybridised_d < d_nr:
                        orbital.orbital_type = orbital_type
                        hybridised_d += 1
                        orbital.orbital_nr = orbital_nr
                        orbital_nr += 1

    def excite(self):
        try:
            assert self.is_excitable()
        except AssertionError:
            raise StructureError('violated_bonding_laws')
        electron_nr = self.count_electrons()
        electron_ids = []

        for orbital_set in ATOM_PROPERTIES.orbital_order:
            if orbital_set in self.orbital_sets:
                for orbital in self.orbital_sets[orbital_set].orbitals:
                    for i in range(orbital.electron_nr):
                        electron_id = orbital.empty_orbital()
                        electron_ids.append(electron_id)

                    if electron_nr > 0:
                        orbital.fill_orbital(electron_ids.pop())
                        electron_nr -= 1

    def get_lone_pairs(self):
        lone_pairs = []
        for orbital in self.orbitals:
            if orbital.electron_nr == 2:
                if orbital.electrons[0].atom == orbital.electrons[1].atom:
                    lone_pairs.append(orbital.electrons)

        return lone_pairs

    def get_lone_pair_nr(self):
        lone_pair_nr = 0
        for orbital in self.orbitals:
            if orbital.electron_nr == 2:
                if orbital.electrons[0].atom == orbital.electrons[1].atom:
                    lone_pair_nr += 1
        return lone_pair_nr

    def get_lone_electrons(self):
        lone_electrons = 0
        for orbital in self.orbitals:
            if orbital.electron_nr == 1:
                lone_electrons += 1

        return lone_electrons

    def find_bonding_orbitals(self):
        self.bonding_orbitals = []

        for orbital in self.orbitals:
            if orbital.electron_nr == 1:
                self.bonding_orbitals.append(orbital)

    def count_electrons(self):
        electron_nr = 0
        for orbital in self.orbitals:
            electron_nr += orbital.electron_nr

        return electron_nr

    def count_orbitals(self):
        orbital_nr = len(self.orbitals)

        return orbital_nr

    def is_excitable(self):
        electron_nr = self.count_electrons()
        orbital_nr = self.count_orbitals()

        if orbital_nr >= electron_nr:
            return True
        else:
            return False

    def drop_electrons(self):

        lone_orbitals = []

        for orbital_set in ATOM_PROPERTIES.orbital_order:
            if orbital_set in self.orbital_sets:
                for orbital in self.orbital_sets[orbital_set].orbitals:
                    if orbital.electron_nr == 1:
                        if not orbital.electrons[0].aromatic:
                            lone_orbitals.append(orbital)

        while len(lone_orbitals) > 1 and (lone_orbitals[0].orbital_type != lone_orbitals[-1].orbital_type or
                                          lone_orbitals[0].orbital_nr != lone_orbitals[-1].orbital_nr):
            receiver_orbital = lone_orbitals[0]
            donor_orbital = lone_orbitals[-1]

            moved_electron = donor_orbital.electrons[0]

            donor_orbital.remove_electron(moved_electron)
            receiver_orbital.add_electron(moved_electron)
            receiver_orbital.electrons[1].set_orbital(receiver_orbital)

            del lone_orbitals[0]
            del lone_orbitals[-1]

    def print_shell(self):
        for orbital in self.orbitals:
            print(orbital)
            print(orbital.electrons)
            