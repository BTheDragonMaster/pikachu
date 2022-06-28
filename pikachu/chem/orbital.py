from pikachu.chem.electron import Electron


class OrbitalSet:
    def __init__(self, atom, shell_nr, orbital_type):
        self.atom = atom
        self.shell_nr = shell_nr

        self.orbital_type = orbital_type
        self.orbitals = []
        self.define_orbitals()
        self.capacity = len(self.orbitals) * 2

    def __repr__(self):
        return f'{self.shell_nr}{self.orbital_type}'

    def define_orbitals(self):
        if self.orbital_type == 's':
            self.append_s_orbital()
        if self.orbital_type == 'p':
            self.append_p_orbitals()
        if self.orbital_type == 'd':
            self.append_d_orbitals()
        if self.orbital_type == 'f':
            self.append_f_orbitals()

    def append_s_orbital(self):
        self.orbitals.append(Orbital(self.atom, self.shell_nr, 's'))

    def append_p_orbitals(self):
        self.orbitals.append(Orbital(self.atom, self.shell_nr, 'p', 1))
        self.orbitals.append(Orbital(self.atom, self.shell_nr, 'p', 2))
        self.orbitals.append(Orbital(self.atom, self.shell_nr, 'p', 3))

    def append_d_orbitals(self):
        self.orbitals.append(Orbital(self.atom, self.shell_nr, 'd', 1))
        self.orbitals.append(Orbital(self.atom, self.shell_nr, 'd', 2))
        self.orbitals.append(Orbital(self.atom, self.shell_nr, 'd', 3))
        self.orbitals.append(Orbital(self.atom, self.shell_nr, 'd', 4))
        self.orbitals.append(Orbital(self.atom, self.shell_nr, 'd', 5))

    def append_f_orbitals(self):
        self.orbitals.append(Orbital(self.atom, self.shell_nr, 'f', 1))
        self.orbitals.append(Orbital(self.atom, self.shell_nr, 'f', 2))
        self.orbitals.append(Orbital(self.atom, self.shell_nr, 'f', 3))
        self.orbitals.append(Orbital(self.atom, self.shell_nr, 'f', 4))
        self.orbitals.append(Orbital(self.atom, self.shell_nr, 'f', 5))
        self.orbitals.append(Orbital(self.atom, self.shell_nr, 'f', 6))
        self.orbitals.append(Orbital(self.atom, self.shell_nr, 'f', 7))

    def fill_orbitals(self, electrons, electron_nr):
        while electrons > 0:
            for orbital in self.orbitals:
                if electrons > 0:
                    orbital.fill_orbital(electron_nr)
                    electrons -= 1
                    electron_nr += 1
                else:
                    break


class Orbital:
    subtype_dict = {'p': {1: 'x',
                          2: 'y',
                          3: 'z'},
                    'd': {1: 'z^2',
                          2: 'zx',
                          3: 'yz',
                          4: 'xy',
                          5: 'x^2-y^2'},
                    'f': {1: 'z^3-3/5zr^2',
                          2: 'x^3-3/5xr^2',
                          3: 'y^3-3/5yr^2',
                          4: 'xyz',
                          5: 'y(x^2-z^2)',
                          6: 'x(z^2-y^2)',
                          7: 'z(x^2-y^2)'}
                    }

    def __init__(self, atom, shell_nr, orbital_type, orbital_nr=None):
        self.shell_nr = shell_nr
        self.orbital_type = orbital_type
        self.orbital_nr = orbital_nr
        self.electron_nr = 0
        self.electrons = []
        self.atom = atom
        self.bond = None
        self.bonding_orbital = None

    def __hash__(self):
        if self.orbital_nr:

            return f'{self.shell_nr}{self.orbital_type}{self.orbital_nr}'

        else:
            return f'{self.shell_nr}{self.orbital_type}'

    def __repr__(self):
        if self.orbital_nr:

            return f'{self.shell_nr}{self.orbital_type}{self.orbital_nr}'

        else:
            return f'{self.shell_nr}{self.orbital_type}'

    def set_electron_nr(self):
        self.electron_nr = len(self.electrons)

    def set_bond(self, bond, bonding_orbital):
        self.bond = bond
        self.bonding_orbital = bonding_orbital

    def remove_bond(self):
        self.bond = None
        self.bonding_orbital = None

    def fill_orbital(self, electron_id):
        """
        """
        assert self.electron_nr < 2

        self.electrons.append(Electron(electron_id, self.shell_nr, self.orbital_type,
                                       self.orbital_nr, 0.5, self.atom))
        self.set_electron_nr()

        if self.electron_nr == 2:
            self.electrons[0].pair(self.electrons[1])

    def empty_orbital(self):
        """
        """
        assert self.electron_nr > 0

        electron_id = self.electrons[-1].id

        del self.electrons[-1]
        self.set_electron_nr()

        if self.electron_nr == 1:
            self.electrons[0].unpair()
            
        return electron_id

    def add_electron(self, electron):
        assert self.electron_nr < 2

        self.electrons.append(electron)
        self.set_electron_nr()

        if self.electron_nr == 2:
            self.electrons[0].pair(self.electrons[1])

    def remove_electron(self, electron):
        assert electron in self.electrons

        self.electrons.remove(electron)
        self.set_electron_nr()

        if self.electron_nr == 1:
            self.electrons[0].unpair()
