class Electron:
    def __init__(self, shell_nr, orbital_type, orbital_nr, spin, atom):
        self.shell_nr = shell_nr
        self.orbital_type = orbital_type
        self.orbital_nr = orbital_nr
        self.spin = spin
        self.atom = atom
        self.paired = False
        self.partner = None
        self.aromatic = False

    def __repr__(self):

        if self.aromatic:
            aromatic_string = '*'
        else:
            aromatic_string = ''
        if self.orbital_nr:
            return f'{self.atom}_{self.shell_nr}{self.orbital_type}{self.orbital_nr}_{self.spin}{aromatic_string}'
        else:
            return f'{self.atom}_{self.shell_nr}{self.orbital_type}_{self.spin}{aromatic_string}'

    def __hash__(self):
        if self.orbital_nr:
            return f'{self.atom}_{self.shell_nr}{self.orbital_type}{self.orbital_nr}_{self.spin}'
        else:
            return f'{self.atom}_{self.shell_nr}{self.orbital_type}_{self.spin}'

    def set_orbital(self, orbital):
        self.orbital_type = orbital.orbital_type
        self.orbital_nr = orbital.orbital_nr
        self.orbital = orbital

    def set_paired(self):
        self.paired = True

    def set_unpaired(self):
        self.paired = False

    def set_aromatic(self):
        self.aromatic = True

    def set_unaromatic(self):
        self.aromatic = False

    def pair(self, electron):
        self.spin = 0.5
        electron.spin = -0.5

        self.set_paired()
        electron.set_paired()
        self.partner = electron
        electron.partner = self

    def unpair(self):
        self.set_unpaired()
        self.partner = None
        self.spin = 0.5
