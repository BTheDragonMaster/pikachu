class AromaticSystem:
    def __init__(self, aromatic_system_id):
        self.id = aromatic_system_id
        self.atoms = set()
        self.electrons = set()
        self.bonds = set()

    def add_atom(self, atom):
        self.atoms.add(atom)
        atom.aromatic_system = self

    def add_electron(self, electron):
        self.electrons.add(electron)

    def remove_atom(self, atom):
        self.atoms.remove(atom)
        atom.aromatic_system = None

    def remove_electron(self, electron):
        self.electrons.remove(electron)

    def add_bond(self, bond):
        self.bonds.add(bond)

    def remove_bond(self, bond):
        self.bonds.remove(bond)
