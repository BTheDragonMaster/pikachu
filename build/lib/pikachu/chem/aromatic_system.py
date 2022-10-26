class AromaticSystem:
    def __init__(self, aromatic_system_id, atoms):
        self.id = aromatic_system_id
        self.atoms = set(atoms)
        for atom in self.atoms:
            atom.aromatic_system = self
        self.electrons = set()
        self.bonds = set()
        self.set_bonds()
        for bond in self.bonds:
            bond.aromatic_system = self

        self.set_electrons()

    def __hash__(self):
        return self.id

    def __eq__(self, other):
        return self.id == other.id

    def __repr__(self):
        return f"Aromatic system {self.id}, members: {', '.join([atom.__repr__() for atom in self.atoms])}"

    def set_bonds(self):
        for atom_1 in self.atoms:
            for atom_2 in self.atoms:
                if atom_1 != atom_2:
                    bond = atom_1.get_bond(atom_2)
                    if bond:
                        self.bonds.add(bond)
                        bond.aromatic_system = self

    def get_contributed_electrons(self, atom):
        contributed_electrons = []
        for electron in self.electrons:
            if electron.atom == atom:
                contributed_electrons.append(electron)

        return contributed_electrons

    def set_electrons(self):
        for atom in self.atoms:

            p_orbital = atom.get_orbitals('p')[0]
            electrons_participate_in_system = True
            for electron in p_orbital.electrons:
                if electron.atom not in self.atoms:
                    electrons_participate_in_system = False

            if electrons_participate_in_system:
                electrons = p_orbital.electrons[:]
                for electron in electrons:
                    p_orbital.remove_electron(electron)
                    self.electrons.add(electron)

    def add_atom(self, atom):
        self.atoms.add(atom)
        atom.aromatic_system = self

    def relocalise_electrons(self):
        for atom in self.atoms:
            p_orbital = atom.get_orbitals('p')[0]
            electrons = self.get_contributed_electrons(atom)
            for electron in electrons:
                p_orbital.add_electron(electron)
                self.remove_electron(electron)

    def add_electron(self, electron):
        self.electrons.add(electron)

    def remove_atom(self, atom):
        self.atoms.remove(atom)
        atom.aromatic_system = None

    def remove_electron(self, electron):
        self.electrons.remove(electron)

    def add_bond(self, bond):
        self.bonds.add(bond)
        bond.aromatic_system = self

    def remove_bond(self, bond):
        self.bonds.remove(bond)
        bond.aromatic_system = None
