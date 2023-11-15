from typing import List

from pikachu.chem.atom import Atom
from pikachu.chem.electron import Electron


class AromaticSystem:
    """
    Class for storing aromatic systems, including bonds, atoms, and involved electrons

    Attributes:
        id: int, unique identifier of aromatic system
        atoms: set, atoms involved in aromatic system
        electrons: set, electrons shared in aromatic system
        bonds: set, bonds involved in aromatic system
    """
    def __init__(self, aromatic_system_id: int, atoms: List[Atom]):
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

    def __eq__(self, other: "AromaticSystem"):
        return self.id == other.id

    def __repr__(self):
        return f"Aromatic system {self.id}, members: {', '.join([atom.__repr__() for atom in self.atoms])}"

    def set_bonds(self):
        """
        Find and store bonds that are part of the aromatic system from involved atoms
        """
        for atom_1 in self.atoms:
            for atom_2 in self.atoms:
                if atom_1 != atom_2:
                    bond = atom_1.get_bond(atom_2)
                    if bond:
                        self.bonds.add(bond)
                        bond.aromatic_system = self

    def get_contributed_electrons(self, atom: Atom):
        """
        Returns electrons that an atom contributes to an aromatic system

        Parameters
        ----------
        atom: Atom instance

        """
        contributed_electrons = []
        for electron in self.electrons:
            if electron.atom == atom:
                contributed_electrons.append(electron)

        return contributed_electrons

    def set_electrons(self):
        """
        Move electrons from aromatic atoms from their p-orbital to the aromatic system
        """
        for atom in self.atoms:

            p_orbital = atom._get_orbitals('p')[0]
            electrons_participate_in_system = True
            for electron in p_orbital.electrons:
                if electron.atom not in self.atoms:
                    electrons_participate_in_system = False

            if electrons_participate_in_system:
                electrons = p_orbital.electrons[:]
                for electron in electrons:
                    p_orbital.remove_electron(electron)
                    self.add_electron(electron)

    def relocalise_electrons(self):
        """
        Move electrons from the aromatic system back to their p-orbitals - done for kekulisation purposes

        """
        for atom in self.atoms:
            p_orbital = atom._get_orbitals('p')[0]
            electrons = self.get_contributed_electrons(atom)
            for electron in electrons:
                p_orbital.add_electron(electron)
                self.remove_electron(electron)

    def add_electron(self, electron: Electron):
        """
        Add electron to aromatic system

        Parameters
        ----------
        electron: Electron instance

        """
        self.electrons.add(electron)

    def remove_electron(self, electron: Electron):
        """
        Remove electron from aromatic system

        Parameters
        ----------
        electron: Electron instance
        """
        self.electrons.remove(electron)
