from typing import TYPE_CHECKING, Optional, List, Tuple, Dict

import copy

from pikachu.chem.lone_pair import LonePair
from pikachu.chem.atom_properties import ATOM_PROPERTIES
from pikachu.chem.bond_properties import BOND_PROPERTIES
from pikachu.errors import StructureError
from pikachu.chem.shell import Shell
from pikachu.math_functions import Vector

if TYPE_CHECKING:
    from pikachu.chem.bond import Bond
    from pikachu.chem.aromatic_system import AromaticSystem
    from pikachu.chem.electron import Electron
    from pikachu.chem.structure import Structure
    from pikachu.chem.orbital import OrbitalSet, Orbital


class Atom:
    """
    Class to store an atom

    Attributes:
        type: str, atom type as represented in the periodic table. * for wildcard atom
        nr: int, atom index. Must be unique within a structure
        chiral: Optional[str], 'clockwise' or 'counterclockwise',
            indicates chirality according to order of neighbouring atoms
        charge: int, charge of the atom
        aromatic: bool, True if atom is part of an aromatic system, False if otherwise
        shells: dict of {shell_nr: Shell, ->}
    """
    #
    # def __new__(cls, atom_type: str, atom_nr: int, chiral: Optional[str], charge: int, aromatic: bool):
    #     self = super().__new__(cls)  # Must explicitly create the new object
    #     # Aside from explicit construction and return, rest of __new__
    #     # is same as __init__
    #     self.type = atom_type
    #     self.nr = atom_nr
    #     self.chiral = chiral
    #     self.charge = charge
    #     self.aromatic = aromatic
    #     self.shells = {}
    #
    #     return self  # __new__ returns the new object
    #
    # def __getnewargs__(self):
    #     # Return the arguments that *must* be passed to __new__
    #     return self.type, self.nr, self.chiral, self.charge, self.aromatic

    def __init__(self, atom_type: str, atom_nr: int, chiral: Optional[str], charge: int, aromatic: bool) -> None:

        self.type: str = atom_type
        self.nr: int = atom_nr
        self.chiral: Optional[str] = chiral
        self.charge: int = charge
        self.aromatic: bool = aromatic

        self.pyrrole: bool = False
        self.furan: bool = False
        self.thiophene: bool = False

        self.shells: Dict[int, "Shell"] = {}
        self.lone_pairs: List[LonePair] = []
        self.draw = AtomDrawProperties()
        self.annotations = AtomAnnotations()
        self.hybridisation: str = ''

        self.bonds: List["Bond"] = []
        self.neighbours: List["Atom"] = []
        self.drawn_neighbours: List["Atom"] = []
        self.aromatic_system: Optional["AromaticSystem"] = None
        self.connectivity: Optional[Tuple[str]] = None

        self.shell_nr: int = ATOM_PROPERTIES.element_to_shell_nr[self.type]

    def __eq__(self, atom: "Atom") -> bool:
        if type(atom) == Atom:
            return self.nr == atom.nr
        else:
            return False

    def __hash__(self) -> int:
        return self.nr

    def __repr__(self) -> str:
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

    def copy(self) -> "Atom":
        """
        Creates and returns a copy of the current atom

        """
        atom_copy = Atom(self.type, self.nr, self.chiral, self.charge, self.aromatic)
        atom_copy.hybridisation = self.hybridisation
        atom_copy.pyrrole = self.pyrrole
        atom_copy.furan = self.furan
        atom_copy.thiophene = self.thiophene
        # TODO: Make proper copy object for AtomDrawProperties
        atom_copy.draw = AtomDrawProperties()
        atom_copy.draw.colour = self.draw.colour
        atom_copy.annotations = self.annotations.copy()
        if self.connectivity is None:
            atom_copy.connectivity = None
        else:
            connectivity: List[str] = []

            for connection in self.connectivity:
                connectivity.append(connection)

            atom_copy.connectivity = tuple(connectivity)

        for neighbour in self.neighbours:
            atom_copy.neighbours.append(neighbour)

        # Note: when copying a structure, neighbours need to be refreshed again
        # after all atoms have been changed!

        return atom_copy

    def _set_lone_pairs(self) -> None:
        """
        Find and store lone pairs associated with the atom
        """
        lone_pairs: List[Tuple[Electron]] = self.valence_shell.get_lone_pairs()

        for i, electrons in enumerate(lone_pairs):
            lone_pair = LonePair(self, self.nr * (i + 10000))
            for electron in electrons:
                lone_pair.add_electron(electron)

            self.lone_pairs.append(lone_pair)

    def get_bond(self, atom: "Atom") -> Optional["Bond"]:
        """
        Return None if there exists no bond between atom and self; Bond instance if a bond does exist

        Parameters
        ----------
        atom: Atom instance

        """
        for bond in self.bonds:
            if bond.atom_1 == atom or bond.atom_2 == atom:
                return bond

        return None

    def get_ring_index(self, structure: "Structure") -> Optional[int]:
        cycles = structure.cycles.all_cycles

        for i, cycle in enumerate(cycles):
            if self in cycle:
                return i

        return None

    def get_ring(self, structure: "Structure") -> Optional[List["Atom"]]:
        cycles = structure.sssr

        for i, cycle in enumerate(cycles):
            if self in cycle:
                return cycle

        return None

    def is_inside_ring(self, structure: "Structure") -> bool:
        cycles = structure.cycles.all_cycles

        for i, cycle in enumerate(cycles):
            if self in cycle:
                return True

        return False

    def _set_neighbours(self, structure: "Structure") -> None:
        """
        Set atom neighbours from a list of atoms

        Parameters
        ----------
        structure: Structure instance

        """
        self.neighbours = structure.graph[self]

    # TODO: Move method to AtomDrawProperties
    def set_drawn_neighbours(self) -> None:
        """
        Store a list of neighbours that will be drawn in visualisation
        """
        self.drawn_neighbours = []
        for neighbour in self.neighbours:
            if neighbour.draw.is_drawn:
                self.drawn_neighbours.append(neighbour)

    def _remove_neighbour(self, neighbour: "Atom") -> None:
        """
        Remove a neighbouring atom.
        """
        self.neighbours.remove(neighbour)

    # TODO: Move method to AtomDrawProperties
    def get_drawn_neighbours(self) -> List["Atom"]:
        """
        Returns all neighbours of the atom that are explicitly drawn in the visualisation
        """
        drawn_neighbours = []
        for neighbour in self.neighbours:
            if neighbour.draw.is_drawn:
                drawn_neighbours.append(neighbour)

        return drawn_neighbours

    def has_neighbour(self, atom_type: str) -> bool:
        """
        Return True if atom has a neighbour of a certain atom type, False if not

        Parameters
        ----------
        atom_type: str, atom type
        """
        for neighbour in self.neighbours:
            if neighbour.type == atom_type:
                return True

        return False

    def set_connectivity(self) -> None:
        """
        Sets the connectivity of an atom (sorted tuple of str, with each string a concatenation of atom type and
            bond type)

        Example: ('C_single', 'O_double')
        """
        self.connectivity = self.get_connectivity()

    def get_connectivity(self) -> Tuple[str]:
        """
        Returns a sorted tuple of str, representing atom connectivity

        Example: ('C_single', 'O_double')
        """
        connectivity = []

        for bond in self.bonds:
            for atom in bond.neighbours:
                if atom.type != 'H' and atom != self:
                    connectivity.append(f'{atom.type}_{bond.type}')

        connectivity = tuple(sorted(connectivity))
        return connectivity

    def _has_similar_connectivity(self, substructure_connectivity: Tuple[str]) -> bool:
        """
        Return True if all bond-atom combinations of the child atom are present in the parent atom, False otherwise

        Parameters
        ----------
        substructure_connectivity: tuple of str, with each str representing the connectivity of an atom
            Example: ('C_single', 'O_double')

        """
        parent_connectivity_copy = list(copy.copy(self.connectivity))
        substructure_connectivity_copy = list(copy.copy(substructure_connectivity))

        same_connectivity = True

        for atom in substructure_connectivity_copy:
            if atom in parent_connectivity_copy:
                parent_connectivity_copy.remove(atom)
            else:
                same_connectivity = False

        return same_connectivity

    def _add_electron_shells(self) -> None:
        """
        Fill all shells with electrons and excite electrons in valence shell
        """

        # Generate empty shells depending on atom type
        self.__make_shells()

        nr_double_bonds: int = 0
        nr_single_bonds: int = 0

        for bond in self.bonds:
            if bond.type == 'double':
                nr_double_bonds += 1
            elif bond.type == 'single':
                nr_single_bonds += 1

        # Deals with nitro-groups defined incorrectly in SMILES strings
        if self.type == 'N' and self.charge == 0 and nr_double_bonds == 2 and nr_single_bonds == 1:
            oxygen_bonds: List["Bond"] = []
            oxygens: List["Atom"] = []

            for bond in self.bonds:
                neighbour = bond.get_connected_atom(self)
                if bond.type == 'double' and neighbour.type == 'O':
                    oxygens.append(neighbour)
                    oxygen_bonds.append(bond)

            if len(oxygens) >= 1:
                oxygen = oxygens[0]
                bond = oxygen_bonds[0]
                bond.type = 'single'
                bond.set_bond_summary()
                oxygen.charge = -1
                self.charge = 1

                if oxygen.shells:
                    oxygen._add_electron_shells()

        # Fill them with electrons
        self.__fill_shells()

        # Check if the total number of electrons can be distributed such that each orbital contains one electron each
        is_excitable = self.valence_shell.is_excitable()

        # Can't excite carbon if it has more than 4 electrons
        if (self.type == 'C' or self.type == 'B') and is_excitable:
            self._excite()
        else:
            # Assigns a value to all bonds: 2 for double bond, 1 for single. The aromatic value depends on the atom
            bond_weights: List[int] = []
            aromatic_count = 0

            for bond in self.bonds:
                if bond.type == 'aromatic':
                    aromatic_count += 1

                bond_weights.append(BOND_PROPERTIES.bond_type_to_weight[bond.type])

            # If odd number of aromatic bonds (such as central atoms in Trp), only add 1 'extra' bond for the
            # three outgoing aromatic bonds

            nr_h_bonds = 0
            for bond in self.bonds:
                if bond.get_connected_atom(self).type == 'H':
                    nr_h_bonds += 1

            nr_non_h_bonds = sum(bond_weights) + int(aromatic_count / 2)

            if self.pyrrole or self.furan or self.thiophene or self.__is_aromatic_nitrogen():
                nr_non_h_bonds -= 1

            # TODO: Does this work for all atoms? Doesn't for carbon. Should this be made general?

            bonding_electrons = self.__get_bonding_electrons()

            if nr_non_h_bonds > bonding_electrons:
                if is_excitable:
                    self._excite()

                elif nr_h_bonds:

                    nr_non_h_bonds -= nr_h_bonds
                    if nr_non_h_bonds > bonding_electrons:
                        if is_excitable:
                            self._excite()
                        else:
                            raise StructureError('violated_bonding_laws')
                else:
                    raise StructureError('violated_bonding_laws')

    def __get_bonding_electrons(self) -> int:
        """
        Returns the number of unpaired electrons in the valence shell
        """
        nr_unpaired_electrons = 0
        for orbital in self.valence_shell.orbitals:
            if len(orbital.electrons) == 1:
                nr_unpaired_electrons += 1
        return nr_unpaired_electrons

    def __make_shells(self) -> None:
        """
        Create electron shells for an atom
        """
        for i in range(self.shell_nr):
            current_shell: int = i + 1
            self.shells[current_shell]: "Shell" = Shell(self, current_shell)

        self.valence_shell: "Shell" = self.shells[self.shell_nr]

    def in_ring(self, structure: "Structure") -> bool:
        """
        Returns True if atom is part of a ring (system), False if otherwise

        Parameters
        ----------
        structure: Structure instance. Must contain the atom

        """
        assert self in structure.graph
        cycles: List[List["Atom"]] = structure.sssr

        for cycle in cycles:
            if self in cycle:
                return True

        return False

    def _adjacent_to_stereobond(self) -> bool:
        """
        Returns True if atom is adjacent to a stereochemically restrictive bond, False otherwise
        """
        for bond in self.bonds:
            if bond.chiral:
                return True

        return False

    # TODO: raise exception if there are more electrons than fit in the orbitals
    def __fill_shells(self) -> None:
        """
        Fill electron shells based on atom type
        """

        electrons_remaining: int = ATOM_PROPERTIES.element_to_atomic_nr[self.type] - self.charge
        electron_nr: int = 1

        # Iterate over the orbitals in order of them being filled

        for orbital in ATOM_PROPERTIES.orbital_order:
            if electrons_remaining > 0:
                # Find the correct set of orbitals
                shell: int = int(orbital[0])
                orbital_set: "OrbitalSet" = self.shells[shell].orbital_sets[orbital]

                # Place either the maximum number of electrons or the number of electrons remaining in the orbital set

                electrons_to_dump: int = min([electrons_remaining, orbital_set.capacity])
                orbital_set.fill_orbitals(electrons_to_dump, electron_nr)
                electron_nr += electrons_to_dump

                electrons_remaining -= electrons_to_dump
            else:
                # All electrons have been placed and we can break out of the loop
                break
                
    def get_neighbour(self, atom_type: str) -> Optional["Atom"]:
        """
        Returns the first encountered neighbour of type atom_type if such a neighbour exists, None otherwise

        Parameters
        ----------
        atom_type: str, atom type according to periodic table

        """
        for neighbour in self.neighbours:
            if neighbour.type == atom_type:
                return neighbour
        return None

    def get_neighbours(self, atom_type: str) -> List["Atom"]:
        """
        Returns a list of all neighbours of type atom_type

        Parameters
        ----------
        atom_type: str, atom type according to periodic table

        """
        neighbours: List["Atom"] = []
        for neighbour in self.neighbours:
            if neighbour.type == atom_type:
                neighbours.append(neighbour)
        return neighbours

    def _excite(self) -> None:
        """
        Excite electrons in valence shell
        """
        self.valence_shell.excite()

    def get_non_hydrogen_neighbours(self) -> List["Atom"]:
        """
        Returns a list of neighbouring, non-hydrogen atoms
        """
        neighbours: List["Atom"] = []
        for atom in self.neighbours:
            if atom.type != 'H' and atom.type != '*':
                neighbours.append(atom)
        return neighbours

    def get_non_hydrogen_bonds(self) -> List["Bond"]:
        """
        Returns a list of neighbouring bonds that do not neighbour a hydrogen atom
        """
        bonds: List["Bond"] = []
        for bond in self.bonds:
            if bond.atom_1.type != 'H' and bond.atom_2 == self:
                bonds.append(bond)
            elif bond.atom_2.type != 'H' and bond.atom_1 == self:
                bonds.append(bond)
        return bonds

    def _remove_bond(self, bond: "Bond") -> None:
        """
        Remove a bond from a list of neighbouring bonds. Does not break the bond!

        To break a bond, use the 'break_bond' method in the Structure instance instead.

        Parameters
        ----------
        bond: Bond instance, must be a bond neighbouring the atom

        """
        assert bond in self.bonds
        self.bonds.remove(bond)

    # TODO: check if the steric number can be determined more cleanly - after all, electrons have already
    # been dropped at this stage

    def __get_nr_electron_pairs(self) -> int:
        """
        Returns the number of electron pairs that currently sit in excited orbitals

        Returns
        -------
        electron_pair_nr: int, number of electron pairs that can be formed from lone electrons currently in
            excited orbitals.
        """

        valence: int = self.get_valence()
        bonds_accounted_for: int = 0
        electron_nr: int = 0

        for orbital in self.valence_shell.orbitals:

            # Count number of unpaired electrons

            if orbital.electron_nr == 1:
                electron_nr += 1

            # Count number of paired electrons
            elif orbital.electron_nr == 2:
                # Lone pairs
                if orbital.electrons[0].atom == orbital.electrons[1].atom:
                    electron_nr += 2
                # Electrons involved in bonds
                else:
                    bonds_accounted_for += 1

        nr_bonds_to_make: int = valence - bonds_accounted_for

        nr_unbonded_electrons: int = electron_nr - nr_bonds_to_make

        # TODO: Implement radicals; check why this code gets flagged when refreshing an aromatic structure

        if nr_unbonded_electrons % 2 != 0:
            pass
            # print("Warning! Rogue electron.")
            # print(self)
            # print(self)
            # print(bond_nr)
            # print(self.bonds)
            # print(bonds_to_make)
            # print(bonds_accounted_for)
            # print(electron_nr)
            # self.valence_shell.print_shell()

        electron_pair_nr: int = int(nr_unbonded_electrons / 2)

        return electron_pair_nr

    def _drop_electrons(self) -> None:
        """
        Drop excited electrons that did not end up forming pairs with electrons from bonded atoms to lower orbitals
        """
        if self.valence_shell.get_lone_electrons() > 1:
            self.valence_shell.drop_electrons()

    def get_valence(self) -> int:
        """
        Returns the valence of the atom based on its neighbours and attached bonds
        """

        valence: int = 0
        aromatic_bond_nr: int = 0

        for bond in self.bonds:
            if bond.type == 'single':
                valence += 1
            elif bond.type == 'double':
                valence += 2
            elif bond.type == 'triple':
                valence += 3
            elif bond.type == 'quadruple':
                valence += 4
            elif bond.type == 'aromatic':
                aromatic_bond_nr += 1

        if aromatic_bond_nr == 2:
            # Pyrrole, furan and thiophene systems have a valence of 2 + the number of non-aromatic bonds
            if self.pyrrole or self.furan or self.thiophene or self.__is_aromatic_nitrogen():

                valence += 2
            elif self.aromatic:
                oxygen: Optional["Atom"] = None
                for bond in self.bonds:
                    connected_atom: "Atom" = bond.get_connected_atom(self)

                    if bond.type == 'double' and connected_atom.type == 'O':
                        oxygen = connected_atom

                # If a double-bonded oxygen can provide resonance, the valence of the atom is 2 +
                # the number of non-aromatic bonds
                if oxygen is not None and oxygen._resonance_possible(self):
                    valence += 2
                # Otherwise, the valence of the atom is 3 + the number of non-aromatic bonds
                else:
                    valence += 3

            # If the atom was not labelled as aromatic in the SMILES, it might still be aromatic.
            # In this case, label as a valence of 3 + the number of non-aromatic bonds
            else:
                valence += 3
        elif aromatic_bond_nr == 3 and self.type == 'C':
            valence += 4
        # TODO: Check if this should be generalised for other atom types
        elif aromatic_bond_nr == 3 and self.type == 'N':
            if self.charge == 1:
                valence += 4
            else:
                valence += 3

        return valence

    def _is_promotable(self) -> bool:
        """
        Returns True if electrons in the valence shell can be promoted to a d-orbital, False if otherwise
        """
        promotable: bool = False
        for orbital_set in self.valence_shell.orbital_sets:
            if 'd' in orbital_set:
                promotable = True

        return promotable

    def __is_aromatic_nitrogen(self) -> bool:
        """
        Returns True if the current nitrogen atom is marked as aromatic in the input, False otherwise
        """
        if self.type == 'N' and len(self.bonds) == 3 and self.aromatic and self.charge == 0:
            return True
        return False

    def _resonance_possible(self, neighbour: "Atom") -> bool:
        """
        Return True if the current atom can participate in resonance with a neighbouring aromatic system,
            False otherwise

        Parameters
        ----------
        neighbour: Atom instance
        """
        if self.type == 'O' and len(self.bonds) == 1 and self.bonds[0].type == 'double' and neighbour.aromatic:
            return True
        return False

    def _promote_lone_pair_to_p_orbital(self) -> None:
        """
        Promote the electrons of a single lone pair in an aromatic system to a p-orbital,
            such that they can be delocalised appropriately to participate in the system
        """
        # Promotion of lone pairs to p-orbitals can only be done starting from a sp3-hybridised system
        assert self.hybridisation == 'sp3'

        # Remove hybridisation
        self.valence_shell.dehybridise()

        p_orbitals: List["Orbital"] = []
        sp2_orbitals: List["Orbital"] = []

        for orbital in self.valence_shell.orbitals:
            if orbital.electron_nr == 2:
                # Any electrons that are already involved in bonds will localise to sp2 orbitals
                # Example: the bond between an aromatic carbon and a pyrrole nitrogen
                if orbital.electrons[0].atom != orbital.electrons[1].atom and \
                        (orbital.orbital_type == 's' or orbital.orbital_type == 'p'):
                    sp2_orbitals.append(orbital)
                # Any electrons that are not involved in bonds can be delocalised to p orbitals
                # Example: one of the lone pairs of a furan oxygen
                elif orbital.electrons[0].atom == orbital.electrons[1].atom == self and \
                        (orbital.orbital_type == 's' or orbital.orbital_type == 'p'):
                    p_orbitals.append(orbital)
            else:
                if orbital.orbital_type == 's' or orbital.orbital_type == 'p':
                    sp2_orbitals.append(orbital)

        assert len(p_orbitals) >= 1

        # Should more than one p-orbital be found, make sure we only delocalise the electrons of 1 to a p-orbital,
        # and the rest to sp2 orbitals

        if len(p_orbitals) > 1:
            for i in range(0, len(p_orbitals) - 1):
                sp2_orbitals.append(p_orbitals[i])

        p_orbital = p_orbitals[-1]

        # Change one selected sp3 orbital to a p orbital

        p_orbital.orbital_type = 'p'
        p_orbital.orbital_nr = 1

        # Change the other sp3 orbitals to sp2 orbitals

        for i, orbital in enumerate(sp2_orbitals):
            orbital.orbital_type = 'sp2'
            orbital.orbital_nr = i + 1

        # Change the hybridisation to sp2

        self.hybridisation = 'sp2'

        # Store the new orbitals in the electrons of the valence shell

        for orbital in self.valence_shell.orbitals:
            for electron in orbital.electrons:
                if electron.atom == self:
                    electron.set_orbital(orbital)

    # TODO: Change orbital type to Enum
    def _get_orbitals(self, orbital_type: str) -> List["Orbital"]:
        """
        Return all orbitals of a certain orbital type

        Parameters
        ----------
        orbital_type: str, orbital type
        """
        orbitals: List["Orbital"] = []
        for orbital in self.valence_shell.orbitals:
            if orbital.orbital_type == orbital_type:
                orbitals.append(orbital)

        return orbitals

    def _get_hybrid_orbitals(self, orbital_type: str) -> List["Orbital"]:
        """
        Return all orbitals that are a hybrid orbital with the given orbital type

        Parameters
        ----------
        orbital_type: str, orbital type or subtype
            e.g. given orbital type 'p', function will return any orbitals containing 'p', such as 'sp' and 'sp2'
        """
        orbitals: List["Orbital"] = []
        for orbital in self.valence_shell.orbitals:
            if orbital_type in orbital.orbital_type:
                orbitals.append(orbital)

        return orbitals

    # TODO: Check when this function is needed
    def _promote_pi_bonds_to_d_orbitals(self) -> None:
        """
        Promote pi-bonds from sp[x]d[y] orbitals to d-orbitals,
            e.g. SF4, which contains 5 electron pairs and 4 bonded atoms
        """

        if self._is_promotable() and 'd' in self.hybridisation:
            
            donor_orbitals: List["Orbital"] = []
            receiver_orbitals: List["Orbital"] = []

            for orbital in self.valence_shell.orbitals:
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
                donor_orbital._remove_bond()

    # TODO: Check when this function is called
    def _promote_pi_bond_to_d_orbital(self) -> None:
        """
        Promote pi-bonds from sp[x]d[y] orbitals to d-orbitals
        """
        assert self._is_promotable()

        donor_orbitals = []
        receiver_orbitals = []
        for orbital in self.valence_shell.orbitals:
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
        donor_orbital._remove_bond()

        self.valence_shell.dehybridise()

        self.hybridise()

    def _reset_hybridisation(self) -> None:
        """
        Reset the hybridisation of the atom. Called after updating an atom's direct neighbourhood
        """
        self.valence_shell.dehybridise()
        self.hybridise()

    def _get_nr_implicit_hydrogens(self) -> int:
        """
        Returns the number of implicit hydrogens for the atom given its atom type and bonds to non-H atoms

        """
        nr_hydrogens: int = 0

        # These are the only atom types for which hydrogens can be implicit in SMILES
        if self.type in ['B', 'C', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I']:

            valence = self.get_valence()
            if valence in ATOM_PROPERTIES.element_to_valences[self.type]:
                nr_hydrogens = 0
            else:
                max_bonds = self.valence_shell.get_lone_electrons()
                nr_hydrogens = max_bonds - valence

        return nr_hydrogens

    def _add_bond(self, bond: "Bond") -> None:
        """
        Add neighbouring bond to the atom. Does not rearrange electrons or make a bond.

        To make a new bond, use the 'make_bond' method of class Structure instead

        Parameters
        ----------
        bond: Bond instance
        """
        self.bonds.append(bond)

    def hybridise(self) -> None:
        """
        Hybridise the atom and valence shell
        """

        hybridisation = self.__get_hybridisation()
        self.valence_shell.hybridise(hybridisation)
        self.__set_hybridisation()

    def __set_hybridisation(self) -> None:
        """
        Set the hybridisation for the atom
        """

        self.hybridisation = 's'
        for orbital in self.valence_shell.orbitals:
            if orbital.orbital_type in {'sp', 'sp2', 'sp3', 'sp3d', 'sp3d2'}:
                self.hybridisation = orbital.orbital_type
                break

    # TODO: Use enums for hybridisation
    def __get_hybridisation(self) -> str:
        """
        Returns the hybridisation of the atom based on steric number
        """
        steric_number = self.__get_nr_electron_pairs() + len(self.bonds)
        hybridisation = None
        if steric_number in ATOM_PROPERTIES.steric_nr_to_hybridisation:
            hybridisation = ATOM_PROPERTIES.steric_nr_to_hybridisation[steric_number]

        return hybridisation


class AtomDrawProperties:
    def __init__(self, x=0, y=0):
        self.rings: List[int] = []
        self.original_rings = []
        self.anchored_rings = []
        self.is_bridge_atom: bool = False
        self.is_bridge: bool = False
        self.bridged_ring = None
        self.is_drawn: bool = True
        self.has_hydrogen: bool = False
        self.positioned: bool = False
        self.previous_position: "Vector" = Vector(0, 0)
        self.position: "Vector" = Vector(x, y)
        self.angle: Optional[float] = None
        self.force_positioned: bool = False
        self.connected_to_ring: bool = False
        self.draw_explicit: bool = False
        self.drawn_neighbours: List["Atom"] = []
        self.previous_atom: Optional["Atom"] = None
        self.colour: str = 'black'

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


class AtomAnnotations:
    def __init__(self):

        self.annotations = set()

    def copy(self):
        annotation_copy = AtomAnnotations()
        for annotation in self.annotations:
            annotation_copy.add_annotation(annotation, self.get_annotation(annotation))

        return annotation_copy

    def add_annotation(self, name, default):

        assert getattr(self, name, 'zahar') == 'zahar'
        setattr(self, name, default)
        self.annotations.add(name)

    def set_annotation(self, name, value):
        getattr(self, name)
        setattr(self, name, value)

    def get_annotation(self, name):
        return getattr(self, name)

    def has_annotation(self, name):
        if getattr(self, name, None) is not None:
            return True
        return False

    def print_annotations(self):
        for annotation in self.annotations:
            print(annotation)
