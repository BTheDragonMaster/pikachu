import math
from pikachu.errors import StructureError


def get_permissible_double_bond_number(aromatic_bonds):
    non_doubleable = set()
    for aromatic_bond in aromatic_bonds:
        for bond_1 in aromatic_bond.atom_1.bonds:
            if bond_1 not in aromatic_bonds and bond_1.type == 'double':
                non_doubleable.add(aromatic_bond)
        for bond_2 in aromatic_bond.atom_2.bonds:
            if bond_2 not in aromatic_bonds and bond_2.type == 'double':
                non_doubleable.add(aromatic_bond)

        if aromatic_bond.atom_1.hybridisation == 'sp3':
            non_doubleable.add(aromatic_bond)

        if aromatic_bond.atom_2.hybridisation == 'sp3':
            non_doubleable.add(aromatic_bond)

    permissible = list(set(aromatic_bonds) - non_doubleable)
    aromatic_stretches = get_neighbouring_bonds(permissible)

    nr_permissible = sum([int(math.ceil(len(aromatic_stretch) / 2)) for aromatic_stretch in aromatic_stretches])

    return nr_permissible


def is_aromatic(atom_set):
    double_bonds = set()
    aromatic_bonds = set()
    sp3 = []

    for atom_1 in atom_set:
        if atom_1.hybridisation == 'sp3':
            has_lone_pair = False

            for orbital in atom_1.valence_shell.orbitals:
                if len(orbital.electrons) == 2:
                    has_lone_pair = True
                    break
            if has_lone_pair:
                sp3.append(atom_1)
            else:
                return False

        for bond in atom_1.bonds:
            connected_atom = bond.get_connected_atom(atom_1)

            if bond.type == 'double' and connected_atom not in atom_set and connected_atom.type not in {'O', 'S'} and connected_atom.charge <= 0:
                return False

        for atom_2 in atom_set:
            if atom_1 != atom_2:
                bond = atom_1.get_bond(atom_2)
                if bond:
                    if bond.type == 'double':
                        double_bonds.add(bond)
                    elif bond.type == 'aromatic':
                        aromatic_bonds.add(bond)

    if len(aromatic_bonds) == len(atom_set):
        permissible_nr = get_permissible_double_bond_number(aromatic_bonds)

        pi_electrons = (permissible_nr + len(sp3)) * 2

        # Make exception for heme
        if pi_electrons == 16:
            pyrroles = []
            nitrogens = []
            for atom in atom_set:
                if atom.pyrrole:
                    pyrroles.append(atom)
                elif atom.type == 'N':
                    nitrogens.append(atom)
            if len(pyrroles) == 2 and len(nitrogens) == 2:
                pi_electrons = 18
        
        if pi_electrons % 4 != 2:
            raise StructureError('aromaticity')
    elif not aromatic_bonds:
        pi_electrons = (len(sp3) + len(double_bonds)) * 2
    else:

        # This happens when parts of the aromatic system are not defined
        neighbouring_aromatic_bonds = get_neighbouring_bonds(aromatic_bonds)
        if len(neighbouring_aromatic_bonds) == 2:

            if len(neighbouring_aromatic_bonds[0]) == 1 and len(neighbouring_aromatic_bonds[1]) == 1:
                aromatic_bond_1 = neighbouring_aromatic_bonds[0][0]
                aromatic_bond_2 = neighbouring_aromatic_bonds[1][0]

                if aromatic_bond_1.atom_1.hybridisation == aromatic_bond_1.atom_2.hybridisation == aromatic_bond_2.atom_1.hybridisation == aromatic_bond_2.atom_2.hybridisation == 'sp2':
                    pi_electrons = (2 + len(sp3) + len(double_bonds)) * 2

                else:
                    return False

            else:
                return False
        elif len(neighbouring_aromatic_bonds) == 1:
            inaccessible_aromatic_bonds = []
            for aromatic_bond in neighbouring_aromatic_bonds[0]:
                if aromatic_bond.atom_1.hybridisation != 'sp2' or aromatic_bond.atom_2.hybridisation != 'sp2':
                    inaccessible_aromatic_bonds.append(aromatic_bond)

            bond_nr = len(neighbouring_aromatic_bonds[0]) - len(inaccessible_aromatic_bonds)

            pi_electrons = (int(math.ceil(bond_nr / 2)) + len(double_bonds) + len(sp3)) * 2

        else:
            return False

    if pi_electrons % 4 == 2:
        if pi_electrons > 2:
            return True

        if len(atom_set) == 3:
            return True

    else:
        return False


def get_neighbouring_bonds(bonds):
    bond_groups = []
    for bond in bonds:
        bond_groups.append([bond])

    bond_group_nr = len(bond_groups)
    previous_bond_group_nr = -1

    while bond_group_nr != previous_bond_group_nr:
        previous_bond_group_nr = bond_group_nr
        indices_to_remove = None
        new_bond_group = None

        for i, bond_group_1 in enumerate(bond_groups):
            bond_group_1_found = False
            for j, bond_group_2 in enumerate(bond_groups):
                bond_group_2_found = False
                if i != j:
                    for bond_1 in bond_group_1:
                        bond_1_found = False
                        for bond_2 in bond_group_2:
                            if bond_1 != bond_2 and bond_1.bond_is_neighbour(bond_2):
                                new_bond_group = bond_group_1 + bond_group_2
                                indices_to_remove = [i, j]
                                bond_group_1_found = True
                                bond_group_2_found = True
                                bond_1_found = True
                                break
                        if bond_1_found:
                            break
                if bond_group_2_found:
                    break
            if bond_group_1_found:
                break

        if new_bond_group:
            indices_to_remove.sort(reverse=True)
            for index in indices_to_remove:
                bond_groups.pop(index)
            bond_groups.append(new_bond_group)

        bond_group_nr = len(bond_groups)

    return bond_groups