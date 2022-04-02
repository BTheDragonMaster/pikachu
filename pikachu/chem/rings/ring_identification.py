import math


def check_aromatic(atom_set):
    aromatic = True
    for atom in atom_set:
        if atom.hybridisation == 'sp2':
            pass
        else:
            aromatic = False
            break

    if aromatic:
        pi_electron_nr = 0
        for atom in atom_set:
            for orbital in atom.valence_shell.orbitals.values():
                if orbital.orbital_type == 'p':
                    for electron in orbital.electrons:
                        if electron.atom == atom:
                            pi_electron_nr += 1
                        # This should ensure that pi-electrons are only counted if they are part of the ring system
                            
                        elif electron.atom not in atom_set:
                            aromatic = False
                        

        if not pi_electron_nr % 4 == 2:
            aromatic = False

    return aromatic

def electrons_can_delocalise(resonance_atoms, sp3):
    matches = []
    used_resonance_atoms = set()
    if len(resonance_atoms) != len(sp3):
        return False

    for sp3_atom in sp3:
        for resonance_atom in resonance_atoms:
            if not resonance_atom in used_resonance_atoms and \
                    set(sp3_atom.neighbours).intersection(sp3_atom.neighbours):
                matches.append((sp3_atom, resonance_atom))
                used_resonance_atoms.add(resonance_atom)
                break

    if len(matches) == len(sp3) == len(resonance_atoms):
        return True
    else:
        return False


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
            if len(atom_1.lone_pairs) > 0:
                sp3.append(atom_1)
            else:
                return False

        for atom_2 in atom_set:
            if atom_1 != atom_2:
                bond = atom_1.get_bond(atom_2)
                if bond:
                    if bond.type == 'double':
                        double_bonds.add(bond)
                    elif bond.type == 'aromatic':
                        aromatic_bonds.add(bond)

    pi_electrons = None

    if len(aromatic_bonds) == len(atom_set):
        permissible_nr = get_permissible_double_bond_number(aromatic_bonds)
        pi_electrons = (permissible_nr + len(sp3)) * 2
    elif not aromatic_bonds:
        pi_electrons = (len(sp3) + len(double_bonds)) * 2
    else:
        # This happens when parts of the aromatic system are not defined
        neighbouring_aromatic_bonds = get_neighbouring_bonds(aromatic_bonds)
        if len(neighbouring_aromatic_bonds) == 2:
            if len(neighbouring_aromatic_bonds[0]) == 1 and len(neighbouring_aromatic_bonds[1]) == 1:
                pi_electrons = (2 + len(sp3) + len(double_bonds)) * 2

            else:
                return False
        elif len(neighbouring_aromatic_bonds) == 1:

            pi_electrons = (int(math.ceil(len(neighbouring_aromatic_bonds[0]) / 2)) + len(double_bonds) + len(sp3)) * 2
        else:
            return False

    print(pi_electrons)
    
    if pi_electrons % 4 == 2:
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

def check_five_ring(atom_set):
    assert len(atom_set) == 5

    sp2_hybridised = []
    sp3_hybridised_lone_pair = []

    aromatic = False
    heteroatom = None

    for atom in atom_set:
        if atom.hybridisation == 'sp2':
            sp2_hybridised.append(atom)
        elif atom.hybridisation == 'sp3':
            if atom.calc_electron_pair_nr() > 0:
                sp3_hybridised_lone_pair.append(atom)

    if len(sp2_hybridised) == 4 and len(sp3_hybridised_lone_pair) == 1:

        pi_electron_nr = 0
        for atom in sp2_hybridised:
            for orbital in atom.valence_shell.orbitals.values():
                if orbital.orbital_type == 'p':
                    for electron in orbital.electrons:
                        if electron.atom == atom:
                            pi_electron_nr += 1
        if pi_electron_nr % 4 == 0:
            aromatic = True
            heteroatom = sp3_hybridised_lone_pair[0]

    return aromatic, heteroatom