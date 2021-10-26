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
            for orbital_name in atom.valence_shell.orbitals:
                orbital = atom.valence_shell.orbitals[orbital_name]
                if orbital.orbital_type == 'p':
                    for electron in orbital.electrons:
                        if electron.atom == atom:
                            pi_electron_nr += 1

        if not pi_electron_nr % 4 == 2:
            aromatic = False

    return aromatic


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
            for orbital_name in atom.valence_shell.orbitals:
                orbital = atom.valence_shell.orbitals[orbital_name]
                if orbital.orbital_type == 'p':
                    for electron in orbital.electrons:
                        if electron.atom == atom:
                            pi_electron_nr += 1
        if pi_electron_nr % 4 == 0:
            aromatic = True
            heteroatom = sp3_hybridised_lone_pair[0]

    return aromatic, heteroatom