from pikachu.errors import SmilesError
from pikachu.chem.chirality import get_chiral_permutations


def compare_matches(match_1, match_2):
    matching = True
    for key in match_1:
        if key not in match_2:
            matching = False
            break

        if match_1[key] != match_2[key]:
            matching = False
            break

    return matching


def compare_all_matches(matches):
    matching_pairs = set()

    for i, match_1 in enumerate(matches):
        for j, match_2 in enumerate(matches):
            if i != j:
                if compare_matches(match_1, match_2):
                    matching_pairs.add(tuple(sorted([i, j])))

    matches_to_remove = set()

    for matching_pair in matching_pairs:
        matches_to_remove.add(matching_pair[1])

    matches_to_remove = sorted(list(matches_to_remove), reverse=True)

    for match_to_remove in matches_to_remove:
        del matches[match_to_remove]


def check_same_chirality(atom_1, atom_2, match):

    equivalent_atom_list = []
    for atom in atom_1.neighbours:
        if atom.type == 'H':
            for atom_2 in atom_2.neighbours:
                if atom_2.type == 'H':
                    equivalent_atom_list.append(atom_2)
                    break
        else:
            equivalent_atom_list.append(match[atom])

    permutation = equivalent_atom_list[:]

    if len(equivalent_atom_list) != 4:
        lone_pairs = atom_2.lone_pairs

        try:
            assert len(equivalent_atom_list) + len(lone_pairs) == 4
        except AssertionError:
            raise SmilesError('chiral centre')

        permutation += lone_pairs

    chiral_permutations = get_chiral_permutations(permutation)

    if atom_1.chiral == atom_2.chiral:
        if tuple(permutation) in chiral_permutations:
            return True
        else:
            return False
    else:
        if tuple(permutation) in chiral_permutations:
            return False
        else:
            return True
