def get_chiral_permutations(order):
    permutations = [tuple(order),
                    (order[0], order[3], order[1], order[2]),
                    (order[0], order[2], order[3], order[1]),
                    (order[1], order[0], order[3], order[2]),
                    (order[1], order[2], order[0], order[3]),
                    (order[1], order[3], order[2], order[0]),
                    (order[2], order[0], order[1], order[3]),
                    (order[2], order[3], order[0], order[1]),
                    (order[2], order[1], order[3], order[0]),
                    (order[3], order[0], order[2], order[1]),
                    (order[3], order[1], order[0], order[2]),
                    (order[3], order[2], order[1], order[0])]

    return permutations


def find_chirality_from_nonh(neighbours, order, chirality):
    permutations = get_chiral_permutations(neighbours)
    for permutation in permutations:

        if tuple(permutation[:3]) == tuple(order):
            return chirality

    if chirality == 'counterclockwise':
        return 'clockwise'
    else:
        return 'counterclockwise'
