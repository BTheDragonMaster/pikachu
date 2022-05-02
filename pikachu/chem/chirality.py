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


def same_chirality(order_1, order_2):
    permutations = get_chiral_permutations(order_1)
    if tuple(order_2) in permutations:
        return True
    else:
        return False


def get_chiral_permutations_lonepair(order):
    permutations = [tuple(order),
                    (order[1], order[2], order[0]),
                    (order[2], order[0], order[1])]

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
