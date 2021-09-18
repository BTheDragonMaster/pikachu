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
