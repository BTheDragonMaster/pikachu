from pikachu.fingerprinting.ecfp_4 import ECFP


def get_jaccard_index(structure_1, structure_2, fingerprinting_depth=2):
    ecfp_1 = ECFP(structure_1, iterations=fingerprinting_depth)
    ecfp_2 = ECFP(structure_2, iterations=fingerprinting_depth)

    jaccard_index = len(ecfp_1.fingerprint.intersection(ecfp_2.fingerprint)) / len(ecfp_1.fingerprint.union(ecfp_2.fingerprint))

    return jaccard_index


def get_jaccard_from_ecfp(ecfp_1, ecfp_2):
    jaccard_index = len(ecfp_1.fingerprint.intersection(ecfp_2.fingerprint)) / len(
        ecfp_1.fingerprint.union(ecfp_2.fingerprint))
    jaccard_distance = 1 - jaccard_index

    return jaccard_distance


def get_jaccard_distance(structure_1, structure_2, fingerprinting_depth=2):
    jaccard_index = get_jaccard_index(structure_1, structure_2, fingerprinting_depth=fingerprinting_depth)
    jaccard_distance = 1 - jaccard_index

    return jaccard_distance


def get_jaccard_matrix(name_to_compound, fingerprinting_depth=2):
    name_to_ecfp = {}
    for name, compound in name_to_compound.items():
        ecfp = ECFP(compound, iterations=fingerprinting_depth)
        name_to_ecfp[name] = ecfp

    matrix = {}
    for name_1, ecfp_1 in name_to_ecfp.items():
        if name_1 not in matrix:
            matrix[name_1] = {}
        for name_2, ecfp_2 in name_to_ecfp.items():
            jaccard_distance = get_jaccard_from_ecfp(ecfp_1, ecfp_2)
            matrix[name_1][name_2] = jaccard_distance

    return matrix
