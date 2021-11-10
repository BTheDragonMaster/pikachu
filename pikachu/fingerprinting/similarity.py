from pikachu.fingerprinting.ecfp_4 import ECFP


def get_jaccard_index(structure_1, structure_2, fingerprinting_depth=2):
    ecfp_1 = ECFP(structure_1, iterations=fingerprinting_depth)
    ecfp_2 = ECFP(structure_2, iterations=fingerprinting_depth)

    jaccard_index = ecfp_1.fingerprint.intersection(ecfp_2.fingerprint) / ecfp_1.fingerprint.union(ecfp_2.fingerprint)

    return jaccard_index


def get_jaccard_distance(structure_1, structure_2, fingerprinting_depth=2):
    jaccard_index = get_jaccard_index(structure_1, structure_2, fingerprinting_depth=fingerprinting_depth)
    jaccard_distance = 1 - jaccard_index

    return jaccard_distance