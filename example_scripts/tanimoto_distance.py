from pikachu.fingerprinting.similarity import get_jaccard_matrix
from pikachu.general import read_smiles, svg_from_structure

from sys import argv
import os

import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from sklearn.manifold import TSNE

PROTEINOGENIC = {
    "alanine",
    "cysteine",
    "aspartate",
    "glutamate" "aspartic acid",
    "glutamic acid",
    "phenylalanine",
    "glycine",
    "histidine",
    "isoleucine",
    "lysine",
    "leucine",
    "methionine",
    "asparagine",
    "proline",
    "glutamine",
    "arginine",
    "serine",
    "threonine",
    "valine",
    "tryptophan",
    "tyrosine",
}


def plot_tanimoto_distances(matrix):
    figure(figsize=(10, 10), dpi=80)
    sorted_compounds = sorted(list(matrix.keys()))
    distances = []

    for compound_1 in sorted_compounds:
        distance = []
        for compound_2 in sorted_compounds:
            distance.append(matrix[compound_1][compound_2])

        distances.append(distance)

    tsne = TSNE(n_components=2, metric="precomputed", random_state=6)
    results = tsne.fit(distances)

    coords = results.embedding_

    plt.subplots_adjust(bottom=0.1)
    plt.scatter(coords[:, 0], coords[:, 1], marker="o")
    for label, x, y in zip(sorted_compounds, coords[:, 0], coords[:, 1]):
        plt.annotate(
            label,
            xy=(x, y),
            xytext=(0, 20),
            textcoords="offset points",
            ha="right",
            va="bottom",
            bbox=dict(boxstyle="round,pad=0.5", fc="grey", alpha=0.2),
            arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=0"),
        )

    plt.savefig("tanimoto.svg")


def get_coors(coor_dict):
    """Gets x and y coordinates from a dict with coordinates as values

    Keyword arguments:
        coor_dict -- dict, values are lists of floats, with the x and y
            coordinate as first and second float.

    Returns:
        x_list -- list, x coordinates from the input dictionary.
        y_list -- list, y coordinates from the input dictionary.
    """
    sorted_compounds = sorted(list(coor_dict.keys()))
    x = []
    y = []
    for compound in sorted_compounds:

        x.append(coor_dict[compound][0])
        y.append(coor_dict[compound][1])

    return sorted_compounds, x, y


def make_tsne_matrix(matrix):
    """Creates dissimilarity matrix and tsne coordinates, and saves the latter.

    Keyword arguments:
        list_mols -- list of lists, inner list contains the molecule id and
            the fingerprint of the molecule.

    Returns:
        tsne_dict -- dict, keys are molecule_ids, values are lists with the
            coordinates of that point on the tsne plot.
    """
    sorted_compounds = sorted(list(matrix.keys()))
    distances = []

    for compound_1 in sorted_compounds:
        distance = []
        for compound_2 in sorted_compounds:
            distance.append(matrix[compound_1][compound_2])

        distances.append(distance)

    tsne_dict = {}

    res_tsne = TSNE(n_components=2, metric="precomputed", random_state=0).fit_transform(
        distances
    )

    for i in range(len(sorted_compounds)):
        tsne_dict[sorted_compounds[i]] = list(res_tsne[i, :])

    print(tsne_dict)
    print(len(tsne_dict.keys()))

    return tsne_dict


def plot_tsne(tsne_vals):
    """Makes a plot of the t-SNE results and saves the plot.

    Keyword arguments:
        tsne_vals -- dict, keys are molecule_ids, values are lists with the
            coordinates of that point on the tsne plot.
        save_folder -- string, folder where to save the t-SNE plot.
        *sub_dict -- optional dicts, values are lists of floats, with the x
            and y coordinate as first and second float. These coordinates are
            plotted in a different color on the t-SNE plot.
    """

    sorted_compounds, x_coors, y_coors = get_coors(tsne_vals)
    print(sorted_compounds, x_coors, y_coors)
    plt.scatter(x_coors, y_coors)

    for label, x, y in zip(sorted_compounds, x_coors, y_coors):
        plt.annotate(
            label,
            xy=(x, y),
            xytext=(0, 20),
            textcoords="offset points",
            ha="right",
            va="bottom",
            bbox=dict(boxstyle="round,pad=0.5", fc="grey", alpha=0.2),
            arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=0"),
        )

    plt.savefig("tanimoto.svg")


#   plt.show()


def parse_smiles(tbd_file):
    name_to_compound = {}
    with open(tbd_file, "r") as tbd:
        tbd.readline()
        for line in tbd:
            line = line.strip()
            if line:
                compound_name, smiles = line.split("\t")
                try:
                    structure = read_smiles(smiles)
                except Exception:
                    print(f"Couldn't convert {compound_name}.")
                if structure:
                    name_to_compound[compound_name] = structure
                else:
                    print(f"Couldn't convert {compound_name}.")
    return name_to_compound


def write_network(matrix, out_file):
    with open(out_file, "w") as out:
        out.write("Compound 1\tDistance\tCompound 2\n")
        for compound_1 in matrix:
            for compound_2 in matrix[compound_1]:
                if compound_1 != compound_2:
                    out.write(
                        f"{compound_1.title()}\t{1 - matrix[compound_1][compound_2]}\t{compound_2.title()}\n"
                    )


def is_amino_acid(name_to_compound, out_dir):
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    out_file = os.path.join(out_dir, "substrate_identities.txt")
    with open(out_file, "w") as out:
        out.write("Substrate\tamino acid\n")
        for name, structure in name_to_compound.items():
            is_amino_acid = False
            is_beta_amino_acid = False
            amino_acid = read_smiles("NCC(O)=O")
            beta_amino_acid = read_smiles("NCCC(O)=O")
            if structure.find_substructures(amino_acid):
                is_amino_acid = True
            elif structure.find_substructures(beta_amino_acid):
                is_beta_amino_acid = True

            if is_amino_acid:
                if name.lower() in PROTEINOGENIC:
                    out.write(f"{name.title()}\tproteinogenic\n")
                else:
                    out.write(f"{name.title()}\tamino_acid\n")
            elif is_beta_amino_acid:
                out.write(f"{name.title()}\tbeta_amino_acid\n")
            else:
                out.write(f"{name.title()}\tno_amino_acid\n")

            svg_from_structure(structure, os.path.join(out_dir, f"{name}.svg"))


if __name__ == "__main__":
    tbd_file = argv[1]
    out_file = argv[2]
    out_2 = argv[3]
    name_to_compound = parse_smiles(tbd_file)
    matrix = get_jaccard_matrix(name_to_compound)
    write_network(matrix, out_file)
    is_amino_acid(name_to_compound, out_2)

# plot_tanimoto_distances(matrix)
# tsne_dict = make_tsne_matrix(matrix)
# plot_tsne(tsne_dict)
