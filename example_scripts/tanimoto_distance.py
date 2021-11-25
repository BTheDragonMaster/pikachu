from pikachu.fingerprinting.similarity import get_jaccard_matrix
from pikachu.general import read_smiles

from sys import argv

import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from sklearn.manifold import TSNE


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
    plt.scatter(
        coords[:, 0], coords[:, 1], marker='o'
        )
    for label, x, y in zip(sorted_compounds, coords[:, 0], coords[:, 1]):
        plt.annotate(
            label,
            xy=(x, y), xytext=(0, 20),
            textcoords='offset points', ha='right', va='bottom',
            bbox=dict(boxstyle='round,pad=0.5', fc='grey', alpha=0.2),
            arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'))

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

    res_tsne = TSNE(n_components=2, metric='precomputed',
                    random_state=0).fit_transform(distances)

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
            xy=(x, y), xytext=(0, 20),
            textcoords='offset points', ha='right', va='bottom',
            bbox=dict(boxstyle='round,pad=0.5', fc='grey', alpha=0.2),
            arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'))

    plt.savefig("tanimoto.svg")
 #   plt.show()


def parse_smiles(tbd_file):
    name_to_compound = {}
    with open(tbd_file, 'r') as tbd:
        for line in tbd:
            line = line.strip()
            if line:
                compound_name, smiles = line.split('\t')
                structure = read_smiles(smiles)
                if structure:
                    name_to_compound[compound_name] = structure
                else:
                    print(f"Couldn't convert {compound_name}.")
    return name_to_compound

#def get_jaccard_matrix(name_to_compound):
#    matrix = {}
 #   for name_1, compound_1 in name_to_compound.items():
  #      if name_1 not in matrix:
  #          matrix[name_1] = {}
  #      for name_2, compound_2 in name_to_compound.items():
  #          jaccard_distance = get_jaccard_distance(compound_1, compound_2)
  #          matrix[name_1][name_2] = jaccard_distance

 #   return matrix

if __name__ == "__main__":
    tbd_file = argv[1]
    name_to_compound = parse_smiles(tbd_file)
    matrix = get_jaccard_matrix(name_to_compound)
   # plot_tanimoto_distances(matrix)
    tsne_dict = make_tsne_matrix(matrix)
    plot_tsne(tsne_dict)
