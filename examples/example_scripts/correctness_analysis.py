from sys import argv
from pikachu.general import svg_from_smiles


def parse_smiles(tbd_file):
    name_to_compound = {}
    with open(tbd_file, 'r') as tbd:
        for line in tbd:
            line = line.strip()
            if line:
                compound_name, smiles = line.split('\t')
                if compound_name not in name_to_compound:
                    name_to_compound[compound_name] = []
                name_to_compound[compound_name].append(smiles)
    return name_to_compound


if __name__ == "__main__":
    in_file = argv[1]
    name_to_compounds = parse_smiles(in_file)
    for name, compounds in name_to_compounds.items():
        for i, compound in enumerate(compounds):
            out_file = f"{name}_{i}.svg"
            svg_from_smiles(compound, out_file)
