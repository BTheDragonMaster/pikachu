from sys import argv


def parse_coconut(coconut_file):
    smiles_strings = []
    with open(coconut_file, 'r') as coconut:
        for line in coconut:
            line = line.strip()
            smiles, identifier = line.split()
            smiles_strings.append(smiles)

    return smiles_strings


def write_smiles(smiles, out_file):
    with open(out_file, 'w') as out:
        for smi in smiles:
            out.write(f"{smi}\n")


if __name__ == "__main__":
    smiles = parse_coconut(argv[1])
    write_smiles(smiles, argv[2])
