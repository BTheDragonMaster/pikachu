from sys import argv

def smiles_from_chembl(chembl_file, smiles_out_file):
    chembl_smiles = []
    with open(chembl_file, 'r') as chembl:
        chembl.readline()
        for line in chembl:
            line = line.strip()
            smiles = line.split('\t')[1]
            chembl_smiles.append(smiles)

    chembl_smiles.sort(key=lambda x: len(x))

    with open(smiles_out_file, 'w') as smiles_out:
        for smiles in chembl_smiles:
            smiles_out.write(f'{smiles}\n')

def smiles_from_chembl(chembl_file, smiles_out_file):


if __name__ == "__main__":
    smiles_from_chembl(argv[1], argv[2])

