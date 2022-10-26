def smiles_from_npatlas_tabular(npatlas_file):
    npaid_to_smiles = {}
    with open (npatlas_file, 'r') as npatlas:
        npatlas.readline()
        for line in npatlas:
            line = line.strip()
            compound_info = line.split('\t')
            print(compound_info)
            npaid = compound_info[0]
            smiles = compound_info[10]
            npaid_to_smiles[npaid] = smiles

    return npaid_to_smiles