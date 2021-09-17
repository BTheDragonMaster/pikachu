from groupdefiner import GroupDefiner

from sys import argv

def parse_smiles(smiles_file):
    name_to_group = {}

    with open(smiles_file, 'r') as smiles:
        smiles.readline()
        for line in smiles:
            line = line.strip()
            smiles_string, name, index = line.split('\t')
            index = int(index)
            name_to_group[name] = GroupDefiner(name, smiles_string, index)

    return name_to_group


   # smiles = open(smiles_file, 'r')
   # smiles.close()

if __name__ == "__main__":
    smiles_file = argv[1] #argv = ['parsers.py', 'data.txt']
    name_to_group = parse_smiles(smiles_file)
    for name, group in name_to_group.items():
        print(name, group.atom_1)


