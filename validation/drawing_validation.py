from sys import argv

from pikachu.general import read_smiles
from pikachu.general import smiles_to_molfile
from pikachu.drawing.drawing import Options

from rdkit.Chem.rdmolfiles import MolFromMolFile, MolToMolFile
from rdkit.Chem import MolToSmiles, MolFromSmiles


def validate_pikachu(smiles_file, finetune=True, strict_mode=False):
    correct = 0
    incorrect = 0
    total = 0

    incorrect_smiles = []

    options = Options()
    options.finetune = finetune
    options.strict_mode = strict_mode

    with open(smiles_file, "r") as smiles_f:
        with open("failed_smiles.smi", "w") as failed_smiles:
            for line in smiles_f:
                smiles = line.strip()

                if total % 100 == 0:
                    print(f"Correct SMILES: {correct}")
                    print(f"Incorrect SMILES: {incorrect}")
                    print(f"Total SMILES: {total}")

                total += 1

                try:
                    # pikachu_structure = read_smiles(smiles)
                    rdkit_structure = MolFromSmiles(smiles)

                    MolToMolFile(rdkit_structure, "temp_rdkit.mol")

                    smiles_to_molfile(smiles, "temp.mol", options=options)

                    # MolFileWriter(pikachu_structure, 'temp.mol', drawing_options=options).write_mol_file()

                    rdkit_structure_pikachu = MolFromMolFile("temp.mol")

                    rdkit_smiles_original = MolToSmiles(rdkit_structure)
                    rdkit_smiles_pikachu = MolToSmiles(rdkit_structure_pikachu)

                    if rdkit_smiles_original == rdkit_smiles_pikachu:
                        correct += 1
                    else:
                        incorrect += 1
                        incorrect_smiles.append(
                            (smiles, rdkit_smiles_original, rdkit_smiles_pikachu)
                        )
                        print("Original:", smiles)
                        print("RDKit:   ", rdkit_smiles_original)
                        print("PIKAChU: ", rdkit_smiles_pikachu)

                except Exception as e:
                    print(smiles, e)
                    failed_smiles.write(f"{smiles}\t{e}\n")

                if total == 100000:
                    break

    print(f"Correct SMILES: {correct}")
    print(f"Incorrect SMILES: {incorrect}")
    print(f"Total SMILES: {total}")

    return incorrect_smiles


if __name__ == "__main__":
    smiles_file = argv[1]
    out_file = argv[2]
    strict_mode = False
    finetune = False
    if len(argv) > 3:
        strict_mode = bool(int(argv[3]))
    if len(argv) > 4:
        finetune = bool(int(argv[4]))

    incorrect_smiles = validate_pikachu(
        smiles_file, finetune=finetune, strict_mode=strict_mode
    )

    with open(out_file, "w") as out:
        out.write(
            "Original SMILES\tCanonical RDKit SMILES original\tCanonical RDKit SMILES PIKAChU\n"
        )
        for smiles, rdkit_smiles_original, rdkit_smiles_pikachu in incorrect_smiles:
            out.write(f"{smiles}\t{rdkit_smiles_original}\t{rdkit_smiles_pikachu}\n")
