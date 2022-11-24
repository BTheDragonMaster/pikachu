from sys import argv

from rdkit.Chem import MolFromSmiles, MolToSmiles
from pikachu.general import read_smiles, structure_to_smiles, smiles_from_file
import timeout_decorator


@timeout_decorator.timeout(10)
def read_and_write(smiles):
    pikachu_smiles = structure_to_smiles(read_smiles(smiles))
    return pikachu_smiles


def assess_writing_accuracy(smiles_strings):
    wrong_smiles = []
    correct_smiles = 0
    failed_smiles = []
    for i, smiles in enumerate(smiles_strings):
        try:
            pikachu_smiles = read_and_write(smiles)

            rdkit_original = MolToSmiles(MolFromSmiles(smiles))
            rdkit_pikachu = MolToSmiles(MolFromSmiles(pikachu_smiles))

            if pikachu_smiles:
                if rdkit_pikachu != rdkit_original:
                    print("Mistake:")
                    print(f"Original: {smiles}")
                    print(f"PIKAchU:  {pikachu_smiles}")
                    print(f"RDKit:    {rdkit_original}")
                    print(f"rd_pika:  {rdkit_pikachu}")
                    wrong_smiles.append(smiles)
                else:
                    correct_smiles += 1
            else:
                failed_smiles.append(smiles)

        except Exception as e:
            failed_smiles.append(smiles)
            print(smiles, e)

        if i % 1000 == 0:
            print(f"Correct SMILES: {correct_smiles}")
            print(f"Incorrect SMILES: {len(wrong_smiles)}")
            print(f"Failed SMILES: {len(failed_smiles)}")
            print(f"Total SMILES: {i}")

    print(f"Correct SMILES: {correct_smiles}")
    print(f"Incorrect SMILES: {len(wrong_smiles)}")
    print(f"Failed SMILES: {len(failed_smiles)}")
    print(f"Total SMILES: {len(smiles_strings)}")


if __name__ == "__main__":
    smiles = smiles_from_file(argv[1], all=True)
    assess_writing_accuracy(smiles)
