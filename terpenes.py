import pikachu
from pprint import pprint

if __name__ == "__main__":
    F_11 = "C/C1=C\\CC(C)(C)[CH+]CC/C(C)=C/CC1"
    N_11 = "C/C1=C/CC(C)(C)[CH+]CC/C(C)=C/CC1"

    F_11_structure = pikachu.Smiles(F_11).smiles_to_structure()
    N_11_structure = pikachu.Smiles(N_11).smiles_to_structure()

    for bond_nr, bond in F_11_structure.bonds.items():
        if bond.chiral:
            print(bond)
            pprint(bond.chiral_dict)

    for bond_nr, bond in N_11_structure.bonds.items():
        if bond.chiral:
            print(bond)
            pprint(bond.chiral_dict)
            
    
    
