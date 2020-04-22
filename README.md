# PIKACHU

Python-based Informatics Kit for the Analysis of Chemical Units

Step 1: Clone the repository:

```git clone https://git.wur.nl/terlo012/pikachu.git```

Step 2: Navigate to folder containing pikachu.py

Step 3: Open python

Step 4: Import required modules

```
import pikachu
import reactions
from pprint import pprint
```

Step 5: Load your SMILES string of interest and turn it into a structure

```
smiles = pikachu.Smiles("CCCCCCCCCC(=O)N[C@@H](CC1=CNC2=CC=CC=C21)C(=O)N[C@@H](CC(=O)N)C(=O)N[C@@H](CC(=O)O)C(=O)N[C@H]3[C@H](OC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@H](NC(=O)CNC(=O)[C@@H](NC(=O)[C@H](NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)CNC3=O)CCCN)CC(=O)O)C)CC(=O)O)CO)[C@H](C)CC(=O)O)CC(=O)C4=CC=CC=C4N)C")
structure = smiles.smiles_to_structure()
pprint(structure.graph)
```

Step 6: Now, you can manipulate your structure with the reactions in reactions.py (only hydrolysis has been fully tested). First, define the bond you want to break with BondDefiner:

```
peptide_bond = reactions.BondDefiner('peptide_bond', 'C(=O)NC', 0, 2)
```
This function takes a custom bond name, a SMILES string of the bond you want to (in this case) hydrolise, an the atom indices of the two atoms between which the bond needs to be hydrolysed.

Now, find occurrences of this bond in your structure:
```
bonds = reactions.find_bonds(structure, peptide_bond)
```

Then, you can hydrolyse any or all of these bonds with the hydrolysis function. Note: this edits the structure in place, so you might want to make a copy of the original structure first.

```
import copy

product = copy.deepcopy(structure)
for bond in bonds:
    reactions.hydrolyse(bond, product)
```

Finally, you can split the resulting (possibly disconnected) graph(s) into separate graphs and visualise them:

```
products = product.split_disconnected_structures()
for product in products:
    pprint(product)
```