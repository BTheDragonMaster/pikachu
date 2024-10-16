<img src="logo.png" alt="thumbnail" width="600px" />

# INSTALLATION

Python-based Informatics Kit for the Analysis of Chemical Units

Step 1: Make a conda environment:

```bash
conda create -n pikachu python>=3.9
conda activate pikachu
```

Step 2: install pip:

```bash
conda install pip
```

Step 3: Install PIKAChU:

```bash
pip install pikachu-chem
```

# GETTING STARTED

Step 1: Open python or initiate an empty .py file.

Step 2: Import required modules to visualise your first structure:

```Python
from pikachu.general import draw_smiles
```

Step 3: Load your SMILES string of interest and draw it!

```Python
smiles = draw_smiles("CCCCCCCCCC(=O)N[C@@H](CC1=CNC2=CC=CC=C21)C(=O)N[C@@H](CC(=O)N)C(=O)N[C@@H](CC(=O)O)C(=O)N[C@H]3[C@H](OC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@H](NC(=O)CNC(=O)[C@@H](NC(=O)[C@H](NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)CNC3=O)CCCN)CC(=O)O)C)CC(=O)O)CO)[C@H](C)CC(=O)O)CC(=O)C4=CC=CC=C4N)C")
```

Step 4: Play around with the other functions in pikachu.general. For guidance, refer to documentation in the wiki and function descriptors.

# Citation
```
Terlouw, Barbara R., Sophie PJM Vromans, and Marnix H. Medema. "PIKAChU: a Python-based informatics kit for analysing chemical units." Journal of Cheminformatics 14.1 (2022): 34.
```
