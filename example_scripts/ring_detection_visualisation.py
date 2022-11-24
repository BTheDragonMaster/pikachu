from pikachu.general import highlight_subsmiles_multiple, read_smiles

daptomycin = r"CCCCCCCCCC(=O)N[C@@H](CC1=CNC2=CC=CC=C21)C(=O)N[C@H](CC(=O)N)C(=O)N[C@@H](CC(=O)O)C(=O)N[C@H]3[C@H](OC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@H](NC(=O)CNC(=O)[C@@H](NC(=O)[C@H](NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)CNC3=O)CCCN)CC(=O)O)C)CC(=O)O)CO)[C@H](C)CC(=O)O)CC(=O)C4=CC=CC=C4N)C"
vancomycin = r"C[C@H]1[C@H]([C@@](C[C@@H](O1)O[C@@H]2[C@H]([C@@H]([C@H](O[C@H]2OC3=C4C=C5C=C3OC6=C(C=C(C=C6)[C@H]([C@H](C(=O)N[C@H](C(=O)N[C@H]5C(=O)N[C@@H]7C8=CC(=C(C=C8)O)C9=C(C=C(C=C9[C@H](NC(=O)[C@H]([C@@H](C1=CC(=C(O4)C=C1)Cl)O)NC7=O)C(=O)O)O)O)CC(=O)N)NC(=O)[C@@H](CC(C)C)NC)O)Cl)CO)O)O)(C)N)O"

s = read_smiles(vancomycin)
print(s.aromatic_cycles)
daptomycin_substructures = [
    "c1ccccc1",
    "c1c[nH]c2c1cccc2",
    "C1CNCCNCCNCCNCCNCCNCCNCCNCCNCCOC1",
]
vancomycin_substructures = [
    "C1CCCCO1",
    "C(Cc1ccc(O)cc1)CNC",
    "CNC(c1ccccc1)",
    "CNCCNC(c1ccccc1)",
    "NC(c1ccccc1)",
]  # "C(Cc1ccc2cc1)CNCCNC(c3cc(O2)cc4c3)CNC(c5cccc6c5)CNC(Cc7ccc(O4)cc7)CNCc8c6cccc8"]

highlight_subsmiles_multiple(
    daptomycin_substructures,
    daptomycin,
    colours=["blue", "hot pink", "blue"],
    visualisation="svg",
    out_file="daptomycin_rings.svg",
    check_chiral_centres=False,
)
highlight_subsmiles_multiple(
    vancomycin_substructures,
    vancomycin,
    colours=["blue", "pink", "pink", "pink", "pink"],
    visualisation="svg",
    out_file="vancomycin_rings.svg",
    check_chiral_centres=False,
)
