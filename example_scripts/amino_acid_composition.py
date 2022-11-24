from pikachu.general import (
    smiles_from_file,
    highlight_subsmiles_single,
    highlight_subsmiles_all,
    highlight_subsmiles_multiple,
    highlight_substructure,
)

tryptophan = r"N[C@H](C=O)CC1=CNC2=CC=CC=C21"
d_asparagine = r"N[C@@H](C=O)CC(N)=O"
asparagine = r"N[C@H](C=O)CC(N)=O"
aspartate = r"N[C@H](C=O)CC(O)=O"
threonine = r"N[C@H](C=O)[C@@H](C)[O]"
glycine = r"O=CCN"
ornithine = r"O=C[C@H](CCCN)N"
d_alanine = r"O=C[C@@H](C)N"
d_serine = r"O=C[C@@H](CO)N"
glutamate = r"O=C[C@H](CCC(O)=O)N"
kynurenine = r"O=C[C@H](CC(C1=CC=CC=C1N)=O)N"
d_leucine = r"O=C[C@H](N)CC(C)C"
d_tyrosine = r"C1=CC(=CC=C1C[C@H](C(=O))N)O"
# d_hydroxyphenylglycine = r"OC1=CC=C([C@H](C=O)N)C=C1"
d_hydroxyphenylglycine = r"N[C@@H](C=O)C1=CC=C(O)C=C1"
tyrosine = r"C1=CC(=CC=C1C[C@@H](C(=O))N)O"
dihydroxyphenylglycine = r"N[C@H](C(O)=O)C1=CC(O)=CC(O)=C1"

daptomycin = r"CCCCCCCCCC(=O)N[C@@H](CC1=CNC2=CC=CC=C21)C(=O)N[C@H](CC(=O)N)C(=O)N[C@@H](CC(=O)O)C(=O)N[C@H]3[C@H](OC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@H](NC(=O)CNC(=O)[C@@H](NC(=O)[C@H](NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)CNC3=O)CCCN)CC(=O)O)C)CC(=O)O)CO)[C@H](C)CC(=O)O)CC(=O)C4=CC=CC=C4N)C"
vancomycin = r"C[C@H]1[C@H]([C@@](C[C@@H](O1)O[C@@H]2[C@H]([C@@H]([C@H](O[C@H]2OC3=C4C=C5C=C3OC6=C(C=C(C=C6)[C@H]([C@H](C(=O)N[C@H](C(=O)N[C@H]5C(=O)N[C@@H]7C8=CC(=C(C=C8)O)C9=C(C=C(C=C9[C@H](NC(=O)[C@H]([C@@H](C1=CC(=C(O4)C=C1)Cl)O)NC7=O)C(=O)O)O)O)CC(=O)N)NC(=O)[C@@H](CC(C)C)NC)O)Cl)CO)O)O)(C)N)O"

# daptomycin = smiles_from_file('daptomycin.smi')
highlight_subsmiles_single(
    aspartate,
    daptomycin,
    colour="light blue",
    visualisation="svg",
    out_file="daptomycin_aspartate_single.svg",
)
highlight_subsmiles_all(
    aspartate,
    daptomycin,
    colour="light blue",
    visualisation="svg",
    out_file="daptomycin_aspartate_all.svg",
)
highlight_subsmiles_multiple(
    [aspartate, tryptophan],
    daptomycin,
    colours=["red", "blue"],
    visualisation="svg",
    out_file="daptomycin_multiple.svg",
)

highlight_subsmiles_multiple(
    [
        glycine,
        d_alanine,
        d_serine,
        threonine,
        d_asparagine,
        aspartate,
        glutamate,
        kynurenine,
        ornithine,
        tryptophan,
    ],
    daptomycin,
    colours=[
        "red",
        "blue",
        "orange",
        "hot pink",
        "light blue",
        "dark blue",
        "dark red",
        "purple",
        "lime",
        "yellow",
    ],
    visualisation="svg",
    out_file="daptomycin_substructures.svg",
)
highlight_subsmiles_multiple(
    [
        asparagine,
        d_leucine,
        d_hydroxyphenylglycine,
        dihydroxyphenylglycine,
        d_tyrosine,
        tyrosine,
    ],
    vancomycin,
    colours=["red", "blue", "hot pink", "lime", "light blue", "yellow"],
    visualisation="svg",
    out_file="vancomycin_substructures.svg",
)

highlight_substructure(
    aspartate,
    daptomycin,
    search_mode="single",
    colour="light blue",
    visualisation="svg",
    out_file="daptomycin_aspartate_single.svg",
)
highlight_substructure(
    aspartate,
    daptomycin,
    search_mode="all",
    colour="light blue",
    visualisation="svg",
    out_file="daptomycin_aspartate_all.svg",
)
highlight_substructure(
    [aspartate, tryptophan], daptomycin, search_mode="multiple", colour=["red", "blue"]
)
