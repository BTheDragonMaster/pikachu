class BondProperties:
    """
    A class storing various properties of bonds

    Attributes
    ----------
    bond_type_to_weight: dict of {bond type: bond_weight, ->}
        with bond type (str) the type of covalent bond, and bond weight (int)
        the number of times that bond is counted in determining if an atom is
        in its excited state or not
    bond_type_to_symbol: dict of {bond type: SMILES symbol, ->}
        with bond type (str) the type of covalent bond, and SMILES symbol the
        text symbol used in SMILES strings to represent that bond
    bond_type_to_p_orbitals: dict of {bond type: p orbitals, ->}
        with bond type (str) the type of covalent bond, and p orbitals (int)
        the number of p orbitals involved in the formation of that bond
    """

    type_to_dash2d_input = {'single': 1,
                            'double': 2,
                            'triple': 3,
                            'quadruple': 4}

    bond_type_to_weight = {'single': 1,
                           'double': 2,
                           'triple': 3,
                           'quadruple': 4,
                           'aromatic': 1}

    bond_type_to_symbol = {'single': '',
                           'double': '=',
                           'triple': '#',
                           'aromatic': ''}

    bond_type_to_p_orbitals = {'single': 0,
                               'double': 1,
                               'triple': 2,
                               'quadruple': 3,
                               'aromatic': 1}

    bond_type_to_order = {'single': 1,
                          'double': 2,
                          'triple': 3,
                          'aromatic': 4,
                          'quadruple': 5}


BOND_PROPERTIES = BondProperties()