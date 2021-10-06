from hashlib import sha256

class Inchi:
    def __init__(self, inchi):
        self.inchi = inchi
        self.layers = {'version': None,
                       'formula': None,
                       'atom connections': None,
                       'hydrogens': None,
                       'charge': None,
                       'protons': None,
                       'double bond stereochemistry': None,
                       'tetrahedral stereochemistry': None,
                       'allene stereochemistry': None,
                       'stereochemistry information': None,
                       'isotopic layer'}

    def get_layers(self):
        layers = self.inchi.split(r'/')
        version = layers[0].split('InChI=')[-1]
        version = version.split('S')[0]
        self.layers['version'] = version
        self.layers['formula'] = layers[1]

        fixed_h_layer = False

        for layer in layers[2:]:
            prefix = layer[0]

            if prefix == 'f':


            if prefix == 'c':
                self.layers['atom connections'] = layer[1:]
            elif prefix == 'h':
                self.layers['hydrogens'] = layer[1:]


    def inchi_to_structure(self):

