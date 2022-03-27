#!/usr/bin/env python

class SmilesError(Exception):
    error_to_message = {'chiral double bond': "Conflicting double bond stereochemistry.",
                        'invalid smiles': "Invalid smiles.",
                        'bond': "Incorrect bond placement.",
                        'violated_bonding_laws': "Basic bonding laws have been violated.",
                        'chiral centre': "Non-chiral atom defined as chiral."}

    def __init__(self, error_type):
        self.message = self.error_to_message[error_type]


class DrawingError(Exception):
    error_to_message = {'chiral bond ring': "PIKAChU could not correctly draw the cis/trans stereochemistry of a double bond in a cycle.",
                        'chiral center': "Only three elements attached to chiral center, including hydrogens and lone pairs."}

    def __init__(self, error_type):
        self.message = self.error_to_message[error_type]


class ColourError(Exception):
    def __init__(self, colour):
        if type(colour) == str:
            if colour == 'too few colours':
                self.message = f"Pikachu has too few colours to work with."

            else:
                self.message = f"Pikachu is unfamiliar with the colour {colour}. \n\nPikachu is confused. \nIt hurt itself in confusion."



