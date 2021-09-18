#!/usr/bin/env python

class SmilesError(Exception):
    error_to_message = {'chiral double bond': "Conflicting double bond stereochemistry.",
                        'invalid smiles': "Invalid smiles.",
                        'bond': "Incorrect bond placement.",
                        'violated_bonding_laws': "Basic bonding laws have been violated.",
                        'chiral centre': "Non-chiral atom defined as chiral."}

    def __init__(self, error_type):
        self.message = self.error_to_message[error_type]


