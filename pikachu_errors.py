#!/usr/bin/env python

class SmilesError(Exception):
    error_to_message = {'chiral double bond': "Bond stereochemistry incorrectly defined.",
                        'invalid smiles': "Invalid smiles.",
                        'bond': "Incorrect bond placement."}
    def __init__(self, error_type):
        self.message = self.error_to_message[error_type]


