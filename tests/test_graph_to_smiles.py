import os

from joblib import load
import unittest
import logging
logging.basicConfig(level=logging.DEBUG)

from pikachu.smiles.graph_to_smiles import structure_to_smiles
from pikachu.general import read_smiles

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
TEST_DIR = os.path.join(CURRENT_DIR, "test_data")


class TestGraphToSmiles(unittest.TestCase):

    def test_structure_to_smiles(self):

        structure_1 = load(os.path.join(TEST_DIR, "structure_1.struct"))
        structure_2 = read_smiles(structure_to_smiles(structure_1))

        self.assertEqual(len(structure_1.find_substructures(structure_2)), 1)

