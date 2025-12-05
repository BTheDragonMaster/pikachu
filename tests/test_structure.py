import os

from joblib import load
import unittest
import logging
logging.basicConfig(level=logging.DEBUG)

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
TEST_DIR = os.path.join(CURRENT_DIR, "test_data")


class TestStructure(unittest.TestCase):

    def test_find_substructures(self):

        structure_1 = load(os.path.join(TEST_DIR, "structure_1.struct"))
        structure_2 = load(os.path.join(TEST_DIR, "structure_2.struct"))

        self.assertEqual(len(structure_1.find_substructures(structure_2)), 0)

        # Assert there is a match when the stereochemistry is not considered
        self.assertEqual(len(structure_1.find_substructures(structure_2, check_chiral_centres=False)), 1)