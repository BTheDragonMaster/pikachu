import os

from joblib import load
import unittest
import logging
logging.basicConfig(level=logging.DEBUG)

from pikachu.chem.substructure_matching import SubstructureMatch, find_substructures

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
TEST_DIR = os.path.join(CURRENT_DIR, "test_data")


class TestSubstructureMatching(unittest.TestCase):

    def test_find_substructures(self):

        structure_1 = load(os.path.join(TEST_DIR, "structure_1.struct"))
        structure_2 = load(os.path.join(TEST_DIR, "structure_2.struct"))
        # draw_structure(structure_1)
        # draw_structure(structure_2)
        # Assert there is at least one match, even when the stereochemistry is not correct
        self.assertEqual(len(find_substructures(structure_1, structure_2)), 1)


if __name__ == "__main__":
    unittest.main()