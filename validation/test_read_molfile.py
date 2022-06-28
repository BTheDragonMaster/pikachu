import os
from pikachu.chem.molfile.read_molfile import MolFileReader

"""
Unit tests for MolFileReader class
"""


class TestMolFileReader:
    # Define test molfile_path and molefile_str
    test_dir = os.path.split(__file__)[0]
    molfile_path = os.path.join(test_dir, 'temp.mol')
    with open(molfile_path, 'r') as molfile:
        molfile_str = molfile.read()

    def test_constructor(self):
        # Make sure instantiation fails when no args are given
        try:
            MolFileReader()
        except Exception as exception:
            assert type(exception) == ValueError
        # Make sure instantiation fails when both args are given
        try:
            MolFileReader(self.molfile_str, self.molfile_path)
        except Exception as exception:
            assert type(exception) == ValueError
        # Make sure instantiation works with molfile_path given
        assert MolFileReader(molfile=self.molfile_path)
        # Make sure instantiation works with molfile_str given
        assert MolFileReader(molfile_str=self.molfile_str)

    def test_get_molfile_lines(self):
        # Assert that reading from file/str makes no difference
        reader_from_path = MolFileReader(molfile=self.molfile_path)
        reader_from_str = MolFileReader(molfile_str=self.molfile_str)
        lines_from_path = reader_from_path.get_molfile_lines()
        lines_from_str = reader_from_str.get_molfile_lines()
        assert lines_from_path == lines_from_str
