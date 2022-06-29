from pikachu.general import read_smiles
from pikachu.drawing.drawing import Drawer


class TestDrawer:
    def test_set_r_group_indices_subscript(self):
        # Test that only R group indices are returned as subscript
        dummy_structure = read_smiles("CC")
        drawer = Drawer(dummy_structure)
        test_str = ['Xe',
                    'C',
                    'O',
                    'R',
                    'R1',
                    'X23',
                    'Z54']
        expected_str = ['Xe',
                        'C',
                        'O',
                        'R',
                        'R₁',
                        'X₂₃',
                        'Z₅₄']
        for index in range(len(test_str)):
            result = drawer.set_r_group_indices_subscript(test_str[index])
            expected = expected_str[index]
            assert result == expected
