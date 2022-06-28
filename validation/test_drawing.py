from pikachu.general import read_smiles
from pikachu.drawing.drawing import Drawer


class TestDrawer:
    def test_set_r_group_indices_subscript(self):
        # Test that only R group indices are returned as subscript
        dummy_structure = read_smiles("CC")
        drawer = Drawer(dummy_structure)
        test_str = ['c1',
                    'C1',
                    'C',
                    'O',
                    '13',
                    'R',
                    'R1',
                    'X23',
                    'Z54']
        expected_str = ['c1',
                        'C1',
                        'C',
                        'O',
                        '13',
                        'R',
                        'R₁',
                        'X₂₃',
                        'Z₅₄']
        for index in range(len(test_str)):
            result = drawer.set_r_group_indices_subscript(test_str[index])
            expected = expected_str[index]
            assert result == expected
