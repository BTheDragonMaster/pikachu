"""
repl.py
=======
This module contains the REPL for SMARTS.
"""
from sys import argv

from smarts import read_smarts

def main() -> None:
    """
    The main function for the REPL.
    """
    smarts_str = argv[1]

    smarts_obj = read_smarts(smarts_str)

    print(smarts_obj)

    exit(0)

if __name__ == "__main__":
    main()
