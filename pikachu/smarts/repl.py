"""
repl.py
=======
This module contains the REPL for SMARTS.
"""
from sys import argv

from tokens import EOL
from lexing import Lexer 
from parsing import Parser 

def main() -> None:
    """
    The main function for the REPL.
    """
    src = argv[1]

    for node in Parser(Lexer(src)):
        print(node)

        if node == EOL:
            break 

    exit(0)

if __name__ == "__main__":
    main()
