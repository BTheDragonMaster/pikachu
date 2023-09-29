"""
compiling.py
============
This module contains the compiler for compiling a parsed SMARTS into a graph.
"""
import typing as ty
from collections import defaultdict

from nodes import (
    NestingStatement,
    AtomStatement,
    BondStatement,
    ConnectivityStatement
)
from parsing import AST

class SmartsAtom:
    """
    A node in a SMARTS structure graph.
    """
    def __init__(self, symbol: str) -> None:
        """
        Initialize the SMARTS atom.
        """
        self.symbol = symbol

    def __repr__(self) -> str:
        """
        Returns the representation of the SMARTS atom.

        Returns
        -------
        str
            The representation of the SMARTS atom.
        """
        return f"{self.__class__.__name__}({self.symbol!r})"
    
class SmartsBond:
    """
    An edge in a SMARTS structure graph.
    """
    def __init__(self, symbol: str) -> None:
        """
        Initialize the SMARTS bond.
        """
        self.symbol = symbol
    
    def __repr__(self) -> str:
        """
        Returns the representation of the SMARTS bond.

        Returns
        -------
        str
            The representation of the SMARTS bond.
        """
        return f"{self.__class__.__name__}({self.symbol!r})"

class Smarts:
    """
    A SMARTS object.
    """
    def __init__(self) -> None:
        """
        Initialize the SMARTS object.
        """
        graph = defaultdict(lambda: defaultdict(list))

    def __repr__(self) -> str:
        """
        Returns the representation of the SMARTS object.

        Returns
        -------
        str
            The representation of the SMARTS object.
        """
        return f"{self.__class__.__name__}()"

class Compiler:
    """
    A compiler compiles a parsed SMARTS into a graph.
    """
    def __init__(self, ast: NestingStatement) -> None:
        """
        Initialize the compiler with an AST.

        Parameters
        ----------
        ast: NestingStatement
            The AST to compile.
        """
        self.ast = ast

    def compile(self) -> Smarts:
        """
        Compiles the AST into a graph.
        
        Returns
        -------
        Smarts
            The compiled SMARTS.
        """
        smarts = Smarts()

        # TODO: Compile the AST into a graph.

        return smarts

def compile_smarts(ast: AST) -> Smarts:
    """
    Compiles a parsed SMARTS into a graph.

    Parameters
    ----------
    ast: AST
        The AST to compile.

    Returns
    -------
    Smarts
        The compiled SMARTS.
    """
    compiler = Compiler(ast)

    smarts = compiler.compile()

    return smarts