"""
smarts.py
=========
This module contains the entry point for the SMARTS parser.
"""
import typing as ty

from parsing import parse_smarts 
from compiling import Smarts, compile_smarts

def read_smarts(smarts: str) -> Smarts:
    """
    Parses and compiles a SMARTS string into a graph.

    Parameters
    ----------
    smarts: str
        The SMARTS string to parse and compile.
    
    Returns
    -------
    Smarts
        The compiled SMARTS.
    """
    ast = parse_smarts(smarts)

    smarts = compile_smarts(ast)

    return smarts