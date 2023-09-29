"""
parsing.py
=========
This module contains the parser for parsing SMARTS.
"""
import typing as ty

from tokens import EOL, BOND_TOKENS, Token
from lexing import Lexer
from nodes import (
    Literal, 
    Expression,
    Statement,
    NestingStatement,
    AtomStatement,
    BondStatement,
    ConnectivityStatement
)

class Parser:
    """
    A parser parses SMARTS into an AST.
    """
    def __init__(self, lexer: Lexer) -> None:
        """
        Initialize the parser with a lexer.

        Parameters
        ----------
        lexer: Lexer
            The lexer to use for parsing.
        """
        self.lex = iter(lexer.lex())
        self.curr_token = None 
        self.next_token = None 
        self.update() # Initialize self.next_token
        self.update() # Initialize self.curr_token

    def update(self) -> None:
        """
        Updates the current and next token.
        """
        self.curr_token = self.next_token
        self.next_token = next(self.lex)

    def __iter__(self) -> "Parser":
        """
        Returns the parser as an iterator.
        """
        return self
    
    def __next__(self) -> Statement:
        """
        Returns the next statement.
        
        Returns
        -------
        Statement
            The next statement.
        """
        while self.curr_token.name != EOL:
            if statement := self.parse_node():
                return statement 
        else:
            return EOL
        
    def is_next(self, tokens: ty.List[Token]) -> None:
        """
        Checks if the next token is one of the given tokens.
        
        Parameters
        ----------
        tokens: ty.List[Token]
            The tokens to check.
        
        Raises  
        ------
        ValueError
            If the next token is not one of the given tokens.
        """
        if self.next_token.name not in tokens:
            expected_tokens = [str(t) for t in tokens]
            msg = "Expected one of [" + ", ".join(expected_tokens) + f"] but got {self.next_token}"
            raise ValueError(msg)
        
        self.update()

    def parse_literal(self) -> Expression:
        """
        Parses a literal.
        
        Returns
        -------
        Expression
            The parsed literal.
        """
        return Literal(self.curr_token.name, self.curr_token.val)
    
    def parse_node(self) -> Statement:
        """
        Entry point for parsing.

        Returns
        -------
        Statement
            The parsed statement.
        """
        return (
            self.parse_nesting() 
            or self.parse_atom_without_props() 
            or self.parse_atom_with_props()
            or self.parse_bond()
            or self.parse_connectivity()
            or self.parse_unknown()
        )
    
    def parse_unknown(self) -> None:
        """
        Parses an unknown node.

        Raises
        ------
        ValueError
            If the token is unknown.
        """
        raise ValueError(f"Unknown token {self.curr_token}")
    
    def parse_nesting(self) -> NestingStatement:
        """
        Parses a nesting node.
        
        Returns
        -------
        NestingStatement
            The parsed nesting node.
        """
        if self.curr_token == Token.LPAREN:
            self.update()

            nested = []
            while self.curr_token != Token.RPAREN:
                nested.append(self.parse_node())
            self.update()

            return NestingStatement(nested)
        
    def parse_atom_without_props(self) -> AtomStatement:
        """
        Parses an atom node without properties.
        
        Returns
        -------
        AtomStatement
            The parsed atom node.
        """
        if self.curr_token in [Token.ATOM_SYMBOL, Token.ATOM_ANY]:
            atom = self.parse_literal()
            self.update() 
            
            return AtomStatement(atom)
        
    def parse_atom_with_props(self) -> AtomStatement:
        """
        Parses an atom node with properties.

        Returns
        -------
        AtomStatement
            The parsed atom node.
        """
        if self.curr_token == Token.LSQUARE:
            self.is_next([Token.ATOM_SYMBOL, Token.ATOM_ANY]) 

            atom = self.parse_literal()
            self.update()

            props = dict()
            while self.curr_token != Token.RSQUARE:
                # TODO: Parse atom properties.
                self.update()
            self.update()

            return AtomStatement(atom)

    def parse_bond(self) -> BondStatement:
        """
        Parses a bond node.
        
        Returns
        -------
        BondStatement
            The parsed bond node.
        """
        if self.curr_token in BOND_TOKENS:
            bond = self.parse_literal()

            self.update() 
            
            return BondStatement(bond)
        
    def parse_connectivity(self) -> ConnectivityStatement:
        """
        Parses a connectivity node.
        
        Returns
        -------
        ConnectivityStatement
            The parsed connectivity node.
        """
        if self.curr_token == Token.MOD_INT:
            connectivity = self.parse_literal()
            self.update() 
            
            return ConnectivityStatement(connectivity)

class AST:
    """
    An AST represents a parsed SMARTS.
    """
    def __init__(self, root: ty.List[Statement]) -> None:
        """
        Initialize the AST with a root node.

        Parameters
        ----------
        root: ty.List[Statement]
            The root node of the AST.
        """
        self.root = root

def parse_smarts(src: str) -> AST:
    """
    Parses the given SMARTS into an AST.

    Parameters
    ----------
    src: str
        The SMARTS to parse.

    Returns
    -------
    AST
        The parsed AST.
    """
    lexer = Lexer(src)
    parser = Parser(lexer)

    nodes = []
    for node in parser:
        if node != EOL:
            nodes.append(node)
        else:
            break 
    
    return AST(nodes)