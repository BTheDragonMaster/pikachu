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
    BondStatement
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
            or self.parse_atom() 
            or self.parse_bond()
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
        
    def parse_atom(self) -> AtomStatement:
        """
        Parses an atom node.
        
        Returns
        -------
        AtomStatement
            The parsed atom node.
        """
        if self.curr_token == Token.ATOM_SYMBOL:
            atom = self.parse_literal()

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