"""
nodes.py
======
This module contains the AST nodes for parsing SMARTS.
"""
import typing as ty

from tokens import BOND_TOKENS, Token

class Node:
    """
    A node is part of the abstract syntax tree (AST).
    """
    ...

class Statement(Node):
    """
    A statement is a node that can be executed.
    """
    ...

class Expression(Node):
    """
    An expression is a node that can evaluated.
    """
    ...

class Literal(Expression):
    """
    A literal is an expression that represents a value.
    """
    def __init__(self, type: str, value: str) -> None:
        """
        Initialize the literal with its type and value.

        Parameters
        ----------
        type: str
            The type of the literal.
        value: str
            The value of the literal.
        """
        self.type = type
        self.value = value
        
    def __repr__(self) -> str:
        """
        Returns the representation of the literal.
        
        Returns
        -------
        str
            The representation of the literal.
        """
        return f"{self.__class__.__name__}({self.type}, {self.value})"
    
    def type_cast(self) -> ty.Any:
        """
        Type casts the literal to the given type.

        Returns
        -------
        ty.Any
            The type casted value.

        Raises
        ------
        ValueError
            If the literal type is not supported.
        """
        if self.type in [Token.ATOM_SYMBOL, Token.ATOM_ANY]:
            return self.value # TODO: Convert to atom type.
        
        elif self.type in BOND_TOKENS:
            return self.value # TODO: Convert to bond type.
        
        elif self.type == Token.MOD_INT:
            return int(self.value)
        
        else:
            raise ValueError(f"Unknown literal type: {self.type}")
        
    def eval(self) -> ty.Any:
        """
        Evaluates the literal.

        Returns
        -------
        ty.Any
            The evaluated value.
        """
        return self.type_cast()
    
class NestingStatement(Statement):
    """
    A nesting statement is a statement that contains other statements.
    """
    def __init__(self, statements: ty.List[Statement]) -> None:
        """
        Initialize the nesting statement with its statements.

        Parameters
        ----------
        statements : ty.List[Statement]
            The statements of the nesting statement.
        """
        self.statements = statements

    def __repr__(self) -> str:
        """
        Returns the representation of the nesting statement.
        
        Returns
        -------
        str
            The representation of the nesting statement.
        """
        return f"{self.__class__.__name__}({self.statements})"
    
class AtomStatement(Statement):
    """
    Represents an atom.
    """
    def __init__(self, atom: Literal) -> None:
        """
        Initialize the atom statement with its atom.

        Parameters
        ----------
        atom : Literal
            The atom of the atom statement.
        """
        self.atom = atom.eval()

    def __repr__(self) -> str:
        """
        Returns the representation of the atom statement.
        
        Returns
        -------
        str
            The representation of the atom statement.
        """
        return f"{self.__class__.__name__}({self.atom})"
    
class BondStatement(Statement):
    """
    Represents a bond that connects two atoms.
    """
    def __init__(self, bond: Literal) -> None:
        """
        Initialize the bond statement with its bond.

        Parameters
        ----------
        bond: Literal
            The bond of the bond statement.
        """
        self.bond = bond.eval()

    def __repr__(self) -> str:
        """
        Returns the representation of the bond statement.
        
        Returns
        -------
        str
            The representation of the bond statement.
        """
        return f"{self.__class__.__name__}({self.bond})"
    
class ConnectivityStatement(Statement):
    """
    Represents a non-linear connectivity between atoms.
    """
    def __init__(self, connectivity: Literal) -> None:
        """
        Initialize the connectivity statement with its connectivity.

        Parameters
        ----------
        connectivity: Literal
            The connectivity of the connectivity statement.
        """
        self.connectivity = connectivity.eval()

    def __repr__(self) -> str:
        """
        Returns the representation of the connectivity statement.
        
        Returns
        -------
        str
            The representation of the connectivity statement.
        """
        return f"{self.__class__.__name__}({self.connectivity})"
