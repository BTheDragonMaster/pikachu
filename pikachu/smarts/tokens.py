"""
tokens.py
=========
This module contains the tokens used by the lexer and parser for parsing SMARTS.
"""
import re 
import typing as ty 
from dataclasses import dataclass 
from enum import Enum 

EOL     = "EOL"
ILLEGAL = "ILLEGAL"

@dataclass
class TokenInfo:
    """
    Token information.
    
    Parameters
    ----------
    name : str
        Name of the token.  
    val : str
        Value of the token. 
    col_s : int
        Column where the token starts in the SMARTS string.  
    col_e : int
        Column where the token ends in the SMARTS string.     
    """
    name: str
    val: str 
    col_s: int 
    col_e: int 

    def __eq__(self, other: ty.Union["TokenInfo", str]) -> bool:
        """
        Checks if other TokenInfo or string representation of TokenInfo is equal to self.
        
        Parameters
        ----------
        other : Union[TokenInfo, str]
            Other TokenInfo or string representation of TokenInfo.
        
        Returns
        -------
        bool
            True if equal, False otherwise.
        """
        if isinstance(other, str):
            return self.name == other 
        
        else:
            return self.name == other.name
        
    def __hash__(self) -> int:
        """
        Returns hash of TokenInfo.
        
        Returns
        -------
        int
            Hash of TokenInfo.
        """
        return id(self.name)
    
    def __repr__(self) -> str:
        """
        Returns string representation of TokenInfo.
        
        Returns
        -------
        str
            String representation of TokenInfo.
        """
        return f"\"{self.name}\" (\"{self.val}\") from col {self.col_s} to {self.col_e}"
    
class Token(Enum):
    """
    All tokens used by the lexer and parser for parsing SMARTS.

    More information on the tokens can be found at:
    https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html
    """
    # Atom qualifier for symbol: 
    # https://www.johndcook.com/blog/2016/02/04/regular-expression-to-match-a-chemical-element/ 
    ATOM_SYMBOL = re.compile((
        r"A[cglmrstu]"
        r"|B[aehikr]?"
        r"|C[adeflmnorsu]?"
        r"|D[bsy]|E[rsu]"
        r"|F[elmr]?"
        r"|G[ade]"
        r"|H[efgos]?"
        r"|I[nr]?"
        r"|Kr?"
        r"|L[airuv]"
        r"|M[dgnot]"
        r"|N[abdeiop]?"
        r"|Os?"
        r"|P[abdmortu]?"
        r"|R[abefghnu]"
        r"|S[bcegimnr]?"
        r"|T[abcehilm]"
        r"|U(u[opst])?"
        r"|V"
        r"|W"
        r"|Xe"
        r"|Yb?"
        r"|Z[nr]"
    ))

    # Other atom qualifiers.
    ATOM_ANY                = re.compile(r"(\*)")   
    ATOM_AROMATIC           = re.compile(r"(a)")
    ATOM_ALIPHATIC          = re.compile(r"(A)")
    ATOM_DEGREE             = re.compile(r"(D)")
    ATOM_HCOUNT_EXPL        = re.compile(r"(H)")
    ATOM_HCOUNT_IMPL        = re.compile(r"(h)")
    ATOM_IN_RING            = re.compile(r"(R)")
    ATOM_IN_RING_SIZE       = re.compile(r"(r)")
    ATOM_VALENCE            = re.compile(r"(v)")
    ATOM_CONN               = re.compile(r"(X)")
    ATOM_CONN_RING          = re.compile(r"(x)")

    # Atom modifiers.
    MOD_AMPERSAND           = re.compile(r"(&)")
    MOD_AT                  = re.compile(r"(@)")
    MOD_CARET               = re.compile(r"(\^)")
    MOD_COMMA               = re.compile(r"(,)")
    MOD_DOLLAR              = re.compile(r"(\$)")
    MOD_INT                 = re.compile(r"(\d+)")

    # Bonds.
    BOND_SINGLE_OR_MIN      = re.compile(r"(-)")
    BOND_DOUBLE             = re.compile(r"(=)")
    BOND_TRIPLE             = re.compile(r"(\#)")
    BOND_AROMATIC_OR_COLON  = re.compile(r"(:)")
    BOND_UP                 = re.compile(r"(/)")
    BOND_DOWN               = re.compile(r"(\\)")
    BOND_RING               = re.compile(r"(@)")
    BOND_ANY                = re.compile(r"(~)")
    DISCONNECTED            = re.compile(r"(\.)")

    # Brackets.
    LPAREN                  = re.compile(r"(\()")
    RPAREN                  = re.compile(r"(\))")
    LSQUARE                 = re.compile(r"(\[)")
    RSQUARE                 = re.compile(r"(\])")

    def __eq__(self, other: ty.Union["Token", str]) -> bool:
        """
        Checks if other Token or string representation of Token is equal to self.

        Parameters
        ----------
        other : Union[Token, str]
            Other Token or string representation of Token.
        
        Returns
        -------
        bool
            True if equal, False otherwise.
        """
        if isinstance(other, str):
            return self.name == other 
        
        else:
            return self.name == other.name
        
    def __hash__(self) -> int:
        """
        Returns hash of Token.
        
        Returns
        -------
        int
            Hash of Token.
        """
        return id(self.name)
    
    def __str__(self) -> str:
        """
        Returns string representation of Token.
        
        Returns
        -------
        str
            String representation of Token.
        """
        return self.name
    
    def __repr__(self) -> str:
        """
        Returns string representation of Token.
        
        Returns
        -------
        str
            String representation of Token.
        """
        return self.name
    
BOND_TOKENS = [
    Token.BOND_SINGLE_OR_MIN,
    Token.BOND_DOUBLE,
    Token.BOND_TRIPLE,
    Token.BOND_AROMATIC_OR_COLON,
    Token.BOND_UP,
    Token.BOND_DOWN,
    Token.BOND_RING,
    Token.BOND_ANY
]
