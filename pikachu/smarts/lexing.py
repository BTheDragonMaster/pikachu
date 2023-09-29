"""
lexing.py
========
This module contains the lexer for parsing SMARTS.
"""
import typing as ty

from tokens import EOL, ILLEGAL, Token, TokenInfo

class Lexer:
    """
    Lexer for tokenizing SMARTS string.
    """
    def __init__(self, src: str) -> None:
        """
        Initialize the lexer with the source string.

        Parameters
        ----------
        src: str
            The source string to be tokenized.
        """
        self.src = src
        self.pos = 0
        self.offset = 0

    def lex(self) -> ty.Generator[TokenInfo, None, None]:
        """
        Tokenizes the SMARTS string.

        Yields
        ------
        TokenInfo
            Token information object that contains information about token. 
        """
        while self.pos < len(self.src):
            for token_id in Token:
                if match := token_id.value.match(self.src, self.pos):
                    col_s, col_e = self.pos, match.end(0)
                    self.pos = col_e

                    yield TokenInfo(
                        token_id.name, 
                        match.group(0),
                        col_s - self.offset, 
                        col_e - self.offset
                    )
                    break 

            else:
                yield TokenInfo(
                    ILLEGAL,
                    self.src[self.pos],
                    self.pos - self.offset, 
                    self.pos - self.offset + 1
                )
                self.pos += 1 

        else:
            # Yield EOL twice as the parser always looks ahead two tokens.
            col = self.pos - self.offset + 1
            yield TokenInfo(EOL, "\x00", col, col)
            yield TokenInfo(EOL, "\x00", col, col)
