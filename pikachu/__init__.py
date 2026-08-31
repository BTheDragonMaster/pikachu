from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("pikachu-chem")
except PackageNotFoundError:
    __version__ = "unknown"