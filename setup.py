from setuptools import setup, find_packages

VERSION = '1.0.9'
DESCRIPTION = 'PIKACHU: Python-based Informatics Kit for Analysing CHemical Units'
LONG_DESCRIPTION = 'An easy-to-use cheminformatics kit with few dependencies.'

setup(
    name="pikachu-chem",
    version=VERSION,
    author="Barbara Terlouw",
    author_email="barbara.terlouw@wur.nl",
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    packages=find_packages(),
    install_requires=['matplotlib'],
)
