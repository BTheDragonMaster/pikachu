
from distutils.core import setup

setup(name='pikachu',
      version='1.0',
      authors=["Barbara Terlouw"],
      packages=["pikachu"],
      # Add dependencies here!
      install_requires=["scipy", "numpy"]
      # scripts=["bin/caretta-app", "bin/caretta-cli"] add scripts to run from the command-line here
)
