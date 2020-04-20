import structure_inference.smiles.Smiles as Smiles
import structure_inference.structure.Structure as Structure
import structure_inference.bond.Bond as Bond
import structure_inference.atom.Atom as Atom
import structure_inference.electron.Electron as Electron

import importlib
import pkgutil


def import_submodules(package, recursive=True):
    """ Import all submodules of a module, recursively, including subpackages

    :param package: package (name or actual module)
    :type package: str | module
    :rtype: dict[str, types.ModuleType]
    """
    if isinstance(package, str):
        package = importlib.import_module(package)
    results = {}
    for loader, name, is_pkg in pkgutil.walk_packages(package.__path__):
        full_name = package.__name__ + '.' + name
        results[full_name] = importlib.import_module(full_name)
        if recursive and is_pkg:
            results.update(import_submodules(full_name))
    return results

#import_submodules(__name__)
