from pikachu.errors import StructureError
from pikachu.chem.chirality import get_chiral_permutations


class SubstructureMatch:
    def __init__(self, child, parent):
        self.atoms = {}
        self.bonds = {}

        self.child = child
        self.parent = parent

        self.parent_atom_to_bonds = {}

        self.initialise_match()

        self.current_child_atom = None
        self.current_parent_atom = None

        self.current_attempted_path = []
        self.failed_attempted_paths = set()

        self.active = True

    def initialise_match(self):
        for atom in self.child.graph:
            if atom.type != 'H':
                self.atoms[atom] = None

        for bond in self.child.bonds.values():
            if not bond.has_neighbour('H'):
                self.bonds[bond] = None

        for atom in self.parent.graph:
            if atom.type != 'H':
                self.parent_atom_to_bonds[atom] = []

    def initialise_seed(self, child_seed, parent_seed):
        self.atoms[child_seed] = parent_seed
        self.current_child_atom = child_seed
        self.current_parent_atom = parent_seed
        self.current_attempted_path.append(parent_seed)

    def match(self, child_seed, parent_seed):
        self.initialise_seed(child_seed, parent_seed)
        counter = 0
        while None in self.bonds.values() and self.active:
            counter += 1
            next_child_bond = self.get_next_child_bond()
            if next_child_bond:
                next_child_atom, next_parent_atom, next_parent_bond = self.find_next_matching_atoms(next_child_bond)

                if next_child_atom and next_parent_atom and next_parent_bond:
                    self.add_match(next_child_atom, next_parent_atom, next_child_bond, next_parent_bond)
                else:
                    self.failed_attempted_paths.add(tuple(self.current_attempted_path))
                    next_child_atom, next_parent_atom, next_child_bond, next_parent_bond = self.traceback()
                    if next_child_atom and next_parent_atom and next_child_bond and next_parent_bond:
                        self.add_match(next_child_atom, next_parent_atom, next_child_bond, next_parent_bond)
                    else:
                        self.active = False
            else:
                self.active = False

    def traceback(self):
        # Iterate over the current attempted path in reverse

        for i in range(len(self.current_attempted_path) - 1, 0, -1):

            current_parent_atom = self.current_attempted_path[i]
            previous_parent_atom = self.current_attempted_path[i - 1]

            if current_parent_atom == 'hop' or previous_parent_atom == 'hop':
                self.failed_attempted_paths.add(tuple(self.current_attempted_path[:i]))
                continue

            current_child_atom = None
            previous_child_atom = None

            for child_atom, parent_atom in self.atoms.items():
                if current_parent_atom == parent_atom:
                    current_child_atom = child_atom
                if previous_parent_atom == parent_atom:
                    previous_child_atom = child_atom

            assert current_child_atom and previous_child_atom

            # This can happen because of structure hopping: in this case, you want to move on to the next atom,
            # as no alternative candidate bonds exists from the previous child atom

            if not self.child.bond_exists(current_child_atom, previous_child_atom):
                continue
            else:
                next_child_bond = self.child.bond_lookup[current_child_atom][previous_child_atom]
                self.remove_match(current_child_atom, next_child_bond)

                new_path = self.current_attempted_path[:i]
                candidate_bonds = self.parent_atom_to_bonds[previous_parent_atom]

                for bond in candidate_bonds:
                    # Make sure the bond is not already currently matched to another bond
                    if bond not in self.bonds.values():
                        next_parent_atom = bond.get_connected_atom(previous_parent_atom)
                        next_child_atom = current_child_atom

                        # Make sure that we did not try this route before

                        if tuple(new_path + [next_parent_atom]) not in self.failed_attempted_paths:
                            self.current_child_atom = previous_child_atom
                            self.current_parent_atom = previous_parent_atom
                            next_parent_bond = bond
                            self.current_attempted_path = new_path
                            return next_child_atom, next_parent_atom, next_child_bond, next_parent_bond

                self.failed_attempted_paths.add(tuple(new_path))

        return None, None, None, None

    def add_match(self, next_child_atom, next_parent_atom, next_child_bond, next_parent_bond):
        self.atoms[next_child_atom] = next_parent_atom
        self.bonds[next_child_bond] = next_parent_bond

        self.current_child_atom = next_child_atom
        self.current_parent_atom = next_parent_atom

        self.current_attempted_path.append(next_parent_atom)

    def remove_match(self, child_atom, child_bond):

        self.bonds[child_bond] = None
        remove_atom_from_match = True
        for child_bond, parent_bond in self.bonds.items():
            if parent_bond:
                if child_bond.atom_1 == child_atom or child_bond.atom_2 == child_atom:
                    remove_atom_from_match = False

        if remove_atom_from_match:
            self.atoms[child_atom] = None

    def find_next_matching_atoms(self, child_bond):
        candidate_parent_bonds = self.get_candidate_parent_bonds()
        for parent_bond in candidate_parent_bonds:
            if parent_bond.bond_summary == child_bond.bond_summary:
                next_child_atom = child_bond.get_connected_atom(self.current_child_atom)
                next_parent_atom = parent_bond.get_connected_atom(self.current_parent_atom)

                # Make sure the atom matching is possible: either the match must already exist, or if it doesn't,
                # both parent atom and child atom must be free to match to one another.
                if self.atoms[next_child_atom] == next_parent_atom or \
                        (not self.atoms[next_child_atom] and next_parent_atom not in self.atoms.values()):

                    if next_parent_atom.potential_same_connectivity(next_child_atom.connectivity):
                        if parent_bond not in self.parent_atom_to_bonds[self.current_parent_atom]:
                            self.parent_atom_to_bonds[self.current_parent_atom].append(parent_bond)

                        return next_child_atom, next_parent_atom, parent_bond

        return None, None, None

    # TODO: Set explicit tag so hydrogen matching can also be done.

    def get_candidate_parent_bonds(self):
        candidate_parent_bonds = []
        for parent_bond in self.current_parent_atom.bonds:
            # Ignore bonds that are attached to a hydrogen
            # Only consider bonds that haven't been matched to another bond already
            if parent_bond not in self.bonds.values() and not parent_bond.has_neighbour('H'):
                # makes sure a matching attempt wasn't made before between parent bond and child bond.
                next_parent_atom = parent_bond.get_connected_atom(self.current_parent_atom)
                attempted_path = tuple(self.current_attempted_path + [next_parent_atom])
                if attempted_path not in self.failed_attempted_paths:
                    candidate_parent_bonds.append(parent_bond)

        return candidate_parent_bonds

    def get_next_child_bond(self):
        for child_bond in self.current_child_atom.bonds:
            if not child_bond.has_neighbour('H') and not self.bonds[child_bond]:
                return child_bond

        # Happens if one or multiple cycles need closing somewhere.

        if None not in self.atoms.values():
            for child_bond, parent_bond in self.bonds.items():
                if not parent_bond:

                    self.current_child_atom = child_bond.atom_1
                    self.current_parent_atom = self.atoms[self.current_child_atom]

                    self.current_attempted_path.append('hop')

                    self.current_attempted_path.append(self.current_parent_atom)

                    return child_bond

        # Happens if we jump to a different bond elsewhere, but still connected to currently matched atoms

        else:

            # Make sure we only consider bonds that are connected to an atom that is currently matched

            for child_bond, parent_bond in self.bonds.items():
                if not parent_bond:
                    if self.atoms[child_bond.atom_1] or self.atoms[child_bond.atom_2]:
                        if self.atoms[child_bond.atom_1]:
                            self.current_child_atom = child_bond.atom_1
                        elif self.atoms[child_bond.atom_2]:
                            self.current_child_atom = child_bond.atom_2

                        self.current_parent_atom = self.atoms[self.current_child_atom]
                        self.current_attempted_path.append('hop')
                        self.current_attempted_path.append(self.current_parent_atom)

                        return child_bond

        return None


def compare_matches(match_1, match_2):
    matching = True
    for key in match_1:
        if key not in match_2:
            matching = False
            break

        if match_1[key] != match_2[key]:
            matching = False
            break

    return matching


def is_same_match(match_1, match_2):
    matching = True
    for atom in match_1.atoms:
        if atom not in match_2.atoms:
            matching = False
            break

        if match_1.atoms[atom] != match_2.atoms[atom]:
            matching = False
            break

    if matching:
        for bond in match_1.bonds:
            if bond not in match_2.bonds:
                matching = False
                break

            if match_1.bonds[bond] != match_2.bonds[bond]:
                matching = False
                break

    return matching


def compare_all_matches(matches):  # refactor to 'filter_duplicate_matches'
    matching_pairs = set()

    for i, match_1 in enumerate(matches):
        for j, match_2 in enumerate(matches):
            if i != j:
                if compare_matches(match_1, match_2):
                    matching_pairs.add(tuple(sorted([i, j])))

    matches_to_remove = set()

    for matching_pair in matching_pairs:
        matches_to_remove.add(matching_pair[1])

    matches_to_remove = sorted(list(matches_to_remove), reverse=True)

    for match_to_remove in matches_to_remove:
        del matches[match_to_remove]


def filter_duplicate_matches(matches):  # refactor to 'filter_duplicate_matches'
    matching_pairs = set()

    for i, match_1 in enumerate(matches):
        for j, match_2 in enumerate(matches):
            if i != j:
                if is_same_match(match_1, match_2):
                    matching_pairs.add(tuple(sorted([i, j])))

    matches_to_remove = set()

    for matching_pair in matching_pairs:
        matches_to_remove.add(matching_pair[1])

    matches_to_remove = sorted(list(matches_to_remove), reverse=True)

    for match_to_remove in matches_to_remove:
        del matches[match_to_remove]


def check_same_chirality(atom_1, atom_2, match):
    equivalent_atom_list = []
    for atom in atom_1.neighbours:
        if atom.type == 'H':
            for atom_b in atom_2.neighbours:
                if atom_b.type == 'H':
                    equivalent_atom_list.append(atom_b)
                    break
        else:
            equivalent_atom_list.append(match[atom])

    permutation = equivalent_atom_list[:]

    if len(equivalent_atom_list) != 4:
        lone_pairs = atom_2.lone_pairs

        try:
            assert len(equivalent_atom_list) + len(lone_pairs) == 4
        except AssertionError:
            raise StructureError('chiral centre')

        permutation += lone_pairs

    chiral_permutations = get_chiral_permutations(permutation)

    if atom_1.chiral == atom_2.chiral:
        if tuple(atom_2.neighbours) in chiral_permutations:
            return True
        else:
            return False
    else:
        if tuple(atom_2.neighbours) in chiral_permutations:
            return False
        else:
            return True


def find_substructures(structure, child):
    atom_connectivities_child = child.get_connectivities()
    atom_connectivities_parent = structure.get_substructure_connectivities(atom_connectivities_child)
    # Sort based on the complexity of the connectivity

    connectivities = sorted(list(atom_connectivities_child.keys()),
                            key=lambda x: len(set(x)), reverse=True)

    starting_connectivity = connectivities[0]

    starting_atom = atom_connectivities_child[starting_connectivity][0]

    seeds = []

    for atom in atom_connectivities_parent[starting_connectivity]:
        if atom.type == starting_atom.type:
            seeds.append(atom)

    matches = []

    for seed in seeds:

        match = SubstructureMatch(child, structure)
        match.match(starting_atom, seed)

        if match.active:
            matches.append(match)

    filter_duplicate_matches(matches)

    return matches
