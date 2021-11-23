from pikachu.errors import SmilesError
from pikachu.chem.chirality import get_chiral_permutations


class SubstructureMatch:
    def __init__(self, child, parent, seed_child, seed_parent):
        self.child = child
        self.parent = parent

        self.seed_child = seed_child
        self.seed_parent = seed_parent

        self.bonds = {}
        self.atoms = {}

        self.atoms_reversed = {}

        self.initialise_atom_matches()
        self.initialise_bond_matches()

        self.current_atom_child = self.seed_child
        self.current_atom_parent = self.seed_parent
        self.next_atom_child = None
        self.next_atom_parent = None

        self.current_bond_child = None
        self.current_bond_parent = None

        self.atoms_to_place = set()
        self.bonds_to_place = set()

        self.placed_bonds_and_atoms = [(None, self.current_atom_child)]
        self.attempt_path = []

        self.atom_to_bond_count = {}
        self.initialise_bond_counts()

        self.initialise_atoms_to_place()
        self.initialise_bonds_to_place()

        self.parent_bond_to_visited = {}
        self.parent_bond_to_matching_attempts = {}
        self.parent_atom_to_bonds = {}

        self.failed_attempt_paths = set()

        self.initialise_bonds_to_visited()

        self.active = True
        self.complete = False

    def initialise_bond_counts(self):
        for atom in self.child.get_atoms():
            if atom.type != 'H':
                self.atom_to_bond_count[atom] = len(atom.get_non_hydrogen_neighbours())

    def find_next_matching_atoms(self):
        # Get parent bonds that could potentially still match

        candidate_parent_bonds = self.get_next_bonds_parent()
        if self.current_atom_parent not in self.parent_atom_to_bonds:
            self.parent_atom_to_bonds[self.current_atom_parent] = []
        for bond in candidate_parent_bonds:
            if bond.bond_summary == self.current_bond_child.bond_summary:
                next_child_atom_candidate = self.current_bond_child.get_connected_atom(self.current_atom_child)
                next_parent_atom_candidate = bond.get_connected_atom(self.current_atom_parent)
                # Make sure the atom matching is possible: either the match must already exist, or if it doesn't,
                # both parent atom and child atom must be free to match to one another.
                if self.atoms[next_child_atom_candidate] == next_parent_atom_candidate \
                        or (self.atoms[next_child_atom_candidate] == None and next_parent_atom_candidate not in self.atoms.values()):

                    if next_parent_atom_candidate.potential_same_connectivity(next_child_atom_candidate.connectivity):
                        self.current_bond_parent = bond
                        self.next_atom_child = next_child_atom_candidate
                        self.next_atom_parent = next_parent_atom_candidate
                        if bond not in self.parent_atom_to_bonds[self.current_atom_parent]:
                            self.parent_atom_to_bonds[self.current_atom_parent].append(bond)

    def initialise_bonds_to_visited(self):
        for bond in self.parent.bonds.values():
            if 'H' not in [atom.type for atom in bond.neighbours]:
                self.parent_bond_to_visited[bond] = False
                self.parent_bond_to_matching_attempts[bond] = set()

    def initialise_atoms_to_place(self):
        for atom in self.child.get_atoms():
            if atom.type != 'H' and atom != self.seed_child:
                self.atoms_to_place.add(atom)

    def initialise_bonds_to_place(self):
        for bond in self.child.bonds.values():
            if 'H' not in [atom.type for atom in bond.neighbours]:
                self.bonds_to_place.add(bond)

    def initialise_atom_matches(self):
        for atom in self.child.graph:

            # only match non-hydrogen atoms
            if atom.type != 'H':
                self.atoms[atom] = None

        self.atoms[self.seed_child] = self.seed_parent
        self.atoms_reversed[self.seed_parent] = self.seed_child

    def initialise_bond_matches(self):
        for bond in self.child.bonds.values():
            # only select bonds not connected to any hydrogens
            if 'H' not in [atom.type for atom in bond.neighbours]:
                self.bonds[bond] = None

    def get_next_bonds_parent(self):
        bonds = []
        for parent_bond in self.current_atom_parent.bonds:
            # Ignore bonds that are attached to a hydrogen
            # Only consider bonds that haven't been matched to another bond already
            if parent_bond not in self.bonds.values() and 'H' not in [atom.type for atom in parent_bond.neighbours]:
                if not self.parent_bond_to_visited[parent_bond]:
                    bonds.append(parent_bond)
                # makes sure a matching attempt wasn't made before between parent bond and child bond.
                # Problem with this: what if the bond was correct, but the route there wasn't?
                else:
                    attempt_path = tuple(self.attempt_path + [parent_bond])
                    if attempt_path not in self.failed_attempt_paths:
                        bonds.append(parent_bond)

               # elif self.current_bond_child not in self.parent_bond_to_matching_attempts[parent_bond]:
                #    bonds.append(parent_bond)

        return bonds

    # STILL NOT GOOD ENOUGH! ONLY LOOK AT THE LAST ATTEMPT!
    # STORE ATTEMPT PATHS?

    def traceback(self):
        new_parent_candidate = None
        new_child_candidate = None

        for i in range(len(self.placed_bonds_and_atoms) - 1, -1, -1):

            current_child_bond, current_child_atom = self.placed_bonds_and_atoms[i]
            current_parent_atom = self.atoms[current_child_atom]

            try:
                for option in self.parent_atom_to_bonds[current_parent_atom]:

                    if not self.parent_bond_to_visited[option]:
                        new_parent_candidate = current_parent_atom
                        new_child_candidate = current_child_atom
                        break

            except KeyError:
                pass

            if new_parent_candidate:
                break
            else:

                del self.placed_bonds_and_atoms[i]

                try:
                    del self.attempt_path[i - 1]
                except IndexError:
                    pass

                if self.attempt_path:
                    self.failed_attempt_paths.add(tuple(self.attempt_path))


                self.remove_match(current_child_atom, current_child_bond)

        return new_child_candidate, new_parent_candidate

    def get_next_bonds_child(self):
        candidate_child_bonds = []

        # Only consider bonds that still have to be placed

        for bond in self.current_atom_child.bonds:
            if bond in self.bonds_to_place:
                candidate_child_bonds.append(bond)

        return candidate_child_bonds

    def add_match(self):
        self.atoms[self.next_atom_child] = self.next_atom_parent
        self.atoms_reversed[self.next_atom_parent] = self.next_atom_child

        # It is possible the atom has already been placed! Ring closure.
        if self.next_atom_child in self.atoms_to_place:
            self.atoms_to_place.remove(self.next_atom_child)

        self.atom_to_bond_count[self.next_atom_child] -= 1

        previous_child_atom = self.current_bond_child.get_connected_atom(self.next_atom_child)
        self.atom_to_bond_count[previous_child_atom] -= 1

        self.bonds[self.current_bond_child] = self.current_bond_parent
        self.parent_bond_to_visited[self.current_bond_parent] = True
        self.parent_bond_to_matching_attempts[self.current_bond_parent].add(self.current_bond_child)
        self.bonds_to_place.remove(self.current_bond_child)
        self.current_atom_child = self.next_atom_child
        self.current_atom_parent = self.next_atom_parent
        self.placed_bonds_and_atoms.append((self.current_bond_child, self.next_atom_child, ))
        self.attempt_path.append(self.current_bond_parent)

    def remove_match(self, child_atom, child_bond):

        parent_atom_to_remove = self.atoms[child_atom]

        remove_atom = True

        # Make sure that atom is only removed from a match if all bonds with that atom are removed

        for bond in self.bonds:
            # Only consider a bond if it has matched to a parent bond
            if bond and self.bonds[bond]:
                # Only consider bonds neighbouring the child atom
                if child_atom in bond.neighbours:
                    # Only consider bonds different than the one that is currently being removed
                    if bond != child_bond:
                        remove_atom = False
                        break

        if remove_atom:
            self.atoms[child_atom] = None
            del self.atoms_reversed[parent_atom_to_remove]

        self.atom_to_bond_count[child_atom] += 1

        if child_bond:
            self.bonds[child_bond] = None
            previous_child_atom = child_bond.get_connected_atom(child_atom)
            self.atom_to_bond_count[previous_child_atom] += 1
            self.bonds_to_place.add(child_bond)
        else:
            pass

        self.atoms_to_place.add(child_atom)

    def check_complete(self):
        if not self.atoms_to_place and not self.bonds_to_place:
            self.complete = True
        return self.complete


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


def compare_all_matches(matches): # refactor to 'filter_duplicate_matches'
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


def filter_duplicate_matches(matches): # refactor to 'filter_duplicate_matches'
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
    print(equivalent_atom_list)
    print(atom_2.neighbours)

    if len(equivalent_atom_list) != 4:
        lone_pairs = atom_2.lone_pairs

        try:
            assert len(equivalent_atom_list) + len(lone_pairs) == 4
        except AssertionError:
            raise SmilesError('chiral centre')

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

        match = SubstructureMatch(child, structure, starting_atom, seed)

        # keeps track of bond traversal: if there are bonds left to place, continue trying to match.
        # If the match becomes inactive because no exact matches are found, abort

        while match.bonds_to_place and match.active:

            match.current_bond_parent = None
            match.next_atom_child = None
            match.next_atom_parent = None
            candidate_child_bonds = match.get_next_bonds_child()

            # Are there still outgoing bonds that haven't been traversed from the current child atom?

            if candidate_child_bonds:

                # Choose the first one

                match.current_bond_child = candidate_child_bonds[0]
                match.find_next_matching_atoms()

                if match.current_bond_parent:
                    match.add_match()

                else:
                    attempt_path = tuple(match.attempt_path)
                    match.failed_attempt_paths.add(attempt_path)

                    new_child_candidate, new_parent_candidate = match.traceback()

                    if not new_child_candidate:
                        match.active = False
                    else:
                        match.current_atom_child = new_child_candidate
                        match.current_atom_parent = new_parent_candidate

            else:
                if match.bonds_to_place:

                    # Happens if one or multiple cycles need closing somewhere.

                    if None not in match.atoms.values():
                        bond = list(match.bonds_to_place)[0]

                        match.current_atom_child = bond.atom_1
                        match.current_atom_parent = match.atoms[match.current_atom_child]

                    # Happens if we jump to a different bond elsewhere

                    else:
                        # Make sure we only consider bonds that are connected to an atom that is currently matched

                        for bond in match.bonds_to_place:

                            if match.atoms[bond.atom_1] or match.atoms[bond.atom_2]:
                                # decide which atom to 'look' from

                                if match.atoms[bond.atom_1]:
                                    match.current_atom_child = bond.atom_1
                                else:
                                    match.current_atom_child = bond.atom_2

                                match.current_atom_parent = match.atoms[match.current_atom_child]
                                break

                        else:
                            match.active = False

                else:
                    # This only happens if a match is complete, upon which we escape from the while loop.
                    pass

        if match.active:
            matches.append(match)

    filter_duplicate_matches(matches)

    return matches