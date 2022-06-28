class LonePair:
    def __init__(self, atom, nr):
        self.nr = nr
        self.parent = atom
        self.electrons = []

    def __hash__(self):
        return self.nr
    
    def __repr__(self):
        return f"LonePair_{self.parent.nr}_{self.nr - 10000}"
    
    def add_electron(self, electron):
        self.electrons.append(electron)
