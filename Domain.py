import sys


class Domain:

    def __init__(self, name, parent, start, end):
        self.name = name
        self.parent = parent
        self.start = start
        self.end = end

    def get_motif_format(self):
        return [self.start, self.end, "[]", None, 10, "black", "rgradient:blue", f"arial|8|white|{self.name}"]

    def __str__(self):
        return f"{self.name}, {self.parent}, {self.start}, {self.end}"
