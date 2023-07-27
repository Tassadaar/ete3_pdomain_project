import sys


class Domain:

    def __init__(self, name, c_evalue, start, end,):
        self.name = name
        self.start = start
        self.end = end
        self.score = c_evalue

    def get_motif_format(self):
        return [int(self.start), int(self.end), "()", 100, 10, "black", "rgradient:blue", f"arial|7|white|{self.name}"]

    def __str__(self):
        return f"{self.name}, {self.parent}, {self.start}, {self.end}"
