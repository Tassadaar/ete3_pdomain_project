import sys


class Domain:

    def __init__(self, name, parent, start, end):
        self.name = name
        self.parent = parent
        self.start = start
        self.end = end