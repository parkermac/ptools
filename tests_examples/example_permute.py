"""
permutations - for tic tac toe moves
"""

import itertools
l = list(itertools.permutations(range(4)))
for ll in itertools.permutations(range(4)):
    print ll

