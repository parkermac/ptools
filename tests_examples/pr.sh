#!/bin/bash

cd /Users/PM5/Documents/ptools/tests_examples

# this runs my Enthought-Canopy python 2.7
python ./pr.py > pr0.txt &

# this runs my anaconda python 2.7
/Users/PM5/anaconda/bin/python ./pr.py  > pr1.txt &

# this runs my anaconda python 3.5
/Users/PM5/anaconda3/bin/python ./pr.py  > pr2.txt &

# RESULT: all work when this script is run from TextMate
# and (only!) after I added the "> pr#.txt &"
# they also worked from from the command line.

