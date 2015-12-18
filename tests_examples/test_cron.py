"""
Code to test calling python from cron.

Result: it works! Using a crom entry such as:
(need to specify which python because otherwise the default is 2.7)
(can enter optional arguments in any order)

37 15 * * * /Users/PM5/anaconda3/bin/python
    /Users/PM5/Documents/ptools/tests_examples/test_cron.py -b yo -a hi
    > /Users/PM5/Desktop/test_cron.txt
    
Gives:

a = hi
b = yo

"""
import argparse
import sys

a = sys.path
# for item in a:
#     print(item)
    
# optional command line arguments, can be input in any order
a = 'I am a'
b = 'I am b'
parser = argparse.ArgumentParser()
parser.add_argument('-a', '--arg_a', nargs='?', const=a, type=str, default=a)
parser.add_argument('-b', '--arg_b', nargs='?', const=b, type=str, default=b)
args = parser.parse_args()

print('a = ' + args.arg_a)
print('b = ' + args.arg_b)