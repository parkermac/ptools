"""
Test of calling a function that expects arguments.

This is just setting up a python program that accepts command line arguments.

The actual test is done by test_args_call.py.
"""
# get command line arguments
import argparse
parser = argparse.ArgumentParser()
# optional input arguments
parser.add_argument('-t', '--test', nargs='?', const='hi', type=str, default='hi')
args = parser.parse_args()
print(args.test)
