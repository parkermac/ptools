#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 10:59:22 2017

@author: PM5

Test of using the main() function, from
http://interactivepython.org/courselib/static/thinkcspy/Functions/mainfunction.html

The concept is that if I run the program it executes main(),
but if I import it in another program with a line like
import example_main, it does not execute main() but does make the functions
available.  Cool.
"""

def squareit(n):
    return n * n

def cubeit(n):
    return n*n*n

def main():
    anum = int(input("Please enter a number "))
    print(squareit(anum))
    print(cubeit(anum))

if __name__ == "__main__":
    main()
else:
    # this is my addition
    print(__name__)
    # and what it did was, for example:
    # In: import example_main as em
    # Out: example_main   