#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 15:57:23 2017

@author: PM5

Testing function calls with multiple, variable returns.
"""

def test_fun(only_a=False, only_b=False):
    
    def make_a():
        return 'a'
    def make_b():
        return 'b'
    
    if only_a:
        return make_a()
    elif only_b:
        return make_b()
    else:
        return make_a(), make_b()

print('\nTest multiple returned values')
a, b = test_fun()
print('a = %s' % a)
print('b = %s' % b)

print('\nTest single returned value')
aa = test_fun(only_a=True)
print('a = %s' % aa)

print('\nTest other returned value')
bb = test_fun(only_b=True)
print('b = %s' % bb)

print('\nTest raw returned tuple')
c = test_fun()
print('c = %s' % str(c))

print('\nTest to try to break the call')
d = test_fun(only_a=True, only_b=True)
print('d = %s' % str(d))


