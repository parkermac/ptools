"""
Tests of scope when using modules.

Useage: run from the linux command line as

python module_test.py

RESULT 1: Even though sys is imported by module_for_test, it is NOT available
to the calling function.

RESULT 2: However, the addition to the path made in module_for_test IS
available to the calling function.

RESULT 3: Changing a dict entry in a module changes it in the calling code.
    Unexpected!

RESULT 4: Changing a varible in a module does NOT it in the calling code.
    Expected

RESULT 5: A variable defined in a module is available to all the functions
    in the module even if it is not directly passed to the functions.
"""

from importlib import reload
import module_for_test as mft; reload(mft)

mft.add_to_path()

# Test if a module imported in a module is available everywhere.
try:
    print('\nTEST 1:\n' + sys.path[-1])
except NameError:
    print('\nTEST 1:\n' + 'sys is not available')
    print('sys was imported in the module but not in the calling code')

# Test if a path added in a module is available everywhere.
import sys
print('\nTEST 2:\n' + sys.path[-1])
print('path to alpha was added by module, not calling code')

# Test scope of dict
print('\nTEST 3: test of dict scope')
in_dict = {'a': 'original'}
print(in_dict)
mft.change_dict(in_dict)
print(in_dict)
print('dict changed even though is was not returned!')

# Test sope of variables
print('\nTEST 4: test of variable scope')
x = 2
print('Original')
print(' - x in main = ' + str(x))
print('After change_var_noreturn')
mft.change_var_noreturn(x)
print(' - x in main = ' + str(x))
print('After change_var_return')
x = mft.change_var_return(x)
print(' - x in main = ' + str(x))
print('Test of using a variable defined in a module')
try:
    print('a before call = ' + str(a))
except NameError:
    print('a not found before call')
a = mft.look_for_var()
print('a after call = ' + str(a))