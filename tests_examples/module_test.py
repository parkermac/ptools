"""
Test of calling modules.

RESULT 1: Even though sys is imported by module_for_test, it is NOT available
to the calling function.

RESULT 2: However, the addition to the path made in module_for_test IS
available to the calling function.

RESULT 3: in module_for_test we ARE able to import modules at the start
that are then available to all functions defined in the module.
"""
import module_for_test as mft; reload(mft)

mft.add_to_path()

try:
    print('\nTEST 1:\n' + sys.path[-1])
except NameError:
    print('\nTEST 1:\n' + 'sys is not available')
    
import sys
print('\nTEST 2:\n' + sys.path[-1])