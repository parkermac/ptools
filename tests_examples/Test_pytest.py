"""
This is code to try out the program "pytest", a linux command line tool that
comes with anaconda python.  The motivation came from reading this:
http://katyhuff.github.io/python-testing/

To test this from the linux command line, type:

pytest Test_pytest.py

or, for verbose mode:

pytest Test_pytest.py -v 

Since I already have many pieces of code titled "test_" that are not designed
to be unit tests of module functions (like pytest expects) I will reserve
the name "Test_" for all modules that are designed to be used with pytest.

Note that the functions defined in the test module should all start with "test_"
in order for pytest to look at them.

"""

import numpy as np

def test_mean():
    exp = 2
    obs = np.arange(5).mean()
    assert obs == exp