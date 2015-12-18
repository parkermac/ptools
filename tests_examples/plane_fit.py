"""
Code to test fitting a plane to data.
"""

import numpy as np

x = np.array([0, 1, 0, 1])
y = np.array([0, 0, 1, 1])
z = np.array([1, 1, 3, 3])
A = np.column_stack((np.ones(len(x)), x, y))
coeffs = np.linalg.lstsq(A, z)[0]
print coeffs

# for a line of form z = a + b*x * c*y
# coeffs are [a, b, c]