"""
Code to solve a system of basins and straits.
"""
import numpy as np

import bst_fun as bsf
from importlib import reload
reload(bsf)

v_list = ['s'] # set variables to use (always include 's')
case = 'base' # choose basin and strait setup
fcase = 'base' # choose forcing functions
ND = 1000 # set number of days to run

# get the initial setup
Ba, St = bsf.make_BaSt(case, v_list)
# Note: the entries in dicts Ba and St are implicitly global, and are
# modified by functions, even if we do no return them.

# Create output arrays, with columns in dicts "_oc"
Ba_out, Ba_oc, St_out, St_oc = bsf.make_BaSt_out(Ba, St, ND, v_list)

# run to steady state for initial condition
Ba_out, St_out, TD = bsf.bcalc(Ba, St, Ba_out, Ba_oc, St_out, St_oc, v_list, ND, fcase, ic=True)
# Ba and St are also updated, and are the initial condition for the real run.

# run for real
Ba_out, St_out, TD = bsf.bcalc(Ba, St, Ba_out, Ba_oc, St_out, St_oc, v_list, ND, fcase)

# plotting
bsf.plot_basic(Ba, St, Ba_out, Ba_oc, St_out, St_oc, TD)
