This folder is designed to hold all the modules and files required by code used in the SSMSP (Salish Sea Marine Survival Project) long-term data analysis.  Specifically it is relied upon by:

- obs_ncdc
- river_long
- obs_collias
- obs_ecy

All of the code and files are copied from LiveOcean sources, and ordinarily I would leave it there, but the goal is to make the programs in these folders work on their own.

To use the ssmsp/sfun.py module, you would add these lines to any code in the folders above:

# SSMSP import
import os
import sys
pth = os.path.abspath('../ssmsp')
if pth not in sys.path:
    sys.path.append(pth)
import sfun
from importlib import reload
reload(sfun)

and the same goes for river_class.py and pfun.py.
