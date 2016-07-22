Notes on the pgrid code.  7/21/2016 Parker MacCready

First edit the gridname and a few directory location at the top of gfun.py.  This gridname will then be used by all subsequent code.

Each time you run a piece of code it makes a new grid file with the name altered to indicate what happened.  The names start as: grid_m00_r00_s00_x00.nc with the letters and numbers indicating changes to: mask, river, smoothing, or extras.

You can use plot_grid.py to look at any of the grids.

Suggested order to run the code:

1. make_grid.py

2. make_mask.py

3. carve_rivers.py

4. edit_mask.py

5. smooth_grid.py

6. make_extras.py

7. grid_to_LiveOcean.py

