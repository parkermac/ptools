Notes on the pgrid code.  10/24/2016 Parker MacCready

First edit the gridname and a few directory location at the top of gfun.py.  This gridname will then be used by all subsequent code.

Each time you run a piece of code it makes a new grid file with the name altered to indicate what happened.  The names start as: grid_m00_r00_s00_x00.nc with the letters and numbers indicating changes to: mask, river, smoothing, or extras.

You can use plot_grid.py to look at any of the grids.

Suggested order to run the code:

* make_grid.py

* make_mask.py

* edit_mask.py to get rid of obvious issues like Lake Washington and river channels

* carve_rivers.py

* edit_mask.py for real this time

* smooth_grid.py

* carve_rivers.py again to make sure we did not edit them away

* make_extras.py

* grid_to_LiveOcean.py

