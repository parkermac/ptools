Notes on the pgrid code.  Parker MacCready

This collection of programs is designed to make gridfiles for ROMS.  It works for analytical and realistic cases, and handles complex tasks like mask editing and smoothing.  It also creates other files associated with rivers, nudging to climatology, and vertical grid parameters, all in the form expected by LiveOcean/forcing and ROMS.

First edit the gridname and a few directory locations at the top of gfun.py.  This gridname will then be used by all subsequent code.

In order to keep track of several choices typically made about a grid, we use "dch," a dict of "default choices":
- dch =  gfun.default_choices(Gr)
These are initialized in gfun.default_choices, but you typically override some of them in make_grid.py when doing the initial grid specification.  The choices are saved in a pickle file:
- pickle.dump(dch, open(Gr['gdir'] + 'choices.p', 'wb'))
You can also go back and change things in dch later using Z_edit_dch.py.

Each time you run a piece of code it makes a new grid file with the name altered to indicate what happened.  The names start as: grid_m00_r00_s00_x00.nc with the letters and numbers indicating changes to: mask, river, smoothing, or extras.

You can use plot_grid.py to look at any of the grids.

Suggested order to run the code, and what dch items are used:

* make_grid.py
    analytical
    do_cell_average
    t_dir/t_list
    use_z_offset/z_offset

* make_mask.py
    z_land
    do_cell_average/z_lan_alt
    unmask_coast
    remove_islands

* edit_mask.py (to get rid of obvious issues like Lake Washington and river channels)

* carve_rivers.py

* edit_mask.py (for real this time, perhaps running many times)

* smooth_grid.py
    use_min_depth/min_depth
    fjord_cliff_edges

* carve_rivers.py (again to make sure we did not edit them away)

* make_extras.py
    min_depth (enforced for whole grid)

* grid_to_LiveOcean.py
    nudging_edges
    nudging_days

