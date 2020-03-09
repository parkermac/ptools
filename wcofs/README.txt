This code is designed to do mooring extractions for Alex Kurapov's WCOFS model runs, as part of the COMT 3 project led by Chris Edwards.  Because the model is on a rotated grid and has a different way of storing the output I have had to modify many of the methods used in LiveOcean/x_moor.

===========================================================================================
* mooring_extractor.py does the extraction, working over all the saved daily average files (in LiveOcean_roms/output/) and using mooring locations stored in LiveOcean/x_moor/moor_lists.py under job_name == 'comt3_2014_offshore'.

Input:
- ROMS NetCDF average files in LiveOcean_roms/output/WCOFS_avg_Exp37/zuvts_ORWA_Parker_Exp37_0001.nc which have just basic fields like salt, temp, u, v, and zeta.  These are hand-made cutouts of the larger WCOFS grid from Alex.
- the full grid file (with lon_rho, etc.) is in LiveOcean_roms/output/WCOFS_avg_Exp37/grd_wcofs_large_visc200.nc and you have to look into the mooring_extractor.py code to see the indices that relate the cut-outs to the full grid.

Output:
- LiveOcean_output/moor/WCOFS_avg_Exp37/CA015.nc and etc., one for each mooring.  These are identical in format to those made by LiveOcean/x_moor/mooring_extractor.py and thus are easily viewed using LiveOcean/x_moor/plot_mooring.py.

===========================================================================================
* wcofs_fun.py is a module with a few functions, especially wcofs_lonlat_2_xy() that takes regular lon, lat and returns x, y, the values on the WCOFS rotated grid system.  Since the rotated grid is plaid we then use x, y to do our interpolations following the methods in x_moor.

Note: some of the moorings have masked points, so there is a kludge in wcofs_fun.get_itp_dict() which pushes the interpolant pair to the "west" (in the WCOFS rotated grid), which can result in the rho, u, or v coming from slightly different locations for the same mooring.  Here is the screen output from that kludge:

Warning - masked points for CE015: grid = u
 nudge to the west # 1
 Warning - masked points for MB015: grid = rho
 nudge to the west # 1
 Warning - masked points for MB015: grid = u
 nudge to the west # 1
 Warning - masked points for MB015: grid = v
 nudge to the west # 1
 Warning - masked points for TH015: grid = rho
 nudge to the west # 1
 Warning - masked points for TH015: grid = u
 nudge to the west # 1
 Warning - masked points for TH015: grid = u
 nudge to the west # 2
 Warning - masked points for TH015: grid = v
 nudge to the west # 1
 
showing that MB015 was uniformly shifted west by one point for all grids, but the shifting was inconsistent for CE015 and TH015.  I do not think any of these shallow moorings had current measurements so there should not be a problem.