'''
Importing SPICE kernels + plotting flybys of Europa Clipper, from the following example:
https://spiceypy.readthedocs.io/en/stable/exampleone.html
'''

import numpy as np
import matplotlib.pyplot as plt; plt.interactive(True)
import matplotlib.colors as colors
from matplotlib import cm
import spiceypy as spice

# Import the metakernel
spice.furnsh("/Users/michpark/Sync/Documents/JPL-EUROPA/SPICE-Kernels/21F31v6.tm")

####### FIND FLYBYS #######
# Between start + end time
start_time = spice.str2et('Apr 1, 2030')
end_time = spice.str2et('Aug 31, 2034') # Aug 31, 2034
print("Start time: {:e}, end time: {:e}".format(start_time, end_time))

# Define parameters
step_size = 300 # s
nintvls = int(2*2 + (end_time-start_time)/step_size) # number of intervals to return = 2*n + m/step
cnfine = spice.cell_double(2)
spice.wninsd(start_time, end_time, cnfine) # set time window to search
res = spice.cell_double(2 * nintvls) # Spice cell to store output
# Set LOCMIN - (in distance), refval = adjust = 0.0 in this case

res = spice.gfdist('EUROPA', 'NONE', 'EUROPA CLIPPER', '<', 10000, 0.0, step_size, nintvls, cnfine)

# Europa's radius
radii = spice.bodvrd('EUROPA', 'RADII', 3)[1] # x, y, z
img = plt.imread("europa.jpg")

good_flybys = [2, 17, 19, 21, 23, 27, 29, 31, 32, 40]

####### GROUND TRACKS < 100,000 km #######
plotdir = 'ground_tracks/'
hide_ground_tracks = False
if spice.wncard(res) == 0 or hide_ground_tracks:
    print("No flybys found.")
else:
    for i in range(1): #spice.wncard(res)
        print(i+1)
        # Get start + end time for plot labels
        flyby_time = spice.wnfetd(res, i) # whole flyby time this time
        enter = spice.timout(flyby_time[0], "MON DD, YYYY HR:MN ::TDB")
        exitt = spice.timout(flyby_time[1], "MON DD, YYYY HR:MN ::TDB")

        # For each time (spaced apart by x points), gather lat/long
        times = np.linspace(flyby_time[0], flyby_time[1], 700)
        lat_total = np.array([])
        long_total = np.array([])
        alt_total = np.array([])
        cml_total = np.array([])

        for t in times:
            state, _ = spice.spkezr('EUROPA CLIPPER', t, 'J2000', 'NONE', 'EUROPA')
            intersect, _,  alt = spice.subpnt('NEAR POINT/ELLIPSOID', 'EUROPA', t, 'IAU_EUROPA', 'NONE', 'EUROPA CLIPPER')
            r, lon, lat = spice.reclat(intersect) # rect. coords (radians)
            lon_deg = spice.dpr() * lon # degrees
            lat_deg = spice.dpr() * lat

            # gather data for each flyby
            lat_total = np.append(lat_total, lat_deg)
            long_total = np.append(long_total, 360-((lon_deg+360)%360)) #convert to 0-360 W scale
            alt_total = np.append(alt_total, np.sqrt(np.sum((alt-radii)**2)))

            # Find the Central Meridian Longitude with Jupiter
            intersect_jupiter, _, _ = spice.subpnt('NEAR POINT/ELLIPSOID', 'JUPITER', t, 'IAU_JUPITER', 'NONE', 'EUROPA CLIPPER')
            _, cml, _ = spice.reclat(intersect_jupiter)
            cml = spice.dpr() * cml # degrees
            cml_total = np.append(cml_total, cml)

        


spice.kclear()