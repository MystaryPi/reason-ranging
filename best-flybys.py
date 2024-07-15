'''
Importing SPICE kernels + plotting flybys of Europa Clipper, from the following example:
https://spiceypy.readthedocs.io/en/stable/exampleone.html
'''

import numpy as np
import matplotlib.pyplot as plt; plt.interactive(True)
from mpl_toolkits.mplot3d import Axes3D

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

res = spice.gfdist('EUROPA', 'NONE', 'EUROPA CLIPPER', '<', 100000, 0.0, step_size, nintvls, cnfine)

# Print all flybys (in res)
if spice.wncard(res) == 0:
    print("No flybys found.")
else:
    for i in range(spice.wncard(res)):
        flyby_time = spice.wnfetd(res, i)[0]  # Get the time of closest approach
        calendar_string = spice.timout(flyby_time, "MON DD,YYYY  HR:MN:SC.#### (TDB) ::TDB")
        state, lighttime = spice.spkezr('EUROPA CLIPPER', flyby_time, 'J2000', 'NONE', 'EUROPA')
        distance = spice.vnorm(state[:3])


        # Get the transformation matrix from Clipper ->  Europa
        sc_to_body_matrix = spice.pxform('J2000', 'IAU_EUROPA', flyby_time)

        # Define the boresight vector in instrument frame
        # direction vector: along z axis, pointing down [0,0,-1]
        # velocity vector: assume non rotating?
        boresight = [0, 0, 1, 0, 0, 0]

        # Get Europa's radius
        radii = spice.bodvrd('EUROPA', 'RADII', 3)[1]

        # Find the intersection of the boresight with the target body
        intersect_point, found = spice.surfpv(state, boresight, radii[0], radii[1], radii[2])

        # PRINT INFO
        print(f"Flyby {i+1}:")
        print(f"  Time to enter {distance:.0f} km: {calendar_string}")
        print("State: " + str(state)) # position 1st 3, velocity last 3
        print(f"REASON is pointing at:")

        if found:
            # Convert to lat, long in degrees
            r, lon, lat = spice.reclat(intersect_point)
            lon_deg = spice.dpr() * lon
            lat_deg = spice.dpr() * lat

            print(f"Spacecraft is pointing at:")
            print(f"Latitude: {lat_deg:.6f} degrees")
            print(f"Longitude: {lon_deg:.6f} degrees")
        else:
            print("No intersection found. The spacecraft may not be pointing at the target body.")

        print("---")

####### BEST WORKING SCENARIOS #######
'''
• Within maximum range
• Preferable with HF
• No Jovian noise
• Long pulses
• Pulse train
• Pointing towards ridged plains (leading hemisphere, near anti-Jovian point

With SPICE, we can evaluate when within ___ max range, 
'''


spice.kclear()