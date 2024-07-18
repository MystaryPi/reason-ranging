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
<<<<<<< HEAD
    for i in range(1): #spice.wncard(res)
        print(i+1)
        # Get start + end time for plot labels
        flyby_time = spice.wnfetd(res, i) # whole flyby time this time
        enter = spice.timout(flyby_time[0], "MON DD, YYYY HR:MN ::TDB")
        exit = spice.timout(flyby_time[1], "MON DD, YYYY HR:MN ::TDB")

        # For each time (spaced apart by x points), gather lat/long
        times = np.linspace(flyby_time[0], flyby_time[1], 700)
        lat_total = np.array([])
        long_total = np.array([])
        alt_total = np.array([])

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
    
        # SAVE FILES in ground_tracks directory
        filename = "flyby_{}.pdf".format(i+1)
        fig, ax = plt.subplots(figsize=(36,30))

        # prettify the axes
        ax.set_title("Flyby {}: {} to {}, within 100,000 km".format(str(i+1), enter, exit), fontsize = 15)
        ax.set_xlim(360, 0)
        ax.set_ylim(-57, 57)
        ax.set_xlabel(r'Longitude [$\degree$W]', fontsize = 12)
        ax.set_ylabel(r'Latitude [$\degree$]', fontsize = 12)

        # leading subjovian = 270-360 E == 90-0 W
        import matplotlib.patches as patches
        # Create a Rectangle patch
        rect = patches.Rectangle((0, -57), 90, 114, linewidth=4, edgecolor='#314173', facecolor='none')
        ax.add_patch(rect) 
        plt.tight_layout()

        # colorbar with altitude
        # plot flybys, with colorbar for altitude
        axs = ax.scatter(long_total, lat_total, s=9, c=alt_total, cmap='inferno', vmin=0, vmax=100000)
        axp = ax.imshow(img, extent=[360, 0, -57, 57])
        clb = plt.colorbar(axs, ax=ax, shrink=0.14, pad=0.01, aspect=30)
        clb.set_label('Altitude (km)', rotation=270, labelpad=15)
        
        fig.savefig(plotdir+filename, bbox_inches='tight', dpi=600)
        plt.close()

=======
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
>>>>>>> 2cf2837852507cdd97ef214e87418da91d0f16c8

####### BEST WORKING SCENARIOS #######
'''
• Within maximum range
• Preferable with HF
• No Jovian noise
• Long pulses
• Pulse train
• Pointing towards ridged plains (leading hemisphere, near anti-Jovian point

With SPICE, we can evaluate when within ___ max range, 
800,000 km with HF for 100 pulses 
'''



spice.kclear()