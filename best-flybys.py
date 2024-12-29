'''
Importing SPICE kernels + plotting flybys of Europa Clipper, from the following example:
https://spiceypy.readthedocs.io/en/stable/exampleone.html
'''

import numpy as np
import matplotlib.pyplot as plt; plt.interactive(True)
import matplotlib.colors as colors
from matplotlib import cm
import spiceypy as spice
import scienceplots
plt.interactive(True)

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

# Europa's radius
radii = spice.bodvrd('EUROPA', 'RADII', 3)[1] # x, y, z
img = plt.imread("europa.jpg")

good_flybys = [1,43,23] # i -1  [1,16,18]

# SAVE FILES in ground_tracks directory
plt.style.use(['science', 'no-latex'])
filename = "figs/flybys-paper.pdf"
fig, ax = plt.subplots(3,1,figsize=(8,7))

####### GROUND TRACKS < 100,000 km #######
plotdir = '/Users/michpark/Sync/Documents/JPL-EUROPA/Slide Figures/'
if spice.wncard(res) == 0:
    print("No flybys found.")
else:
    index = 0
    for i in good_flybys: #spice.wncard(res)
        print(i+1)
        # Get start + end time for plot labels
        flyby_time = spice.wnfetd(res, i) # whole flyby time this time
        enter = spice.timout(flyby_time[0], "MON DD, YYYY HR:MN ::TDB")
        exitt = spice.timout(flyby_time[1], "MON DD, YYYY HR:MN ::TDB")

        et_ca = ((flyby_time[1] - flyby_time[0]) / 2) + flyby_time[0]
        # finding closest approach

        # For each time (spaced apart by x points), gather lat/long
        #times = np.linspace(flyby_time[0], flyby_time[1], 700)
        times = np.linspace(et_ca - 28*3600, et_ca + 28*3600, 56*3600)
        lat_total = np.array([])
        long_total = np.array([])
        alt_total = np.array([])
        v_total = np.array([])

        for t in times:
            state, _ = spice.spkezr('EUROPA CLIPPER', t, 'IAU_EUROPA', 'NONE', 'EUROPA')
            intersect, _,  alt = spice.subpnt('NEAR POINT/ELLIPSOID', 'EUROPA', t, 'IAU_EUROPA', 'NONE', 'EUROPA CLIPPER')
            r, lon, lat = spice.reclat(intersect) # rect. coords (radians)
            lon_deg = spice.dpr() * lon # degrees
            lat_deg = spice.dpr() * lat

            # gather data for each flyby
            lat_total = np.append(lat_total, lat_deg)
            long_total = np.append(long_total, 360-((lon_deg+360)%360)) #convert to 0-360 W scale
            alt_total = np.append(alt_total, np.sqrt(np.sum((alt-radii)**2)))

            # Find the RADIAL VELOCITY
            #v_total = np.append(v_total, spice.vnorm(state[3:]))
            pos = state[0:3]/spice.vnorm(state[0:3])
            radial_v = spice.vdot(state[3:],pos)
            #trans_v = np.sqrt(spice.vnorm(state[3:])**2-radial_v**2)
            v_total = np.append(v_total, radial_v) # RADIAL VELOCITY

        '''
        plt.figure()
        t = np.linspace(-28, 28, 56*3600)
        plt.plot(t, v_total)
        plt.xlabel("Hours from closest approach")
        plt.ylabel("Radial Velocity (km/s)")
        
        plt.figure()
        t = np.linspace(-28, 28, 56*3600-1)
        plt.plot(t, np.diff(v_total))
        plt.xlabel("Hours from closest approach")
        plt.ylabel("Difference in Radial Velocity (km/s^2)")
        '''

        # fundamental limit
        # Given the velocity changes you calculated, at what altitude is 
        # the time by which the measurement becomes inaccurate the same as the round trip time? 
        # Determine change in velocity (for limitation on max. # pulses)
        #t = np.linspace(-28, 28, 700)
        '''
        v_ch = 0.0001e3 # km/s^2
        # assume within 28 hs
        time_inacc = 0.05/v_ch # using velocity res = 0.05 m/s
        print("{} seconds until velocity measurement becomes inaccurate".format(time_inacc))
        alt_limit = time_inacc*2.998e8/2
        print("Altitude limit: {} km".format(alt_limit/1e3))

        print("##### Flyby " + str(i+1) + " #####")
        print("START: " + enter)
        print("END: " + exit)
        #print("Altitude: ", alt_total)
        #print("Velocities: ", v_total)
        '''
        
        ax[index].set_title("Flyby {}: {} - {}".format(str(i+1), enter, exitt[-5:]),fontsize=14.5)
        ax[index].set_xlim(360, 0)
        ax[index].set_ylim(-57, 57)
        ax[2].set_xlabel(f'Longitude [$\degree$W]',fontsize=14.5)
        ax[index].set_ylabel(f'Latitude [$\degree$]',fontsize=14.5)
        ax[index].tick_params(axis='both', which='major', labelsize=14)

        # leading subjovian = 270-360 E == 90-0 W
        import matplotlib.patches as patches
        # Create a Rectangle patch
        rect = patches.Rectangle((0, -57), 90, 114, linewidth=4, edgecolor='#314173', facecolor='none')
        ax[index].add_patch(rect) 

        # colorbar with altitude
        # plot flybys, with colorbar for altitude
        axs = ax[index].scatter(long_total, lat_total, s=13, c=alt_total, cmap='inferno', vmin=0, vmax=100000)
        axp = ax[index].imshow(img, extent=[360, 0, -57, 57])
        clb = fig.colorbar(axs, ax=ax[index],shrink=0.9,pad=0.01,aspect=13) #, shrink=0.5, pad=0.01, aspect=13
        clb.set_label('Altitude (km)', rotation=270, labelpad=15,size=14)
        index+=1
        
        fig.tight_layout()
        plt.show()
        fig.savefig(filename)
    

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