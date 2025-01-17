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

good_flybys = [2, 17, 19, 21, 23, 27, 29, 31, 32, 40]

####### IO PHASE #######
def calc_orbphase(obs_v, targ_v):
    """
    function for properly computing orbital phase
    input: Planet to Observer vector (ecliptic coords), planet to Target Vector
    output:  the orbital phase of the target relative to the observer's superior conjunction
    """    
    superior_vector=np.zeros(3)
    #superior conjunction is opposite our observer (i.e. negative vector)
    superior_vector[0]=-obs_v[0]
    superior_vector[1]=-obs_v[1]
    superior_vector[2]=-obs_v[2]    #calculate phase in the x-y plane (slight approximation due to inclinations...)
    sup_phase=np.arctan2(superior_vector[1], superior_vector[0])
    targ_phase=np.arctan2(targ_v[1], targ_v[0])    

    '''
    this value will be between -360.0 and 360.0    
    if phase_diff<0.0: phase_diff=phase_diff+360.0
    restrict to 0 to 360.0 (same as literature)  
    '''
    phase_diff=(targ_phase-sup_phase)*180/np.pi 

    return phase_diff

####### GROUND TRACKS < 100,000 km #######
filename = "figs/jovian-flyby-paper.pdf"
plt.style.use(['science', 'no-latex'])
fig, ax = plt.subplots(1,2,figsize=(8,4))

if spice.wncard(res) == 0:
    print("No flybys found.")
else:
    index = 0
    for i in [1, 23]: #spice.wncard(res)
        print(i+1)
        # Get start + end time for plot labels
        flyby_time = spice.wnfetd(res, i) # whole flyby time this time
        enter = spice.timout(flyby_time[0], "MON DD, YYYY HR:MN ::TDB")
        exitt = spice.timout(flyby_time[1], "MON DD, YYYY HR:MN ::TDB")

        # For each time (spaced apart by x points), gather lat/long
        times = np.linspace(flyby_time[0], flyby_time[1], 700)
        #lat_total = np.array([])
        #long_total = np.array([])
        #alt_total = np.array([])
        cml_total = np.array([])
        io_phase_total = np.array([])

        for t in times:
            state, _ = spice.spkezr('EUROPA CLIPPER', t, 'J2000', 'NONE', 'EUROPA')
            intersect, _,  alt = spice.subpnt('NEAR POINT/ELLIPSOID', 'EUROPA', t, 'IAU_EUROPA', 'NONE', 'EUROPA CLIPPER')
            r, lon, lat = spice.reclat(intersect) # rect. coords (radians)
            lon_deg = spice.dpr() * lon # degrees
            lat_deg = spice.dpr() * lat

            #### gather data for each flyby
            #lat_total = np.append(lat_total, lat_deg)
            #long_total = np.append(long_total, 360-((lon_deg+360)%360)) #convert to 0-360 W scale
            #alt_total = np.append(alt_total, np.sqrt(np.sum((alt-radii)**2)))

            #### Find the Central Meridian Longitude with Jupiter
            intersect_jupiter, _, _ = spice.subpnt('NEAR POINT/ELLIPSOID', 'JUPITER', t, 'IAU_JUPITER', 'NONE', 'EUROPA CLIPPER')
            _, cml, _ = spice.reclat(intersect_jupiter)
            cml = spice.dpr() * cml # degrees
            cml_total = np.append(cml_total, (cml+360)%360)

            ##### Find the Io phase -> angle btwn jupiter, io, clipper
            # Define parameters
            obsrvr = 'EUROPA CLIPPER'
            planet = 'JUPITER'
            frame='IAU_EUROPA'
            target_moon = 'IO'

            # from gregor
            vplanet_obs, lt = spice.spkpos(obsrvr, t, frame, "XLT", planet)
            vplanet_targ, lt = spice.spkpos(target_moon, t, frame,"NONE",planet)
            orbphase0=calc_orbphase(vplanet_obs,vplanet_targ)        
            
            io_phase_total = np.append(io_phase_total, (orbphase0+360)%360)
            
        # plot CML vs Io phase
        ax[index].set_title("Flyby {}: {} - {}".format(str(i+1), enter, exitt[-5:]),fontsize=14.5)
        ax[index].scatter(cml_total, io_phase_total, c='black',s=15)
        
        ax[index].set_xlabel(f'Central meridian longitude [$\degree$]',fontsize=14.5)
        ax[0].set_ylabel(f'Io phase [$\degree$]',fontsize=14.5)

        ax[index].tick_params(axis='both', which='major', labelsize=14)

        import matplotlib.patches as patches
        # create rectangle patches for all of the regions
        rectangles = {'Io-A' : patches.Rectangle((185, 185), 90, 80, linewidth=4, edgecolor='#ac2941', facecolor='#ac2941', alpha=0.5),
                "Io-A'" : patches.Rectangle((200,160), 65, 35, linewidth=4, edgecolor='#de4a6a', facecolor='#de4a6a', alpha=0.5),
                "Io-A''" : patches.Rectangle((300,220), 45, 25, linewidth=4, edgecolor='#bd314a', facecolor='#bd314a', alpha=0.5),
                "Io-B" : patches.Rectangle((75,55), 125, 60, linewidth=4, edgecolor='#205283', facecolor='#205283', alpha=0.5),
                "Io-B'" : patches.Rectangle((130,65), 55, 30, linewidth=4, edgecolor='#83cde6', facecolor='#83cde6', alpha=0.5),
                "Io-C" : patches.Rectangle((270,215), 90, 45, linewidth=4, edgecolor='#ee5a31', facecolor='#ee5a31', alpha=0.5),
                "" : patches.Rectangle((0,215), 45, 45, linewidth=4, edgecolor='#ee5a31', facecolor='#ee5a31', alpha=0.5),
                "Io-D" : patches.Rectangle((20, 95), 210, 30, linewidth=4, edgecolor='#624a73', facecolor='#624a73', alpha=0.5)}

        for r in rectangles:
            ax[index].add_artist(rectangles[r])    
            rx, ry = rectangles[r].get_xy()
            cx = rx + rectangles[r].get_width()/2.0
            cy = ry + rectangles[r].get_height()/2.0

            ax[index].annotate(r, (cx, cy), color='w', weight='bold', 
                        fontsize=12, ha='center', va='center',usetex=False)
        ax[index].set_xlim((0, 360))
        ax[index].set_ylim((0, 360))
        ax[index].grid(True, alpha=0.2)
        index += 1

fig.tight_layout()
plt.show()
fig.savefig(filename)
#plt.close()



spice.kclear()