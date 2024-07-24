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

res = spice.gfdist('EUROPA', 'NONE', 'EUROPA CLIPPER', '<', 500000, 0.0, step_size, nintvls, cnfine)

# Europa's radius
radii = spice.bodvrd('EUROPA', 'RADII', 3)[1] # x, y, z
img = plt.imread("europa.jpg")

good_flybys = [2, 17, 19, 21, 23, 27, 29, 31, 32, 40]

####### GROUND TRACKS < 100,000 km #######
plotdir = 'reason-ranging/'
if spice.wncard(res) == 0:
    print("No flybys found.")
else:
    for i in range(1): #spice.wncard(res)
        print(i+1)
        # Get start + end time for plot labels
        flyby_time = spice.wnfetd(res, i) # whole flyby time this time
        enter = spice.timout(flyby_time[0], "MON DD, YYYY HR:MN ::TDB")
        exit = spice.timout(flyby_time[1], "MON DD, YYYY HR:MN ::TDB")

        # For each time (spaced apart by x points), gather lat/long
        #times = np.linspace(flyby_time[0], flyby_time[1], 700)

        state, _ = spice.spkezr('EUROPA CLIPPER', flyby_time[0], 'J2000', 'NONE', 'EUROPA')
        
        print("##### Flyby " + str(i+1) + " #####")
        print("START: " + enter)
        print("END: " + exit)
        print("State vector: ", state)

        times = np.linspace(flyby_time[0], flyby_time[1], 1000)
        v_total = np.array([])
        alt_total = np.array([])
        for t in times:
            state, _ = spice.spkezr('EUROPA CLIPPER', t, 'J2000', 'NONE', 'EUROPA')
            v = np.linalg.norm([state[3:6]])
            v_total = np.append(v_total, v)
            alt = np.linalg.norm([np.abs(state[0:3]) - radii])
            alt_total = np.append(alt_total, alt)
        
        
        # SAVE FILES in ground_tracks directory
        filename = "velo{}.pdf".format(i+1)
        fig, ax = plt.subplots(figsize=(8,6))
        axs = ax.scatter(alt_total, v_total, c=times)
        clb = plt.colorbar(axs, ax=ax)
        clb.set_label('Time (TDB seconds past J2000)', rotation=270, labelpad=15)

        # prettify the axes
        import matplotlib.ticker as mticker
        ax.set_title("Flyby {}: {} to {}, \nwithin 500,000 km".format(str(i+1), enter, exit), fontsize = 12)
        ax.set_xlabel(r'Altitude (km)', fontsize = 12)
        ax.set_ylabel(r'Velocity (km/s)', fontsize = 12)
        ax.xaxis.set_major_formatter(mticker.FormatStrFormatter('%.2e'))
        ax.locator_params(axis='x', nbins=6)

        fig.savefig(filename, bbox_inches='tight', dpi=600)
        plt.show()
        

        # PLOTS for doppler stuff
        #### constants ####
        c = 2.998e8 # m/s
        hf_pulse_length = 236e-6 
        vhf_pulse_length = 200e-6
        hf_freq = 9e6 #Hz
        vhf_freq = 60e6
        hf_wavelength = c/hf_freq
        vhf_wavelength = c/vhf_freq

        #### VHF ####
        vhf_vel_res = np.array([])
        for idx, d in enumerate(alt_total):
            Niterator = round((2*d*1000/c)/vhf_pulse_length) # divide by pulse length
            sampling_length = Niterator*vhf_pulse_length # s -- multiply by total observed time
            f_d = 2*vhf_freq*v_total[idx]*1000/c
            nyquist_satisfied = bool(3000 >= 2*f_d)
            #print("Nyquist satisfied (PRF >= 2x Doppler freq.)? -- " + str(nyquist_satisfied))
            doppler_res = 1/sampling_length # Hz
            v_res = doppler_res*vhf_wavelength/2 #m/s
            vhf_vel_res = np.append(vhf_vel_res, v_res)
            
        #### HF ####
        hf_vel_res = np.array([])
        for idx, d in enumerate(alt_total):
            Niterator = round((2*d*1000/c)/hf_pulse_length) # divide by pulse length
            sampling_length = Niterator*hf_pulse_length # s -- multiply by total observed time
            f_d = 2*hf_freq*v_total[idx]*1000/c
            nyquist_satisfied = bool(3000 >= 2*f_d)
            doppler_res = 1/sampling_length # Hz
            v_res = doppler_res*hf_wavelength/2 #m/s
            hf_vel_res = np.append(hf_vel_res, v_res)

        #### PLOT RESULTS ####
        fig, ax = plt.subplots(1, 2, figsize=(12,4))
        ax[0].scatter(alt_total, hf_vel_res, c=v_total)
        axs = ax[1].scatter(alt_total, vhf_vel_res, c=v_total)
        ax[0].set_title("HF", fontsize=12)
        ax[1].set_title("VHF", fontsize=12)

        clb = plt.colorbar(axs, ax=ax)
        clb.set_label('Spacecraft velocity (km/s)', rotation=270, labelpad=15)

        for axes in ax.reshape(-1):
            axes.set_xlabel('Altitude (km)', fontsize=12)
            axes.set_ylabel('Velocity resolution (m/s)', fontsize=12)
            axes.grid(True)
            axes.locator_params(axis='x', nbins=7)
            axes.locator_params(axis='y', nbins=7)
            #axes.xaxis.set_major_formatter(mticker.FormatStrFormatter('%.2e'))
            #axes.set_ylim((0, 5))
        #fig.tight_layout()
        fig.savefig("velres{}.pdf".format(i+1), bbox_inches='tight', dpi=600)
        plt.show()

        #plt.close()
            
spice.kclear()