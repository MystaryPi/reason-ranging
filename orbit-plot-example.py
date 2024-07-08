'''
Importing SPICE kernels + plotting flybys of Europa Clipper, from the following example:
https://spiceypy.readthedocs.io/en/stable/exampleone.html
'''

import numpy as np
import matplotlib.pyplot as plt; plt.interactive(True)
from mpl_toolkits.mplot3d import Axes3D

import spiceypy as spice

# Print out the toolkit version
spice.tkvrsn("TOOLKIT")

# Import the metakernel
spice.furnsh("/Users/michpark/Sync/Documents/JPL-EUROPA/SPICE-Kernels/21F31v6.tm")

step = 4000
# we are going to get positions between these two dates
utc = ['Apr 1, 2030', 'Aug 31, 2034']

# get et values one and two, we could vectorize str2et
etOne = spice.str2et(utc[0])
etTwo = spice.str2et(utc[1])
print("ET One: {}, ET Two: {}".format(etOne, etTwo))

# get times
times = [x*(etTwo-etOne)/step + etOne for x in range(step)]

positions, lightTimes = spice.spkpos('EUROPA CLIPPER', times, 'J2000', 'NONE', 'EUROPA')
positions = positions.T # positions is shaped (4000, 3), let's transpose to (3, 4000) for easier indexing
fig = plt.figure(figsize=(7, 7))
ax  = fig.add_subplot(111, projection='3d')
ax.plot(positions[0], positions[1], positions[2])
plt.title('SpiceyPy Europa Clipper Position Example from Apr 1, 2030 to Aug 30, 2034')
plt.show()