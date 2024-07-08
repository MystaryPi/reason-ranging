import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import spiceypy as spice

# Print out the toolkit version
spice.tkvrsn("TOOLKIT")

# Import the metakernel
spice.furnsh("SPICE-Kernels/21F31v6.tm")
