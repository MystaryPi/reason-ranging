# Long Distance Ranging and Velocity Measurements by REASON: Ephemerides Refinement for Europa Clipper
The ephemerides are primarily constrained by radio science measurements on Earth. Here, we present
how REASON can supplement ranging and velocity measurements during the non-nadir phases of flybys.

- We begin by determining the **maximum altitude** a surface echo can be detected from Europa’s surface when pinged by REASON.
- We provide a list of **best working scenarios** that enhance the detection into the non-nadir phases of Clipper’s flybys, especially a coherent train of compressed pulses.
- We then identify the unique and high-priority **flybys** for the long-distance ranging objective.
- Next, we consider a **phase shift approach**, which uses existing ranging pulses for velocity detection. We then present our anticipated range and velocity resolution, from evaluating range and velocity recovery from iteratively sweeping guesses.

## Jupyter Notebooks
- _Long Distance Ranging... .ipynb_ -- Ranging detection + best working scenarios
- _Phase Shift Approach.ipynb_ - Existing ranging pulses for velocity detection (and range/velocity resolution)
- _best-flybys.py_ -- Identifies best flybys from best working scenarios, specifically if spacecraft points to most reflective region within non-nadir phase
- _jovian-bursts.py_ -- Plots Central Meridian Longitude vs. Io Phase to detect if flybys pass through regions of potential Jovian noise
- _orbit-plot-example.py_ -- Example of plotting Europa Clipper flybys with SPICE

SPICE kernels should be located in a "SPICE-Kernels" directory outside of the repository. 
