# Godiva SN

This is a purely neutronic eigenvalue case displaying how to use the discrete
ordinate (SN) solver of GeN-Foam. It simulates the Godiva experiment,
constituted by a small super-prompt-critical sphere of enriched uranium. The SN
solver is selected in the constant/neutroRegion/neutronicProperties dictionary.
The file constant/neutroRegion/quadratureSet contains a simple quadrature set
with 4 direction per octant. A more complex  (and more computationally
requiring) quadrature set with 16 directions per octant can be found in
constant/neutroRegion/quadratureSet16. A simpler one, with 1 direction per
octant, can be found in constant/neutroRegion/quadratureSet16. The scattering
anistotropy can be changed by changing the legendreMoments flag in
constant/neutroRegion/nuclearData. Of course one should make sure that the
corresponding scattering matrices are provided in the same file. In the current
tutorial, scattering matrix are provided till the 5th moment. This is a
computationally intensive tutorial. It is suggested to run it using the
`Allrun_parallel` bash script on a good computer. In principle, the SN solver
could be used for time-dependent calculations. However, no acceleration
techniques are currently implemented, making the solver particularly slow. An
`Allclean` script is provided to clean up the case.
