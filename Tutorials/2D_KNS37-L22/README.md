# 2D KNS37-L22

This test-case is representative of the test L22 sodium boiling transient that
was performed at the KNS-37 facility at Karlsruhe in the 1980s. It involves a
single mock-up assembly with 37 electrically heated pins undergoing a ULOF,
with subsequent sodium boiling, two-phase flow excursion, power-off, and flow
restoration. For further details, refer to:

F. Huber and A. Kaiser and K. Mattes and W. Peppler, "Steady state and
transient sodium boiling experiments in a 37-pin bundle", Nuclear Engineering
and Design vol 100, pp. 377-386, 1987,
https://www.sciencedirect.com/science/article/pii/0029549387900872

The computational models consists of a 2-D axial-symmetric wedge representative
of the assembly, while the pump ramp-down is represented via a time-dependent
inlet pressure boundary condition.

The simulation lasts 20s in the -2s to 18s range. The first 5s of simulation
serve the purpose of reaching a steady-state before the pump trip. Pump trip
occurs at 3s, and boiling starts at ~ 9.1s. Flow excursion occurs when the
vapour has radially expanded to occupy the entirety of the available cross-
sectional assembly area, and power is turned off at ~ 12.4s. After that,
the liquid-vapour mixture cools the bundle and condenses until the two-phase
flow ends at ~ 15s.

The plotting script `plot.py` is specific for this transient and allows to
visualize the time evolution of several quantities. To use it, you need to have
matplotlib installed, it works with any Python version >= 2.7. To launch it,
run

``` bash
python plot.py log.GeN-Foam
# or
python3 plot.py log.GeN-Foam
```
