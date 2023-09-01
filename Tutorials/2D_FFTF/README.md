# 2D FFTF

Author: Stefan Radman  
Review and editing: Carlo Fiorina

## Info

This is a model of the Fast Flux Test Facility, a former Sodium Fast
Reactor at the Hanford Site, WA, US, operated in the 1980s. It was in a hybrid
pool-loop configurations, where the the primary pumps and IHXs lie in a loop
outside the vessel.

This case represents the LOFWOS 13 Test performed at the FFTF in order to test
the effectiveness of the Gas Expansion Module (GEM) safety features. In
essence, these were empty assembly wrappers closed at the top and open at the
bottom, partly filled with argon gas, positioned on the periphery of the active
core. In operation, the pressure head at the GEM inlet would compress the
argon gas so that the free surface sodium level would rise above the active
core level. During a ULOF accident, the loss of pressure head would cause the
argon gas to expand and lower the free surface sodium level below the active
core region. From a neutronics perspective, this "uncovers" part of the core
radially, as were there was sodium, now there is argon gas, thus increasing
radial leakage. Needless to say, the FFTF core itself was rather small, ~ 1 m
in diameter, which explains the effectiveness of said feature.

In essence, the LOFWOS 13 test was:

- operate at steady state, at full flow and half the total power;
- trip the pumps;
- watch the GEMs do their magic;


## Calculation details

The model is hydrid as the vessel consists of a 2-degree wedge while the loops
consist of paralellopipes. All the absolute volumes are scaled by a factor
360/2 of the total FFTF primary volume.

From a calculation perspective, the steady state is run for 900 s of model time,
followed by a transient case in which a pointKinetics model is used to model the
the power evolution. Please notice that the power for the steady state is set
directly via the field powerDensity.nuclearFuelPin in 0/fluidRegion. As an
alternative, the `Allrun_powerFromDiffusion` script is provided that runs an
initial diffusion calculations and uses the resulting power profile in the
transient.

The case can run on 5 cores via the `Allrun_parallel` script (suggested) or on
1 core using the `Allrun` script.
On 5 cores of an Intel i5-8600K CPU @ 3.60GHz the calculation durations were:
-	~ 820 s for the steady state
- 	~ 3700 s for the transient

The transient starts at t = 910 s. For convenience, the transient simulation
ends at 1200 s, as most of the dynamics is resolved at that point (even though
a new steady state is not reached). Change the endTimeT variable in the `Allrun`
script to run up to any point you desire.


## GeN-Foam features

From the perspective of GeN-Foam users, this case showcases the usage of two
new features, namely the gapHPowerDensityTable and a time-dependent
momentumSource. Both these features are discussed in
`constant/fluidRegion/phaseProperties`. A third feature is the GEM, which are
currently partially hard coded in the pointKinetics class. I say "partial"
because the GEM Reactivity Map (i.e. reactivity VS free surface sodium height)
is provided manually in `constant/neutroRegion/nuclearData`, as well as the
name of the faceZone over which total flow is kept track of. The hard-coded
part consist of the dependence of the free surface sodium height with respect
to the total integral volume flow through said faceZone (correlation provided
by the experimentalists). Thus, this feature is fairly useless for any other
application for now. Apart from that, the case showcases the use of the novel
decay power model point kinetics (not implemented for liquidFuel yet, though).
In `constant/neutroRegion/reactorState` one can optionally define a decay power
and its evolution in time (via a table or expressions evaluated via the
standard OpenFOAM Function1 feature). The power provided under pTarget is
the initial TOTAL power, and the pointKinetics model will operate only on
the fission power, initially set as total-decay.
The heatExchanger feature is also used to model heat
transfer between the two physically separated mesh-domains representing the
primary and the secondary. Furthermore, the model uses the novel momentumMode
name cellCenteredFaceReconstruction, which in very simple terms is a mix of
the cellCentered mode (i.e. not too diffusive) the and the faceCentered mode
(i.e. not too oscillatory at porous interfaces between regions with different
porosities). The 'porousInterfaceSharpness' factor controls the smoothness of
the velocity discontinuity at said interfaces. In 1-D scenarios, setting it to
1 completely eliminates any form of oscillations at the interface, and the
jumps are as sharp as they should numerically be. However, in 2-D or 3-D
scenarios, a large value of the porousInterfaceSharpness will results
in numerical problems whose nature is not still understood (by me). If the
coefficient is set to 0, a smooth transition will occur, without oscillations.
However, this is not in line with what a sharp porosity changes requires: a
sharp change in velocity. Thus, intermediate values of the coefficient, between
0 and 1, blend between the sharp (potentially unstable) velocity transitions
a smooth (yet oscillation-free) transition. This case uses a
porousInterfaceSharpness = 0.5 (see fvSolution)


## How to run

To run everything, simply launch the `./Allrun` or `./Allrun_parallel` script.
To plot results, use the plot script by passing the log name as an argument,
i.e. either

```bash
python3 plot.py steadyState/log
# or
python3 plot.py transient/log
```


## Notes on the plots

The plotted results consists in the evolution of total power, mass flows for a
variety of faceZones, bulk temperatures for a variety of faceZones and
reactivity contributions of different components in dollars. Experimental
results are marked via 'exp'.


## Limitation of the model

The FFTF core was largely radially thermo-hydraulically heterogeneous, meaning
that neighboring assemblies could have significantly different power-to-flow
ratios, resulting in different outlet temperatures. By representing the core
only via two averaged regions, namely innerCore and outerCore, this
heterogeneity is lost. This can be observed in the comparison of the
experimental temperature (i.e. PIOTA2) against the innerCore temperature. The
PIOTA2 assembly, which belongs to the inner core region of the FFTF, hard
a higher power-to-flow ratio than the region average, so that results do not
compare perfectly in magnitude.

Accurate predictions of the temperature far into the transient are complicated
by the fact that both the primary power and flow are fractions of what they
used to be at steady-state, so that even the slightest deviation form the
experimental power-to-flow ratio will result in potentially large
temperature differences in different regions of the model.

The model is preliminary as many parameters that have a strong impact on
transient evolution (e.g. the gap conductance VS fuel linear power map) have
not been thoroughly investigated, yet the overall trends can be reproduced
reasonably, especially for what concerns the initial minute of the transient,
in which the power evolution is shaped almost exclusively dominated by the
Doppler, fuel axial expansion and GEM reactivity contributions.

Further improvements to the model will be made.
