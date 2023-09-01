# 1D HX

## Description

These test cases showcase the utilization of the novel `heatExchanger` feature.
This feature allows to thermally couple two mesh domains that exists within
the same mesh region, but whose cells are otherwise not connected one another.

In particular, the test case consists of two parallel 1-D channels, both
with their inlet and outlet. Axially, a heated cellZone exists only on one
channel, while two cellZones acting as a heat exchanger primary and secondary
respectively exist on each channel, allowing the heated fluid in the heated
channel to exchange heat with the cooler fluid flow in the unheated channel.

The `onePhase` case showcases the feature for a single phase flow, while the
`twoPhase` case showcases the feature for two-phase flow, in which the
vapour condenses in the heatExchanger cellZone and heats up the fluid on the
other heatExchanger side.

FunctionObjects (defined at the end of the `controlDict` for both cases) allow
to track quantities of interest at the end of each time-step, such as
the liquid inlet and outlet temperature and mass flow on both channels,
allowing to check for energy conservation.

The specific heatExchanger feature options related to its usage are commented
in the `constant/fluidRegion/phaseProperties` dictionary of each case.
