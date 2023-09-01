# 1D MSR Point-Kinetics

This tutorial displays how to use the point kinetics module of GeN-Foam for
MSRs. It is a simple 1-D case with core, hot leg, pump, heat exchanger and
cold leg. The geometry is one dimensional and salt recirculation is simulated
by making use of a cyclic boundary condition between top and bottom boundaries.

Three simulations are performed:
- a first simulation couples energy and fluid dynamics to obtain a steady state
  at the desired power level.
- a second simulation couples energy, fluid dynamics and point kinetics to
  simulate a loss-of-flow.
- a third simulation is run to allow GeN-Foam to recalculate the reactivity
  loss due to recirculation of the delayed neutron precursors (which can be
  used to verify the results vs analytical results).

Please notice that the power for the steady state has been imposed via the
`powerDensity` field specified in the `fluidProperties` sub-dictionary in the
`constant/fluidRegion/phaseProperties` dictionary and models a constant power
density in the core and null in the other regions. This power is then
automatically updated by the point kinetics solver during the transient.

Please notice also that a correct evaluation of the reactivity worth of
delayed neutron precursors in MSRs would require knowledge of the adjoint
flux. In GeN-Foam, the adjoint flux is approximated by the oneGroupFlux.
When fluxes are not calculated via a diffusion calculation, one has to
manually provide the oneGroupFlux in `0/neutroRegion`. In this tutorial,
the oneGroupFlux has been set to 1 in the core and zero elsewhere. This
is done via the `initialOneGroupFluxByZone` keyword in `nuclearData`.

A few python files are provided to plot essential results (for instance:
`python plotPKPower_Temp.py ./transient/log.GeN-Foam`). In addition, a
`.m` file (that can be run using Octave) is provided that calculates expected
results at the end of the transient.

The tutorial has been prepared by Arnaldo Mattioli (Politecnico di Milano) and
revised by Carlo Fiorina (EPFL), Stefan Radman (EPFL).
