# Running GeN-Foam {#RUNNING}

GeN-Foam is launched like any OpenFOAM solver, by executing the following commands in a terminal (after sourcing the OpenFOAM environment):

`GeN-Foam`

or

`mpirun -np (number of processors) GeN-Foam -parallel`

In case of parallel calculations, one should decompose each one of the three meshes using the command

`decomposePar -region (region name)`

where the region name is fluidRegion, neutroRegion or thermalMechanicalRegion. `decomposePar` operates according to the standard OpenFOAM *decomposeParDict* to be placed in *system/(region name)*

N.B.:  all meshes must be decomposed in the same number of domains, and a consistent decomposeParDict must also be present in *system/*

In the tutorials, *Allrun* (and sometimes *Allrun_parallel*) bash scripts are provided that can be used to run the case (one should just type `Allrun` or `Allrun_parallel` in a terminal  ), including mesh decomposition for parallel cases, as well as multiple GeN-Foam runs that are used for instance to: 1) achieve a steady-state; 2) run a transient starting from that steady-state.
