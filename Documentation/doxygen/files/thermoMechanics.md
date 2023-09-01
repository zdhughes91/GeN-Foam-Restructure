
# Thermal-mechanics {#TM}

## Introduction

The thermal-mechanics solver of GeN-Foam is a simple linear elasticity solver that can be used to evaluate thermal deformations in a core. Temperatures are projected from the thermal-hydraulic solver (temperatures of fuels and structures) and the deformation field is used to deform the mesh for neutronics and thermal-hydraulics. In particular, the radial deformation of structures and the axial deformation of fuel are employed. to deform the neutronics mesh. 


## Various properties

<div class="border-box" style='padding:0.1em; margin-left: 4em;  margin-right: 8em;  border: 1px solid gray; background-color:#f2f3fa; color:#05134a'>
<b>The *thermoMechanicalProperties* dictionary</b>

The *thermoMechanicalProperties* dictionary can be found under *constant/thermoMechanicalRegion* and allow to define the thermo-mechanical properties of structures, subdivided according to the cellZones of the thermoMechanicalRegion mesh. 
<br><br>One can find a detailed, commented example in the tutorial 
[3D_SmallESFR](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/tree/master/Tutorials/3D_SmallESFR/rootCase/constant/thermoMechanicalRegion/thermoMechanicalProperties).
</div>
<br>

## Initial and boundary conditions

Besides the standard ones available in OpenFOAM, GeN-Foam includes a *tractionDisplacement* boundary condition that allows to set a pressure or a traction on a boundary. Work is ongoing to adapt to GeN-Foam the contact boundary condition of the OFFBEAT fuel behavior solver \cite SCOLARO2020110416. 


## Discretization and solution

Details for discretization and solution of equations are handled in a standard OpenFOAM way, i.e., through the *fvSolution* and *fvSchemes* dictionaries in *constant/thermoMechanicalRegion*. 

© All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 2021
