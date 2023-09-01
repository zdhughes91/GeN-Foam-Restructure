
# Neutronics {#NEUTRONICS}

## Introduction

The neutronics class (see *neutronics.H*) is a high-level class that contains essential data and variables that are common to various neutronics models. In particular, the neutronics class handles the variables that are included in the *constant/neutroRegion/reactorState* dictionary.

<div class="border-box" style='padding:0.1em; margin-left: 4em;  margin-right: 8em;  border: 1px solid gray; background-color:#f2f3fa; color:#05134a'>
<b>The *reactorState* dictionary</b>

The *reactorState* dictionary is found under the *timeFolder/uniform/* sub-folder. It included essentially 4 keywords:
<UL>
<LI> *keff*  is used in the spatial kinetics solvers as initial guess for keff when doing an eigenvalue calculation. It is then updated automatically at each time step (i.e., at each power iteration) with the calculated value of keff. When performing a transient calculation with the spatial kinetics solvers, *keff* is instead used to divide the neutron source term and is not updated during the simulation. Typically, to run spatial kinetics transient simulations, one first runs an eigenvalue calculation. The resuling *keff* will be the one that makes the reactor critical in a subsequent transient simulation. *keff* is diregarded by the point kinetics sub-solver.
<LI> *pTarget* is used in the spatial kinetics solvers as target power when doing an eigenvalue calculation. It is also used by the point kinetics sub-solver, but only to correctly plot results. As power, GeN-Foam uses what it finds under powerDensity of the neutroRegion, or under the powerDensity of the fluidRegion if it does not find a powerDensity in the neutroRegion. To correctly plot point kinetics result, pTarget must be consistent with the mentioned power densities.
<LI> *externalReactivity* is read both by the spatial kinetics and point kinetics solvers and is used to instantly add (or remove) a certain reactivity at the beginning of a transient.
<LI> *precursorPowers* can be read by the point kinetics sub-solver in case the user wishes to set initial concentrations of precursors. If not found, precursor concentration are initialized so to be in equilibrium with the starting conditions (i.e. a steady state is assumed). 
</UL>

All GeN-Foam neutronics models can be used for liquid-fuel reactors. One can activate this option using the  *liquidFuel* keyword in */system/controlDict*. Of course, in such case one should pay attention to setting proper boundary conditions for the precursors.

A commented  reactorState can be found in
[3D_SmallESFR](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/tree/master/Tutorials/3D_SmallESFR/rootCase/constant/neutroRegion/reactorState). 

NB: Please note that in parallel calculations, the updated *reactorState* can be found in *processor0/constant/neutroRegion/*.
</div>

## Models

Neutronics calculations are performed by classes derived from *neutronics* that contain specific sub-solvers:
* *pointKineticNeutronics* for point kinetics calculations (see *pointKineticNeutronics.H*)
* *diffusionNeutronics* for diffusion calculations (see *diffusionNeutronics.H*)
* *adjointDiffusionNeutronics* for adjoint diffusion calculations (see *adjointDiffusionNeutronics.H*)
* *SP3Neutronics* for diffusion calculations (see *SP3Neutronics.H*)
* *SNNeutronics* for discrete ordinates calculations (see *SNNeutronics.H*)
For the user, the derived classes translate into runtime selectable models. The specific sub-solver to be used in a simulation can be selected at runtime in the *constant/neutroRegion/neutronicsProperties* dictionary. 


<div class="border-box" style='padding:0.1em; margin-left: 4em;  margin-right: 8em;  border: 1px solid gray; background-color:#f2f3fa; color:#05134a'>
<b>The *neutronicsProperties* dictionary</b>

The *neutronicsProperties* dictionary is found under *constant/neutroRegion/* and it can be used to set the type of neutronic simulation by using the following keywords:
<UL>
<LI>  *model*  is used to define what type of simulation needs to be performed. It can be *pointKinetics*, *diffusionNeutronics*, *SP3Neutronics*, **SNNeutronics, *adjointDiffusion*. *adjointDiffusion* has been developed only as an eigenvalue solver. The others can be used for transient calculations. However, the SN transient solver has not been tested. In addition, it is currently not accelrated, thus extremely slow (it can require
hundreds of iterations per time step).
<LI>  *eigenvalueNeutronics* should be set to *true* for eigenvalue calculations, false for transients.
</UL>
One can find detailed, commented examples in most tutorials. See for instance 
[3D_SmallESFR](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/tree/master/Tutorials/3D_SmallESFR/rootCase/constant/neutroRegion/neutronicsProperties) (single phase).
</div>

## Various properties

In GeN-Foam, cross-sections and several other neutronics properties are handled by the *XS.H* class, or by its low-memory version *XSLowMem.H*.

<div class="border-box" style='padding:0.1em; margin-left: 4em;  margin-right: 8em;  border: 1px solid gray; background-color:#f2f3fa; color:#05134a'>
<b>The *nuclearData* dictionaries</b>

The *nuclearData* dictionary can be found under *constant/neutroRegion/*. It contains all basic nuclear properties for the reference reactor state. The other *nuclearData...* files in *constant/neutronics/* should include the cross-sections for perturbed reactor states. In addition, these files include information about the perturbed and reference (*nuclearData*) reactor state. For instance, *nuclearDataFuelTemp* must include *TfuelRef* and *TfuelPerturbed*, which represent the temperatures at which the reference (*nuclearData*) and perturbed  (*nuclearDataFuelTemp*) cross sections have been calculated, respectively. Linear interpolation is performed by GeN-Foam between reference and perturbed reactor states, except for fuel temperature, for which a logarithmic or square root interpolation is provided (depending on the spectrum, which in turns is defined by the keyword *fastNeutrons*).  If no data are provided, the reference cross sections are used. Nuclear data can be generated using any nuclear code.  The [serpentToFoam](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/tree/master/Tools/serpentToFoam/serpent2.1.23) routines provided with GeN-Foam (in the *Tools* folder) is an Octave script that automatically converts Serpent output files into the nuclear data files employed by GeN-Foam. The entry *discFactor* is used only if discontinuity factors have to be used. The term *integralFlux*, is used only if the automatic adjustment of discontinuity factors is performed \cite FIORINA2016212. Nonetheless, these entries should always be present. 
<br><br>One can find detailed, commented examples of nuclearData in the tutorials 
[3D_SmallESFR](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/tree/master/Tutorials/3D_SmallESFR/rootCase/constant/neutroRegion/nuclearData) (for diffusion or SP3),
[Godiva_SN](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/tree/master/Tutorials/Godiva_SN/constant/neutroRegion/nuclearData) (for discrete ordinates) and 
[2D_onePhaseAndPointKineticsCoupling](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/tree/master/Tutorials/2D_onePhaseAndPointKineticsCoupling/rootCase/constant/neutroRegion/nuclearData) (for point kinetics).

N.B.: cross sections must be expressed according to the International System of Units (so m, not cm).

N.B.2: defaultPrec has 1/m3 units except for the adjoint solver that needs 1/m2/s.

N.B.3: the *nuclearData...* files must always be present, even when not parametrizing cross-sections. If no parametrization is needed, the “zone” card must be left “blank” as:
`zones();`
</div>

An additional dictionary is needed to provide the quadrature set when performing discrete ordinate calculations.

<div class="border-box" style='padding:0.1em; margin-left: 4em;  margin-right: 8em;  border: 1px solid gray; background-color:#f2f3fa; color:#05134a'>
<b>The *quadratureSet* dictionary</b>

The *quadratureSet* dictionary is found under *constant/neutroRegion/*. It ccontains the quadrature set for discrete ordinate calculations. 
<br><br>One can find examples of three different quadrature set in the tutorial 
[Godiva_SN](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/tree/master/Tutorials/Godiva_SN/constant/neutroRegion/). 
S4 and S8 chebichev Legendre quadrature sets can be found in [Godiva_SN](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/tree/develop/Tools/chebichevLegendreQuadratureSets/)
</div>

Finally, the *CRMove* dictionary can be used to move control rods. 

<div class="border-box" style='padding:0.1em; margin-left: 4em;  margin-right: 8em;  border: 1px solid gray; background-color:#f2f3fa; color:#05134a'>
<b>The *CRMove* dictionary</b>

The CRMove* dictionary can be found under *constant/neutroRegion/*. It contains input data for control rods movement. Control rods can be moved from the initial position to a new one by selecting initial and final time of the insertion/extraction and the speed of insertion/extraction (positive speed for insertion).
<br><br>One can find a commented example in the tutorial 
[3D_SmallESFR](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/tree/master/Tutorials/3D_SmallESFR/rootCase/constant/neutroRegion/CRmove), though this option is not actually used in the tutorial.
</div>


## Initial and boundary conditions

As in all standard OpenFOAM solvers, initial values (IC) and boundary conditions (BC) should be provided in the “0” folder, or in the folder corresponding to the *startTime* of the simulation, if different than 0. In the case of neutronics, the user can either specify different IC and BC for each one of the energy groups (with fluxes that must be named *fluxStar0*, *fluxStar1*, etc…), or provide the same IC and BC to all fluxes by using the *defaultFlux* field. In case of SP3 calculations, the IC and BC for the second moment can be imposed either for each energy (using fields named *fluxStar20*, *fluxStar21*, etc…), or to all energies by using the *defaultFlux2* field. When both *defaultFlux* and *fluxStar...* are present, the solver gives priority to *fluxStar...*. In case of SN calculations, it is suggested not to modify the boundary conditions and to use the *defaultFlux* file (an example is providedin the Godiva_SN tutorial). When employing the adjoint solver, you will have to add the fields *adjointDefaultPrec* and *adjointDefaultFlux* in your initial time.

In addition to the standard OpenFOAM BC, an albedo boundary condition (see *albedoSP3FvPatchField.H*) is available in GeN-Foam for diffusion and SP3 calculations and can be used according to the following syntax:

```
	type            albedoSP3;
	gamma		0.5; // defined as (1-alpha)/(1+alpha)/2, alpha being the albedo coefficient
	diffCoeffName	Dalbedo;  //not to be changed
	fluxStarAlbedo  fluxStarAlbedo; //not to be changed
	forSecondMoment false;  //true in case it is a condition for a second moment flux (for SP3 calculations)
	value           uniform 1;
```
Please note that the boundary condition needs to be set both for first and second moments in SP3.

IC and BC for precursors do not have to be specified for standard reactors. On the other hand, they should be specified in case of liquid fuel reactors (e.g., Molten Salt Reactors). This is possible by creating a *defaultPrec* field, in case the same conditions apply to all precursor groups, or by creating the fields named *prec0*, *prec1*, etc., in case different conditions must be provided for different precursor groups. 

N.B.: boundary conditions must be applied to *fluxStar...* and not to *flux...* since GeN-Foam solves for these variables. *fluxStar...* represent continuous fluxes, while *flux...* represent the real fluxes. They differ only in case discontinuity factors are employed \cite FIORINA2016212. 


<div class="border-box" style='padding:0.1em; margin-left: 4em;  margin-right: 8em;  border: 1px solid gray; background-color:#f2f3fa; color:#05134a'>
<b>Setting the weighting in point-kinetics calculations</b>

A correct evaluation of the reactivity worth of delayed neutron precursors in MSRs, as well as of the impact of temperatures on reactivities, normally  requires the knowledge of the adjoint flux. In GeN-Foam, the field *oneGroupFlux* is used by the point kinetic solver for weighting temperatures, densities and precursors. When fluxes are not calculated via a spatial neutronics calculation, one has to manually provide the *oneGroupFlux* in *0/neutroRegion*. As an alternative, one can use the *initialOneGroupFluxByZone* keyword in *nuclearData* (see [1D_MSR_pointKinetics](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/blob/master/Tutorials/1D_MSR_pointKinetics/rootCase/constant/neutroRegion/nuclearData)). Please notice that:
<UL>
<LI>  If calcuclated fluxes are available in *neutroRegion*, these will be user to recacluate and overwrite *oneGroupFlux*.
<LI>  If no fluxes are available, the neutronic-subsolver will use the provided *oneGroupFlux*
<LI>  If the *initialOneGroupFluxByZone* keyword in used in *nuclearData*, this will be used to overwrite *oneGroupFlux*
</UL>
</div>



## Discretization and solution

Details for discretization and solution of equations are handled in a standard OpenFOAM way, i.e., through the *fvSolution* and *fvSchemes* dictionaries in *constant/neutroRegion*. 



© All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 2021
