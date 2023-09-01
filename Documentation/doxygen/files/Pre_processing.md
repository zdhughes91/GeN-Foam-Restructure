# Preprocessing {#PREPROCESSING}

Before running GeN-Foam, one has to provide meshes, physical properties, discretization methods (if one does not want to use the default ones) and simulation details.  This section provides a quick overview of the input deck of GeN-Foam. For more details, please refer to the [User manual](@ref USERMAN).


## Meshing

GeN-Foam uses three different meshes for neutronics, thermal-hydraulics and thermal-mechanics. There is no requirement for the three meshes to occupy the same region of space. Consistent mapping of fields is performed and a reference value is given to a field if no correspondence is found in the mesh where its value is being projected from. It follows that the geometry for neutronics can cover only a small part of the overall reactor geometry.  Meshes can be created with every OpenFOAM-compatible tool. Meshes can (should) be divided into zones (cellZones) to allow the use of different physical properties (e.g., cross sections) in different reactor regions. Sometimes, when converting a mesh to the OpenFOAM (polymesh) format, cellSets (and not cellZones) are created. The topoSetDict can be used to convert cellSets to cellZones.

The 3D_SmallESFR tutorial includes an example of mesh generation with Gmsh \cite https://doi.org/10.1002/nme.2579. The folder includes the subfolder *3dESFRMesh* containing three subfolders for the generation of the meshes for neutronics, thermal-hydraulics and thermal-mechanics. To create the 3D ESFR core geometry, one should execute the following command in a terminal: 

`gmsh esfMain.geo`

This will create (after a relatively long time) a geometry and open the Gmsh graphical interface. The geometry must then be meshed and the resulting mesh saved and copied in the root of the ESF3D_small folder. 

By typing in a terminal:

`gmshToFoam (name of mesh file)`

a *polyMesh* folder will be created (or updated) in the folder *constant*. One should then copy this folder in *constant/neutroRegion*, *constant/fluidRegion* or *constant/thermoMechanicalRegion*, and repeat the operation for all the three meshes.

Please notice that the 3D_SmallESFR tutorial already contains the correct *polymesh* folders so that one can avoid the mesh generation step.

N.B. **A dummy mesh must always be present in all physics (region) directories**, even if not solved for. The EMPTY case is already provided with minimal dummy meshes and consistent fields in the “0” folder. Be careful! In case of parallel calculations all your meshes will have to have a number of cells equal or higher than the number of domains you are decomposing your geometry into. In case you need more cells than what available in the EMPTY case, you can run a refineMesh

## Physical properties

All the data for the GeN-Foam simulations can be filled in the following input files (dictionaries):
* *constant/thermoMechanicalRegion/thermoMechanicalProperties* - thermo-mechanical properties of structures, subdivided according to the cellZones of the thermoMechanicalRegion mesh. 
One can find a detailed, commented example in the tutorial 
[3D_SmallESFR](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/tree/master/Tutorials/3D_SmallESFR/rootCase/constant/thermoMechanicalRegion/thermoMechanicalProperties).

* *constant/fluidRegion/g* - gravitational acceleration.

* *constant/fluidRegion/turbulenceProperties* - standard OpenFOAM dictionary to define the turbulence model to be used.
One can find a detailed, commented example in the tutorial 
[3D_SmallESFR](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/tree/master/Tutorials/3D_SmallESFR/rootCase/constant/fluidRegion/turbulenceProperties).


* *constant/fluidRegion/thermophysicalProperties* (for single-phase simulations) - standard OpenFOAM dictionary to define the thermo-physical properties of the coolant.
One can find a detailed, commented example in tutorial 
[3D_SmallESFR](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/tree/master/Tutorials/3D_SmallESFR/rootCase/constant/fluidRegion/thermophysicalProperties) (single phase) 

* *constant/fluidRegion/thermophysicalProperties.(name of fluid)* (for two-phase simulations) - standard OpenFOAM dictionaries to define the thermo-physical properties of various phases. The name of fluid is defined in *constant/fluidRegion/phaseProperties*.
One can find a detailed, commented example in the tutorial 
[1D_boiling (liquid)](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/tree/master/Tutorials/1D_boiling/constant/fluidRegion/thermophysicalProperties.liquid),
[(vapour)](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/tree/master/Tutorials/1D_boiling/constant/fluidRegion/thermophysicalProperties.vapour).


* *constant/fluidRegion/phaseProperties* - large dictionary that can be used to: determined whether the simulation is single-phase or two-phase; set various properties of the phases (beside the thermo-physical properties defined in *constant/fluidRegion/thermophysicalProperties*); set the properties of the sub-scale structures (fuel pins, heat exchangers, etc) in the porous zones, including the possibility to assign a *powerModel* for power production (e.g., nuclear fuel, or constant power) and the *passiveProperties* of another sub-structure that interacts thermally with the fluid (for instance the wrappers in sodium fast reactors).  The name of the porous zones must coincide with that of the cellZones of the fluidRegion mesh. Anisotropic pressure drops can be set by using the keywords *transverseDragModel* (Blasius, GunterShaw, same) and *principalAxis*(localX, localY, localZ) in the sub-dictionary *dragModels.(nameOfPhase).structure.(nameOfCellZones)*. *principalAxis* sets the axis on which the nominal dragModel is used. *transverseDragModel* sets the model to be used on the two directions that are perpendicular to *principalAxis*. If *same* is chosen as *transverseDragModel*, the code will use the nominal model in all directions, but with the possibility of an anisotropic hydraulic diameter. The anisotropy of the hydraulic diameter can be set using the keyword *localDhAnisotrpy* and assign to it a vector of 3 scaling factors (one for each local directions). 
One can find detailed, commented examples in the tutorials 
[3D_SmallESFR](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/tree/master/Tutorials/3D_SmallESFR/rootCase/constant/fluidRegion/phaseProperties) (single phase) and
[1D_boiling](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/tree/master/Tutorials/1D_boiling/constant/fluidRegion/phaseProperties) (two phases).

* *constant/neutroRegion/neutronicsProperties* - dictionary to control how neutronics is solved (point kinetics, diffusion, SP3 or SN), and if it's an eigenvalue calculation or a transient. 
One can find detailed, commented examples in most tutorials. See for instance 
[3D_SmallESFR](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/tree/master/Tutorials/3D_SmallESFR/rootCase/constant/neutroRegion/neutronicsProperties) (single phase).

* *constant/neutroRegion/reactorState* – contains the target power (pTarget) for eigenvalue calculations, the keff that results from the eigenvalue calculations and the external reactivity (i.e., the extra reactivity one can add for instance to simulate a reactivity step). N.B.: keff has no effect on pointKinetics. You can find detailed, commented examples in most tutorials. N.B.2: In point kinetics, pTarget is the initial value used by the point kinetics solver to plot results, but the solver actually scale the powerDensity and flux fields provided by the user. It is up to the user to make sure that pTarget is consistent with the powerDensity and flux fields. A commented  reactorState can be found in
[3D_SmallESFR](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/tree/master/Tutorials/3D_SmallESFR/rootCase/constant/neutroRegion/reactorState) (single phase). Please note that eigenvalue calculations will update the keff value in this dictionary. In parallel calculations, the updated value can be found in *processor0/constant/neutroRegion/reactorState*.

* *constant/neutroRegion/nuclearData* - contains all basic nuclear properties for the reference reactor state. The other *nuclearData...* files in *constant/neutronics/* should include the cross-sections for perturbed reactor states. In addition, these files include information about the perturbed and reference (*nuclearData*) reactor state. For instance, *nuclearDataFuelTemp* must include *TfuelRef* and *TfuelPerturbed*, which represent the temperatures at which the reference (*nuclearData*) and perturbed  (*nuclearDataFuelTemp*) cross sections have been calculated, respectively. Linear interpolation is performed by GeN-Foam between reference and perturbed reactor states, except for fuel temperature, for which a logarithmic or square root interpolation is provided (depending on the spectrum, which in turns is defined by the keyword *fastNeutrons*).  If no data are provided, the reference cross sections are used. Nuclear data can be generated using any nuclear code.  The [serpentToFoam](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/tree/master/Tools/serpentToFoam/serpent2.1.23) routines provided with GeN-Foam (in the *Tools* folder) is an Octave script that automatically converts Serpent output files into the nuclear data files employed by GeN-Foam. The entry *discFactor* is used only if discontinuity factors have to be used. The term *integralFlux*, is used only if the automatic adjustment of discontinuity factors is performed \cite FIORINA2016212. Nonetheless, these entries should always be present. 
One can find detailed, commented examples of nuclearData in the tutorials 
[3D_SmallESFR](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/tree/master/Tutorials/3D_SmallESFR/rootCase/constant/neutroRegion/nuclearData) (for diffusion or SP3),
[Godiva_SN](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/tree/master/Tutorials/Godiva_SN/constant/neutroRegion/nuclearData) (for discrete ordinates) and 
[2D_onePhaseAndPointKineticsCoupling](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/tree/master/Tutorials/2D_onePhaseAndPointKineticsCoupling/rootCase/constant/neutroRegion/nuclearData) (for point kinetics).
One can find examples of the *nuclearData...* files in the tutorial 
[3D_SmallESFR](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/tree/master/Tutorials/3D_SmallESFR/rootCase/constant/neutroRegion) 

* *constant/neutroRegion/quadratureSet* - contains the quadrature set for discrete ordinate calculations. 
One can find examples of three different quadrature set in the tutorial 
[Godiva_SN](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/tree/master/Tutorials/Godiva_SN/constant/neutroRegion/).

* *constant/neutroRegion/CRMove* - contains input data for control rods movement. Control rods can be moved from the initial position to a new one by selecting initial and final time of the insertion/extraction and the speed of insertion/extraction (positive speed for insertion).
One can find a commented example in the tutorial 
[3D_SmallESFR](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/tree/master/Tutorials/3D_SmallESFR/rootCase/constant/neutroRegion/CRmove), though this option is not actually used in the tutorial.



## Initial values and boundary conditions

As in all standard OpenFOAM solvers, initial values (IC) and boundary conditions (BC) should be provided in the “0” folder, or in the time folder corresponding to the *startTime* of the simulation, if different than 0. 


## Discretization and solution

Details for discretization and solution of single-physics equations are handled in a standard OpenFOAM way, i.e., through the *fvSolution* and *fvSchemes* dictionaries in *constant/neutroRegion*, *constant/fluidRegion* or *constant/ thermalMechanicalRegion*. Coupling and simulation details are determined through the *fvSolution* and *controlDict* in the *system* folder. The *controlDict* is significantly extended compared to a standard OpenFOAM controlDict in order to control what solvers are activated and flags that affect the behaviour of GeN-Foam as a whole.



