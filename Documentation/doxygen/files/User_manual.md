# User manual {#USERMAN}

<div class="border-box" style='padding:0.1em; margin-left: 4em;  margin-right: 8em;  border: 1px solid gray; background-color:#f2bfbb; color:#05134a'>
<b>Important note on how to use this manual</b>

This short user manual is meant to provide the user with enough information to understand the logic and various options of GeN-Foam. However, consistent with the Doxygen philosophy, and with the objective of minimizing inconsisitencies between documentation and cource code, detailed usage information of non-trivial sub-solvers, behavioural models, etc., are (or will be) included directly in the header (.H) file of the corresponsing classes. This manual provides links to most of these header files. As an alternative, one can search them by using the search function at the top right of the page. The links provided in the documentation will bring you to the Doxygen page of the corresponding header file. In this page, one can find the mentioned usage information under the section *Classes*. It may happen that the header file was not appropriately formatted for Doxygen at the time of its creation, in which case one may have to directly look at the .H file to find the usage information (work in progress to avoid that).
<br><br> In a similar fashion, to help describe the use of complex dictionaries (i.e., input files), this manual provides links to one or more commented dictionaries that are available in the tutorials.
<br><br> Exact keywords for sub-solvers, models, etc can be found in the corresponding header files. However an easier method to find this names consists in the classical OpenFOAM Banana method: write in the dictionary "Banana" (or any funny word you like), and GeN-Foam will normally give you an error and a list of (typically self-explanatory) valid keywords.
<p>
N.B.: **Users are expected to be already familiar with OpenFOAM and nuclear engineering!**. 
</div>

## Some theoretical background

A fairly general theoretical presentation of GeN-Foam is provided in Ref. \cite FIORINA201524. It is recommended to go through this paper before starting to use GeN-Foam. However, the paper is getting quite old and it is recommended to refer to Ref. \cite FIORINA2016212  for the diffusion solver, Ref. \cite FIORINA2017419 for the SP3 solver, Refs. \cite Fiorina2019DetailedOpenFoam \cite Fiorina2015ApplicationCodes for the thermal-mechanic solver and its use for mesh deformation, Ref. \cite Fiorina2019DetailedOpenFoam for the SN solver, and Refs. \cite Radman2019ADesign \cite RADMAN2021111178 \cite RADMAN2021111422 for single- and two-phase thermal-hydraulics.

Before using this manual, we reccomend to go through the introductury lectures to both OpenFOAM and GeN-Foam that are provided in the folder *Documentation/someUsefulDocumentsAndPResentations*. These lectures are taken from an IAEA e-learning course available at https://elearning.iaea.org/m2/course/view.php?id=1286. The course requires registration and a NUCLEUS account, but it should be available to all IAEA member states. 

## Some practical information

Here below a couple of essential points that make GeN-Foam different than most of the other OpenFOAM-based solvers.

**The multi-region approach**

GeN-Foam employs a  multi-region approach to model different physics using different meshes. This implies that the *0*, *constant* and *system* folders of each case contains multiple folders, one for each physics. In particular, the regions *fluidRegion*, *neutroRegion* and *thermoMechanicalRegion* are employed in GeN-Foam for thermal-hydraulics, neutronics and thermal-mechanics. There is no requirement for the three meshes to occupy the same region of space. Consistent mapping of fields is performed and a reference value is given to a field if no correspondence is found in the mesh where its value is being projected from.

<div class="border-box" style='padding:0.1em; margin-left: 4em;  margin-right: 8em;  border: 1px solid gray; background-color:#f2f3fa; color:#05134a'>
<b>The meshes</b>

A dummy mesh must always be present in all physics (region) directories, even if not solved for. The EMPTY case is already provided with minimal dummy meshes and consistent fields in the “0” folder. Be careful! In case of parallel calculations all your meshes will have to have a number of cells equal or higher than the number of domains you are decomposing your geometry into. In case you need more cells than what available in the EMPTY case, you can run a refineMesh.
</div>


**The multi-zone approach**

In order to assign different properties (for instance, different porous medium properties or different cross sections) to different zones in a mesh, GeN-Foam employs the OpenFOAM concept of cellZone. Each mesh should then be divided in different cellZones. Each cellZone is associated to a name and this name is used in *constant/fluidRegion/phaseProperties*, *constant/neutroRegion/nuclearData*, *constant/themoMechanicalRegion/themoMechanicalProperties* to associate each cellZone with a set of properties. The creation of cellZones is normally allowed by all meshers, though different names are normally used (for instance, *physical entities* in gmsh and *groups* in Salome). In some cases, conversion of the mesh into an OpenFOAM format creates cellSet instead of cellZones. In these case, one can use the topoSet utility to convert cellSets into cellZones.

## The source code

GeN-Foam is an open-source code and makes use in its programming of fairly high-level API, intuitive naming of variables, and frequent comments. As such, we encourage users to consider the code itlsef as an essential part of the documentation. 

The source code is subdivided into 3 main folders:
* main: containing the main  GeN-Foam.C source file (and other files directly employed by it), which is nothing but a fairly complex  coupling loop that calls various functionalities that are found under  *neutronics*, *thermalHydraulics*, and *thermoMechanics*;
* classes: containing the 3 main classes (or sub-libraries of classes)  employed to solve for neutronics, thermal-hydraulics and thermal-mechanics, as well as a class for multi-physics controls;
* include: containing specialized versions of some OpenFOAM base functionalities.

As a general rule, most classes have an *include* folder that is used to store all the .H files that are included via "#include" in the definition of the respective class (normally in the .C file). This is done only to avoid extremely long .C files.


## Content of this manual

The following sections describe the use of the 3 main sub-solvers of GeN-Foam, incuding some theory and useful references:
* [Neutronics](@ref NEUTRONICS)
* [Thermal-hydraulics](@ref TH)
* [Thermal-mechanics](@ref TM)

The following section describes instead the coupling strategy and the general GeN-Foam options.

* [Coupling options and time stepping](@ref COUPLING)




© All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 2021