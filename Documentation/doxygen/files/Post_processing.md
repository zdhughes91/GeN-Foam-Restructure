# Postprocessing {#POSTPROCESSING}

Postprocessing can be performed using paraFoam, the standard post-processing tool used with OpenFOAM.  paraFoam is launched using the command line:

`paraFoam -region (region name)`

where the region name is fluidRegion, neutroRegion or thermalMechanicalRegion (without the parenthesis!).

Please notice that paraFoam is essentially an extension of paraview and it requires having paraview installed. In the openfoam.com distribution, paraview is not distributed with OpenFOAM, but need to be installed separately (see the [ParaView websote](https://www.paraview.org/)). In Ubuntu, it is normally enough to type in the terminal:

`sudo apt-get -y install paraview`

In case of parallel calculations, one should first reconstruct each one of the three meshes using the command

`recontructPar -region (region name)`

where the region name is once again fluidRegion, neutroRegion or thermalMechanicalRegion.

Beside paraFoam, GeN-Foam also creates, in the case folder (or in the processor folder for parallel simulations),  the GeN-Foam.dat file that summurizes few main quantity of interests: time(s), keff(-), power(W), flux0 (m-2s-1), TFuel_Max, TFuel_Avg TFuel_Min, TCladding_Max, TCladding_Avg, TCladding_Min.

useful information is also stored in the log file. The log file can be created by adding the '| tee logFile' command to the launch command, i.e.:

`GeN-Foam | tee logFile`

or

`mpirun -np (number of processors) GeN-Foam -parallel | tee logFile`

Python can effectively be used to extract information from the logFile (several examples available in the tutorials).

In addition, OpenFOAM  allows to use [function objects](https://www.openfoam.com/documentation/guides/latest/doc/guide-function-objects.html) to extract specific information after the simulation.  Also in this case it is essential to indicate which region you want the function object to be applied to, for instance:

`postProcess -func singleGraph -region neutroRegion`

Function objects can also employed at run time via the *controlDict* (see for instance [2D_FFTF](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/blob/master/Tutorials/2D_FFTF/rootCase/system/controlDict)). 

