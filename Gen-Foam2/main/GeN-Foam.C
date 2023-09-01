/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
|                                                                             |
|    Built on OpenFOAM v2212                                                  |
|    Copyright 2011-2016 OpenFOAM Foundation, 2017-2022 OpenCFD Ltd.         |
-------------------------------------------------------------------------------
License
    This file is part of GeN-Foam.

    GeN-Foam is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    GeN-Foam is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    This offering is not approved or endorsed by the OpenFOAM Foundation nor
    OpenCFD Limited, producer and distributor of the OpenFOAM(R)software via
    www.openfoam.com, and owner of the OPENFOAM(R) and OpenCFD(R) trademarks.

    This particular snippet of code is developed according to the developer's
    knowledge and experience in OpenFOAM. The users should be aware that
    there is a chance of bugs in the code, though we've thoroughly test it.
    The source code may not be in the OpenFOAM coding style, and it might not
    be making use of inheritance of classes to full extent.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    GeN-Foam

Description
    Multi-physics solver for nuclear reactor analysis. It couples a multi-scale
    fine/coarse mesh 3-phase (liquid, vapour, porous substructure) sub-solver 
    for thermal-hydraulics, various sub-solvers for neutronics,
    a displacement-based sub-solver for thermal-mechanics. The 
    thermal-hydraulic sub-solver consists of the custom developed FFSEulerFoam 
    solver  (https://gitlab.com/virmodoetiae/FFSEulerFoam). It is capable of
    modelling single and two-phase flows, while modelled fuel types consist
    of either liquid fuel (e.g. MSRs) or fuel pin lattices. For the latter,
    the energy dynamics is represented via a 1.5-D finite difference model.

    Reference publications:
    
    NOTE: these publications do not cover recent multi-phase developments

    Carlo Fiorina, Ivor Clifford, Manuele Aufiero, Konstantin Mikityuk, 2015
    "GeN-Foam: a novel OpenFOAM® based multi-physics solver for 2D/3D transient
    analysis of nuclear reactors", Nuclear Engineering and Design 294, pp. 
    24-37

    Carlo Fiorina, Konstantin Mikityuk, " Application of the new GeN-Foam 
    multi-physics solver to the European Sodium Fast Reactor and verification 
    against available codes", Proceedings of ICAPP 2015, May 03-06, 2015 - 
    Nice (France), Paper 15226

    Authors of this file (and associated .C or included .H files): 
    Carlo Fiorina <carlo.fiorina@outlook.com; carlo.fiorina@epfl.ch;>
    Stefan Radman <stefanradman92@gmail.com; stefan.radman@epfl.ch;>
    EPFL (Switzerland)

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "SquareMatrix.H"
#include "fvMatrixExt.H"
#include "meshToMesh.H"
#include "regionProperties.H"
#include "mergeOrSplitBaffles.H"
#include "volPointInterpolation.H"
#include "fixedGradientFvPatchFields.H"
#include "UPstream.H"

#include "multiphysicsControl.H"
#include "thermalHydraulicsModel.H"
#include "neutronics.H"
#include "thermoMechanics.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #define NO_CONTROL
    #define CREATE_MESH createMeshesPostProcess.H
    
    #include "postProcess.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMeshes.H"
    #include "createFields.H"
    #include "createMeshInterpolators.H"
    #include "createCouplingFields.H"
    #include "createOutput.H"

    Info<< "\nStarting time loop - \n" << endl;

    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s" 
        << nl << endl;

    #include "setDeltaT.H"

    while (runTime.run())
    {
        runTime++;

        #include "setDeltaT.H"
        
        Info << "Time = " << runTime.timeName() << nl << endl;

        while (multiphysics.loop())
        {
            #include "solve.H"
        }

        runTime.write();

        #include "writeOutput.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
