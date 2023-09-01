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

\*---------------------------------------------------------------------------*/

#include "multiphysicsControl.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multiphysicsControl, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiphysicsControl::multiphysicsControl
(
    const Time& runTime,
    fvMesh& THMesh,
    fvMesh& NMesh,
    fvMesh& TMMesh
)
:
    customPimpleControl
    (
        THMesh,
        "PIMPLE"
    ),
    runTime_(runTime),
    thermalHydraulicMesh_(THMesh),
    neutronicMesh_(NMesh),
    thermoMechanicMesh_(TMMesh),
    topLevelDict_
    (
        IOobject
        (
            "fvSolution",
            runTime.system(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    thermalHydraulicDict_(THMesh.solutionDict()),
    neutronicDict_(NMesh.solutionDict()),
    thermoMechanicDict_(TMMesh.solutionDict()),
    tightlyCoupled_(topLevelDict_.get<bool>("tightlyCoupled")),
    timeStepResidual_(topLevelDict_.get<scalar>("timeStepResidual")),
    maxTimeStepIterations_(topLevelDict_.get<label>("maxTimeStepIterations")),
    liquidFuel_
    (
        runTime.controlDict().lookupOrDefault<bool>("liquidFuel", false)
    )
{
    read();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Copied straight from customPimpleControl to modify infos
bool Foam::multiphysicsControl::loop()
{
    read();
 
    ++corr_;
 
    setFirstIterFlag();
 
    if (corr_ == nCorrPIMPLE_ + 1)
    {
        if (!residualControl_.empty() && (nCorrPIMPLE_ != 1))
        {
            Info<< "Outer iterations not converged within "
                << nCorrPIMPLE_ << " iterations" << endl;
        }
 
        corr_ = 0;
        mesh_.data::remove("finalIteration");
        return false;
    }
 
    bool completed = false;
    if (converged_ || customPimpleControl::criteriaSatisfied())
    {
        if (converged_)
        {
            Info<< algorithmName_ << ": converged in " << corr_ - 1
                << " iterations" << endl;
 
            mesh_.data::remove("finalIteration");
            corr_ = 0;
            converged_ = false;
 
            completed = true;
        }
        else
        {
            //- Neutronics and thermoMechanics are solved on the last iteration
            //  only anyway, so print Outer loop iteration info only if solving
            //  any of them
            if (solveFlow_ or solveEnergy_)
                Info<< "Outer iteration " << corr_ << endl;
            storePrevIterFields();
 
            mesh_.data::add("finalIteration", true);
            converged_ = true;
        }
    }
    else
    {
        if (finalIter())
        {
            mesh_.data::add("finalIteration", true);
        }
 
        if (corr_ <= nCorrPIMPLE_)
        {
            if (solveFlow_ or solveEnergy_)
                Info<< "Outer iteration " << corr_ << endl;
            storePrevIterFields();
            completed = false;
        }
    }
 
    return !completed;
}

bool Foam::multiphysicsControl::read()
{
    customPimpleControl::read();
    nCorrPIMPLE_ = topLevelDict_.get<label>("nOuterCorrectors");
    solveFlow_ = 
        runTime_.controlDict().lookupOrDefault<bool>
        (
            "solveFluidMechanics", true
        );
    solveEnergy_ =
        runTime_.controlDict().lookupOrDefault<bool>
        (
            "solveEnergy", false
        );
    solveNeutronics_ =
        runTime_.controlDict().lookupOrDefault<bool>
        (
            "solveNeutronics", false
        );
    solveThermoMechanics_ =
        runTime_.controlDict().lookupOrDefault<bool>
        (
            "solveThermalMechanics", false
        );

    return true;
}

// ************************************************************************* //
