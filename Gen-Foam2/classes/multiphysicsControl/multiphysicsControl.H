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

Class
    Foam::multiphysicsControl

Description
    Class to handle fvSolution dictionaries for the three regions.

    Author of this file (and associated .C or included .H files):  
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com), EPFL (Swizerland)

SourceFiles
    multiphysicsControl.C

\*---------------------------------------------------------------------------*/

#ifndef multiphysicsControl_H
#define multiphysicsControl_H

#include "fvMesh.H"
#include "customPimpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class multiphysicsControl Declaration
\*---------------------------------------------------------------------------*/

class multiphysicsControl
: 
    public customPimpleControl
{
protected :

    //- Ref to runTime
    const Time& runTime_;

    //- Ref to thermalHydraulics mesh
    fvMesh& thermalHydraulicMesh_;

    //- Ref to neutronics mesh
    fvMesh& neutronicMesh_;

    //- Ref to thermoMechanics mesh
    fvMesh& thermoMechanicMesh_;

    IOdictionary        topLevelDict_;
    const dictionary&   thermalHydraulicDict_;
    const dictionary&   neutronicDict_;
    const dictionary&   thermoMechanicDict_;

    //- Solution flags

    //- solveFluidDynamics exists in customPimpleControl as solveFlow_;
    bool solveEnergy_;
    bool solveNeutronics_;
    bool solveThermoMechanics_;

    //- Global loop parameter (i.e. outer iterations)

    bool    tightlyCoupled_;
    scalar  timeStepResidual_;
    label   maxTimeStepIterations_;

    //- liquidFuel flag (true e.g. for MSRs)
    bool liquidFuel_;

public:

    TypeName("multiphysicsControl");

    //- Constructors

        multiphysicsControl
        (
            const Time& runTime,
            fvMesh& THMesh,
            fvMesh& NMesh,
            fvMesh& TMMesh
        );

        //- Destructor
        virtual ~multiphysicsControl() {}

    //- Member Functions

        //- Override PIMPLE so that the total number of outer iterations to
        //  perform is read from the topLevelDict and not from the dictionary
        //  in fluidRegion
        virtual bool read();

        //- Ovveride PIMPLE loop function to modify some Infos (really, 
        //  Stefan?)
        virtual bool loop();

        //- Access

            const dictionary& topLevelDict() const
            {
                return topLevelDict_;
            }

            const dictionary& thermalHydraulicsDict() const
            {
                return thermalHydraulicDict_;
            }

            const dictionary& neutronicsDict() const
            {
                return neutronicDict_;
            }

            const dictionary& thermoMechanicsDict() const
            {
                return thermoMechanicDict_;
            }

            bool solveFluidMechanics() const
            {
                return solveFlow_;
            }

            bool solveEnergy() const
            {
                return solveEnergy_;
            }

            bool solveNeutronics() const
            {
                return solveNeutronics_;
            }

            bool solveThermoMechanics() const
            {
                return solveThermoMechanics_;
            }

            const label& nOuterCorrectors() const
            {
                return nCorrPIMPLE_;
            }

            bool tightlyCoupled() const
            {
                return tightlyCoupled_;
            }

            const scalar& timeStepResidual() const
            {
                return timeStepResidual_;
            }

            const label& maxTimeStepIterations() const
            {
                return maxTimeStepIterations_;
            }

            bool liquidFuel() const
            {
                return liquidFuel_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
