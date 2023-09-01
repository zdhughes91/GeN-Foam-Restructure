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
    Foam::thermalHydraulicsModels::onePhase

Description
    This class encapsulates a single fluid + 1 stationary structure 
    Eulerian solver with a user-defineable regime-map based phyiscs modelling.

    Author: 
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

SourceFiles
    onePhase.C

\*---------------------------------------------------------------------------*/

#ifndef onePhase_H
#define onePhase_H

#include "thermalHydraulicsModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Forward declarations
//class regime;
//class regimeMapModel;

namespace thermalHydraulicsModels
{

/*---------------------------------------------------------------------------*\
                      Class onePhase Declaration
\*---------------------------------------------------------------------------*/

class onePhase
:
    public thermalHydraulicsModel
{
private:

    //- Phases

        //- Structure
        structure structure_;

        //- Fluid
        fluid fluid_;

    //- Phase pairs
        
        //- Fluid2 - structure
        FSPair FSPair_;
    /*
    //- Phase interaction related  
        
        //- table of drag coefficients between the fluid and each
        //  structure region (currently only defined in terms of cellZone)
        volTensorFieldPtrTable Kds_;

        //- table of heat transfer coefficients between the fluid and each
        //  structure region
        volScalarFieldPtrTable htcs_;

        //- Regime map
        autoPtr<regimeMapModel> regimeMap_;
    */

public:

    TypeName("onePhase");

    //- Constructors

        onePhase
        (
            Time& time,
            fvMesh& mesh, 
            customPimpleControl& pimple, 
            fv::options& fvOptions
        );

        //- Destructor
        virtual ~onePhase(){}

    //- Member Functions

        //- Solve thermal-hydraulics
        virtual void correct
        (
            scalar& residual,
            bool solveFluidDynamics, 
            bool solveEnergy
        );

        //- Correct fluid and structure temperatures without
        //  updating the flow regime map
        virtual void correctEnergy(scalar& residual);

        //- Correct pressure, velocities without updating 
        //  the flow regime map
        virtual void correctFluidMechanics(scalar& residual);

        //- Correct models and coupling fields
        virtual void correctModels
        (
            bool solveFluidDynamics, 
            bool solveEnergy
        );

        //- Update and print Courant number
        virtual void correctCourant();

        //- Update continuity errors of the fluid
        void correctContErr();

        //- Print continuity error information
        void printContErr();

        //- Computes the cumulated continuity error
        void calcCumulContErr();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace thermalHydraulicsModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
