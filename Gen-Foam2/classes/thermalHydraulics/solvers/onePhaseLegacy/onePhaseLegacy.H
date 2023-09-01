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
    Foam::thermalHydraulicsModels::onePhaseLegacy

Description
    This class is based on onePhase.H and implements the legacy porous medium 
    Navier Stokes solver with upgrades from the current formulations in onePhase.H.

    Author: 
    
    Mixer: Peter German (peter.german@tamu.edu / grmnptr@gmail.com)

    Authors of the original solvers: 
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)
    Carlo Fiorina (carlo.fiorina@epfl.ch / carlo.fiorina@outlook.com)

SourceFiles
    onePhaseLegacy.C

\*---------------------------------------------------------------------------*/

#ifndef onePhaseLegacy_H
#define onePhaseLegacy_H

#include "thermalHydraulicsModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Forward declarations
class regime;
class regimeMapModel;
class saturationModel;

namespace thermalHydraulicsModels
{

/*---------------------------------------------------------------------------*\
                      Class onePhaseLegacy Declaration
\*---------------------------------------------------------------------------*/

class onePhaseLegacy
:
    public thermalHydraulicsModel
{
public:

    enum class phaseChangeModel
    {
        none,
        heatConductionLimited, 
        individualSideDriven
    };

    static const Enum<phaseChangeModel> phaseChangeModelNames;

private:

    //- Phases

        //- Structure
        structure structure_;

        //- Fluid
        fluid fluid_;

    //- Phase pairs
        
        //- Fluid - structure
        FSPair FSPair_;

    //- Container for the cumulative continuity error
        scalar cumulativeContErr_;

    //- Additional, temporary fields for the solver

        //- Darcy velocity which is the 
        autoPtr<volVectorField> UDarcy_;

        //- Continuous Darcy flux 
        autoPtr<surfaceScalarField> phiDarcy_;

        //- A field that contains a fixed density throughout 
        //- the iteration. In other words, this is a temporary container 
        //- for the current density. 
        volScalarField fixedRho_;

        //- The temperature correction term for the Boussinesq approximation
        volScalarField rhok_;

    //- Additional dimensioned constants for the Boussinesq approximation

        //- Reference temperature
        autoPtr<dimensionedScalar> Tref_;

        //- Thermal expansion coefficient
        autoPtr<dimensionedScalar> beta_;

    //- Member variable that contains the compressibility indicator

        //- User defined (in phaseDict) switch to use incompressible treatment 
        //- with Boussinesq approximation. If true, the equation of state must be 
        //- rhoConst
        bool incompressibleTreatment_;

public:

    TypeName("onePhaseLegacy");

    //- Constructors

        onePhaseLegacy
        (
            Time& time,
            fvMesh& mesh, 
            customPimpleControl& pimple, 
            fv::options& fvOptions
        );

        //- Destructor
        virtual ~onePhaseLegacy(){}

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

        //- Update continuity error of the fluid
        void correctContErr();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace thermalHydraulicsModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
