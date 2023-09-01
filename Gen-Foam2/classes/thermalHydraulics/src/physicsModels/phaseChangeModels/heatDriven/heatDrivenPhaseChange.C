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

#include "heatDrivenPhaseChange.H"
#include "FFPair.H"
#include "addToRunTimeSelectionTable.H"
#include "myOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseChangeModels
{
    defineTypeNameAndDebug(heatDrivenPhaseChange, 0);
    addToRunTimeSelectionTable
    (
        phaseChangeModel,
        heatDrivenPhaseChange,
        phaseChangeModels
    );
}
}

const Foam::Enum
<
    Foam::phaseChangeModels::heatDrivenPhaseChange::mode
>
Foam::phaseChangeModels::heatDrivenPhaseChange::modeNames_
(
    {
        { 
            mode::conductionLimited, 
            "conductionLimited" 
        },
        { 
            mode::onePhaseDriven, 
            "onePhaseDriven" 
        },
        { 
            mode::twoPhaseDriven, 
            "twoPhaseDriven" 
        },
        { 
            mode::mixedDriven, 
            "mixedDriven" 
        }
    }
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseChangeModels::heatDrivenPhaseChange::heatDrivenPhaseChange
(
    FFPair& pair,
    const dictionary& dict
)
:
    phaseChangeModel
    (
        pair,
        dict
    ),
    mode_(modeNames_.get(this->get<word>("mode"))),
    drivingPhaseName_
    (
        (mode_ == heatDrivenPhaseChange::mode::onePhaseDriven) ?
        this->get<word>("drivingPhase") : ""
    )
{
    if (drivingPhaseName_ != "")
    {
        if  
        (
            drivingPhaseName_ != fluid1_.name()
        and drivingPhaseName_ != fluid2_.name()
        )
        {
            FatalErrorInFunction
            << "phaseChangeModel: phase " << drivingPhaseName_ << " not found!"
            << exit(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::phaseChangeModels::heatDrivenPhaseChange::correctInterfacialDmdt() 
{
    //- Interfacial mass transfers for each side
    volScalarField dmdtI1i
    (
        htc1_*iA_*(T1_-iT_)/L_
    );
    volScalarField dmdtI2i
    (
        htc2_*iA_*(T2_-iT_)/L_
    );

    switch (mode_)
    {
        case heatDrivenPhaseChange::mode::conductionLimited : 
            dmdtI_ = dmdtI1i + dmdtI2i;
            break;

        case heatDrivenPhaseChange::mode::twoPhaseDriven :
            dmdtI_ = posPart(dmdtI1i) + negPart(dmdtI2i);
            break;

        case heatDrivenPhaseChange::mode::onePhaseDriven :
            dmdtI_ = 
                (fluid1_.name() == drivingPhaseName_) ?
                dmdtI1i : dmdtI2i;      
            break;

        case heatDrivenPhaseChange::mode::mixedDriven :
            dmdtI_ = 
                (fluid1_.name() == drivingPhaseName_) ?
                dmdtI1i + negPart(dmdtI2i) :
                posPart(dmdtI1i) + dmdtI2i;
            break;

        default : break;
    }
}

// ************************************************************************* //
