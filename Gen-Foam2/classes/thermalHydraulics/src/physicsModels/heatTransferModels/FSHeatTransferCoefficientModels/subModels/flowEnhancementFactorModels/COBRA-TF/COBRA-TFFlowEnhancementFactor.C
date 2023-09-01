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

#include "FSPair.H"
#include "COBRA-TFFlowEnhancementFactor.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace flowEnhancementFactorModels
{
    defineTypeNameAndDebug(COBRA_TF, 0);
    addToRunTimeSelectionTable
    (
        flowEnhancementFactorModel,
        COBRA_TF,
        flowEnhancementFactorModels
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::flowEnhancementFactorModels::COBRA_TF::COBRA_TF
(
    const FSPair& pair,
    const dictionary& dict,
    const objectRegistry& objReg
)
:
    flowEnhancementFactorModel
    (
        pair,
        dict,
        objReg
    ),
    otherFluidPtr_(nullptr)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::flowEnhancementFactorModels::COBRA_TF::value
(
    const label& celli
) const
{
    if (otherFluidPtr_ == nullptr)
    {
        HashTable<const fluid*> fluids(pair_.mesh().lookupClass<fluid>());
        otherFluidPtr_ = 
            (fluids[fluids.toc()[0]]->name() == pair_.fluidRef().name()) ?
            fluids[fluids.toc()[1]] : fluids[fluids.toc()[0]];
    }   

    const scalar& a(pair_.fluidRef().normalized()[celli]);
    const scalar& oa(otherFluidPtr_->normalized()[celli]);
    const scalar& Re1pi(pair_.Re()[celli]);
    scalar Re2pi
    (
        max
        (
            mag
            (
                a*pair_.fluidRef().rho()[celli]*pair_.fluidRef().U()[celli]
            +   oa*otherFluidPtr_->rho()[celli]*otherFluidPtr_->U()[celli]
            )*pair_.structureRef().Dh()[celli]/
            (
                    a*pair_.fluidRef().mu()[celli]
                +   oa*otherFluidPtr_->mu()[celli]
            ),
            10.0
        )
    );

    return min(max(pow(Re2pi/Re1pi, 0.8), 1.0), maxValue_);
}

// ************************************************************************* //
