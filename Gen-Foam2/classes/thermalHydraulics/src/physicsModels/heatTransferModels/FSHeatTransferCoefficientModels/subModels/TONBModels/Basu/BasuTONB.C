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
#include "BasuTONB.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace TONBModels
{
    defineTypeNameAndDebug(Basu, 0);
    addToRunTimeSelectionTable
    (
        TONBModel,
        Basu,
        TONBModels
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::TONBModels::Basu::Basu
(
    const FSPair& pair,
    const dictionary& dict,
    const objectRegistry& objReg
)
:
    TONBModel
    (
        pair,
        dict,
        objReg
    ),
    sigma_(dict.get<scalar>("surfaceTension")),
    contactAngle_
    (
        dict.get<scalar>("contactAngle")*
        2.0*constant::mathematical::pi/360.0
    ),
    A_
    (
        2.0*sigma_/sqr((1.0-exp(-pow(contactAngle_, 3)-0.5*contactAngle_)))
    )
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::TONBModels::Basu::value
(
    const label& celli,
    const scalar& htc2pFCi
) const
{
    const scalar& Tli(pair_.fluidRef().T()[celli]);
    const scalar& Tsati(Tsat_[celli]);

    this->setPtrs();

    scalar deltaTONBsat
    (
        A_*htc2pFCi*Tsati/
        (
            otherFluidPtr_->rho()[celli]*
            mag((*LPtr_)[celli])*
            pair_.fluidRef().kappa()[celli]
        )
    );

    return  
    (
        Tli + 0.25*
        sqr
        (
            sqrt(deltaTONBsat)
        +   sqrt(max(deltaTONBsat+4.0*(Tsati-Tli), 0.0))
        )
    );
}

// ************************************************************************* //
