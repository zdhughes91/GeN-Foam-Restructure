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

#include "FFPair.H"
#include "FFDragFactor.H"
#include "myOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(FFDragFactor, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FFDragFactor::FFDragFactor
(
    const FFPair& pair,
    const dictionary& dict
)
:
    IOdictionary
    (
        IOobject
        (
            "FFDragFactor."+typeName+"."+pair.name()+"."+dict.dictName(),
            pair.mesh().time().timeName(),
            pair.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dict
    ),
    pair_(pair),
    coeffReg_(*this),
    rhoc_(pair_.rhoContinuous()),
    Dhd_(pair_.DhDispersed()),
    magUr_(pair_.magUr())
{
    fdPtr_.reset
    (
        FFDragCoefficientModel::New
        (
            pair_,
            *this,
            coeffReg_
        )
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::FFDragFactor::correctField(volScalarField& Kd) const
{
    const volScalarField& a1(pair_.fluid1());
    const volScalarField& a2(pair_.fluid2());
    forAll(pair_.mesh().cells(), i)
    {
        Kd[i] = 
            (a1[i]+a2[i])*(0.5/Dhd_[i])*rhoc_[i]*magUr_[i]*fdPtr_->value(i);
    }
}


// ************************************************************************* //
