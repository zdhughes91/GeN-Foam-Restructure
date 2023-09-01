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

#include "twoPhaseDragMultiplierModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(twoPhaseDragMultiplierModel, 0);
    defineRunTimeSelectionTable
    (
        twoPhaseDragMultiplierModel, 
        twoPhaseDragMultiplierModels
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoPhaseDragMultiplierModel::twoPhaseDragMultiplierModel
(
    const fvMesh& mesh,
    const dictionary& dict,
    const objectRegistry& objReg
)
:
    IOdictionary
    (
        IOobject
        (
            typeName,
            mesh.time().timeName(),
            objReg,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dict
    ),
    mesh_(mesh),
    mFluidPtr_(nullptr),
    oFluidPtr_(nullptr),
    mKd_(nullptr),
    maxPhi2_
    (
        dict.lookupOrDefault<scalar>("maxPhi2", 1000)
    )
{
    //- Set ptr to multiplier fluid
    mFluidPtr_ = 
        &(
            mesh_.lookupObject<fluid>
            (
                word("alpha."+this->get<word>("multiplierFluid"))
            )
        );
    
    //- Set ptr to other fluid
    HashTable<const fluid*> fluids(mesh_.lookupClass<fluid>());
    wordList fluidNames(fluids.toc());
    oFluidPtr_ = (fluidNames[0] == "alpha."+mFluidPtr_->name()) ?
        fluids[fluidNames[1]] : fluids[fluidNames[0]];

    //- Set ptr to mKd
    const FSPair& pair
    (
        mesh_.lookupObject<FSPair>(mFluidPtr_->name()+".structure")
    );
    mKd_ = &(pair.Kd());
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::twoPhaseDragMultiplierModel::onePhase(const label& celli) const
{
    return (mFluidPtr_->normalized()[celli] > 0.9999);
}

//- Default value
Foam::scalar Foam::twoPhaseDragMultiplierModel::phi2
(
    const label& celli
) const
{
    return 1.0;
}

Foam::tensor Foam::twoPhaseDragMultiplierModel::value
(
    const label& celli
) const
{
    return
        tensor
        (
            min(maxPhi2_, phi2(celli))*((*mKd_)[celli])*
            sqr(max(mFluidPtr_->normalized()[celli], 1e-9))*
            mFluidPtr_->magU()[celli]
        );
}

void Foam::twoPhaseDragMultiplierModel::correctField(volTensorField& KdTotU)
{
    forAll(mesh_.cells(), i)
    {
        KdTotU[i] = this->value(i);
    }
    KdTotU.correctBoundaryConditions();
}

// ************************************************************************* //
