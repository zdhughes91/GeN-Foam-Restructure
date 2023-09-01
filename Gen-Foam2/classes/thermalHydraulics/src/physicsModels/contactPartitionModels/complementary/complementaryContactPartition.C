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
#include "complementaryContactPartition.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace contactPartitionModels
{
    defineTypeNameAndDebug(complementary, 0);
    addToRunTimeSelectionTable
    (
        contactPartitionModel,
        complementary,
        contactPartitionModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::contactPartitionModels::complementary::complementary
(
    const FSPair& pair,
    const dictionary& dict,
    const objectRegistry& objReg
)
:
    contactPartitionModel
    (
        pair,
        dict,
        objReg
    ),
    //complementaryPair_(nullptr),
    complementaryModel_(nullptr)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::contactPartitionModels::complementary::value
(
    const label& celli
) const
{
    if (complementaryModel_ == nullptr)
    {
        HashTable<const contactPartitionModel*> models
        (
            mesh_.lookupClass<contactPartitionModel>()
        );
        complementaryModel_ = 
        (
            (this->type() == models[models.toc()[0]]->type()) ?
            models[models.toc()[1]] : models[models.toc()[0]]
        );
    }

    return 1.0 - complementaryModel_->value(celli);
}

void Foam::contactPartitionModels::complementary::correctField
(
    volScalarField& f
)
{
    /*
    if (complementaryPair_ == nullptr)
    {
        HashTable<const FSPair*> pairs(mesh_.lookupClass<FSPair>());
        complementaryPair_ = 
        &(
            (pair_.name() == pairs[pairs.toc()[0]]->name()) ? 
            pairs[pairs.toc()[1]] : pairs[pairs.toc()[0]]
        );
    }
    f = 1.0-complementaryPair_->f();
    */
    if (complementaryModel_ == nullptr)
    {
        HashTable<const contactPartitionModel*> models
        (
            mesh_.lookupClass<contactPartitionModel>()
        );
        complementaryModel_ = 
        (
            (this->type() == models[models.toc()[0]]->type()) ?
            models[models.toc()[1]] : models[models.toc()[0]]
        );
    }

    forAll(mesh_.cells(), i)
    {
        f[i] = 1.0 - complementaryModel_->value(i);
    }
    f.correctBoundaryConditions();
}

// ************************************************************************* //
