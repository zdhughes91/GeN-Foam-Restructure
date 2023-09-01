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
#include "FSDragFactor.H"
#include "myOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(FSDragFactor, 0);

    typedef HashTable
    <
        autoPtr<FSDragCoefficientModel>,
        word,
        word::hash
    > fdTable;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FSDragFactor::FSDragFactor
(
    const FSPair& pair,
    const dictionary& dict
)
:
    IOdictionary
    (
        IOobject
        (
            "FSDragFactor."+typeName+"."+pair.name()+"."+dict.dictName(),
            pair.mesh().time().timeName(),
            pair.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dict
    ),
    pair_(pair),
    cells_(0),
    isotropic_(true),
    alphas_(pair_.structureRef()),
    alphaN_(pair_.fluidRef().normalized()),
    rho_(pair_.fluidRef().rho()),
    Dh_(pair_.fluidRef().Dh()),
    U_(pair_.fluidRef().U()),
    magU_(pair_.fluidRef().magU()),
    halfAlphaRhoMagUByDh_(0),
    validX_(false),
    validY_(false),
    validZ_(false)
{
    //- Set cells_ in which this model exists. This is done based on the
    //  regions in which this model is defined, with the region names being
    //  in the dict name
    wordList regions(myOps::split<word>(dict.dictName(), ':'));
    forAll(regions, i)
    {
        const labelList& cells(pair_.structureRef().cellLists()[regions[i]]);
        forAll(cells, j)
        {
            cells_.append(cells[j]);
        }
    }

    if (!this->found("type"))
        isotropic_ = false;

    if (isotropic_)
    {
        fdPtr_.reset
        (
            FSDragCoefficientModel::New
            (
                pair_,
                *this,
                pair_.mesh()
            )
        );
    }
    else
    {
        word X("localX");
        word Y("localY");
        word Z("localZ");
        if (this->found(X))
        {
            const dictionary& XDict(this->subDict(X));
            fdPtrs_.set
            (
                X,
                FSDragCoefficientModel::New
                (
                    pair_,
                    XDict,
                    pair_.mesh()
                )
            );
            fdPtrs_[X]->cmpt() = 0;
            validX_ = true;
        }
        if (this->found(Y))
        {
            const dictionary& YDict(this->subDict(Y));
            fdPtrs_.set
            (
                Y,
                FSDragCoefficientModel::New
                (
                    pair_,
                    YDict,
                    pair_.mesh()
                )
            );
            fdPtrs_[Y]->cmpt() = 1;
            validY_ = true;
        }
        if (this->found(Z))
        {
            const dictionary& ZDict(this->subDict(Z));
            fdPtrs_.set
            (
                Z,
                FSDragCoefficientModel::New
                (
                    pair_,
                    ZDict,
                    pair_.mesh()
                )
            );
            fdPtrs_[Z]->cmpt() = 2;
            validZ_ = true;
        }
        if (validX_ or validY_ or validZ_)
            halfAlphaRhoMagUByDh_ = scalarField(cells_.size(), 0);
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::FSDragFactor::correctField(volTensorField& Kd) const
{
    scalar minMagU
    (
        pair_.mesh().solutionDict().subDict("PIMPLE").lookupOrDefault<scalar>
        (
            "minMagU", 
            0
        )
    );
    if (isotropic_)
    {
        forAll(cells_, i)
        {
            const label& celli(cells_[i]);
            tensor& Kdi(Kd[celli]);
            scalar value
            (
                0.5/Dh_[celli]*
                (1.0-alphas_[celli])*
                rho_[celli]*
                max(magU_[celli], minMagU)
                *fdPtr_->value(celli)
            );
            Kdi[0] = value;
            Kdi[4] = value;
            Kdi[8] = value;
        }
    }
    else
    {
        forAll(cells_, i)
        {
            const label& celli(cells_[i]);
            halfAlphaRhoMagUByDh_[i] = 
                0.5/Dh_[celli]*
                (1.0-alphas_[celli])*
                rho_[celli]*
                max(magU_[celli], minMagU);
        }
        if (validX_)
        {
            const autoPtr<FSDragCoefficientModel>& fdXPtr
            (
                fdPtrs_["localX"]
            );
            forAll(cells_, i)
            {
                const label& celli(cells_[i]);
                tensor& Kdi(Kd[celli]);
                Kdi[0] = halfAlphaRhoMagUByDh_[i]*fdXPtr->value(celli);
            }
        }
        if (validY_)
        {
            const autoPtr<FSDragCoefficientModel>& fdYPtr
            (
                fdPtrs_["localY"]
            );
            forAll(cells_, i)
            {
                const label& celli(cells_[i]);
                tensor& Kdi(Kd[celli]);
                Kdi[4] = halfAlphaRhoMagUByDh_[i]*fdYPtr->value(celli);
            }
        }
        if (validZ_)
        {
            const autoPtr<FSDragCoefficientModel>& fdZPtr
            (
                fdPtrs_["localZ"]
            );
            forAll(cells_, i)
            {
                const label& celli(cells_[i]);
                tensor& Kdi(Kd[celli]);
                Kdi[8] = halfAlphaRhoMagUByDh_[i]*fdZPtr->value(celli);
            }
        }
    }
}

// ************************************************************************* //
