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

#include "powerModel.H"
#include "structure.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(powerModel, 0);
    defineRunTimeSelectionTable
    (
        powerModel, 
        powerModels
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::powerModel::powerModel
(
    structure& structureRef,
    const dictionary& dicts
)
:
    IOdictionary
    (
        IOobject
        (
            typeName,
            structureRef.mesh().time().timeName(),
            structureRef.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dicts
    ),
    structure_(structureRef),
    mesh_(structureRef.mesh()),
    cellList_(0),
    cellField_(mesh_.cells().size(), 0),
    iA_
    (
        IOobject
        (
            IOobject::groupName("volumetricArea", typeName),
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar
        (
            "",
            dimArea/dimVol,
            0
        ),
        zeroGradientFvPatchScalarField::typeName
    ),
    alpha_
    (
        IOobject
        (
            IOobject::groupName("alpha", typeName),
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimless, 0),
        zeroGradientFvPatchScalarField::typeName
    )
{
    //- Set the fundamental properties that define a powerModel, its
    //  interfacial area and its volume fraction (alpha). This might be, in
    //  principle, different than the volume fraction of the regions over which
    //  the powerModel is defined. If a volumeFraction keyword is found within
    //  the powerModel dict within a region, set that as the powerModel alpha,
    //  otherwise use the structure alpha for that region
    const wordList& regions(this->toc());
    forAll(regions, i)
    {
        word region(regions[i]);
        const dictionary& regionDict(this->subDict(region));
        
        const labelList& regionCells(structure_.cellLists()[region]);
        forAll(regionCells, j)
        {
            label cellj(regionCells[j]);
            cellList_.append(cellj);
            cellField_[cellj] = 1.0;
        }
        if (regionDict.found("volumeFraction"))
        {
            scalar alpha(regionDict.get<scalar>("volumeFraction"));
            forAll(regionCells, j)
            {
                label cellj(regionCells[j]);
                alpha_[cellj] = alpha;
            }
        }
        else
        {
            const volScalarField& regionAlpha
            (
                structure_.alphaFields()[region]
            );
            forAll(regionCells, j)
            {
                label cellj(regionCells[j]);
                alpha_[cellj] = regionAlpha[cellj];
            }
        }
    } 

    iA_.correctBoundaryConditions();
    alpha_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::powerModel::setInterfacialArea()
{
    const wordList& regions(this->toc());
    forAll(regions, i)
    {
        word region(regions[i]);
        const dictionary& regionDict(this->subDict(region));
        scalar iA(regionDict.get<scalar>("volumetricArea"));
        
        const labelList& regionCells(structure_.cellLists()[region]);
        forAll(regionCells, j)
        {
            label cellj(regionCells[j]);
            iA_[cellj] = iA;
        }
    }
    iA_.correctBoundaryConditions();
}

/*
Foam::volScalarField& Foam::powerModel::initFieldInTable
(
    word name,
    dimensionedScalar value,
    bool readIfPresent,
    bool autoWrite
)
{
    if (!IOFields_.found(name))
    {
        IOFields_.set
        (
            name,
            new volScalarField
            (
                IOobject
                (
                    name,
                    mesh_.time().timeName(),
                    mesh_,
                    (
                        (readIfPresent) ? 
                        IOobject::READ_IF_PRESENT : IOobject::NO_READ
                    ),
                    (
                        (autoWrite) ?
                        IOobject::AUTO_WRITE : IOobject::NO_WRITE
                    )
                ),
                mesh_,
                dimensionedScalar(value.name(), value.dimensions(), 0.0),
                zeroGradientFvPatchScalarField::typeName
            )
        );
    }
    volScalarField& field(*IOFields_[name]);
    
    IOobject fieldHeader
    (
        name,
        mesh_.time().timeName(),
        mesh_,
        IOobject::NO_READ
    );
    if (!fieldHeader.typeHeaderOk<volScalarField>(true))
    {
        forAll(cellList_, i)
        {
            field[cellList_[i]] = value.value();
        }
        field.correctBoundaryConditions();
    }

    return field;
}
*/

// ************************************************************************* //
