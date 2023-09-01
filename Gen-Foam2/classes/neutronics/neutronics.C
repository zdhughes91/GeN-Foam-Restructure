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

#include "neutronics.H"
#include "zeroGradientFvPatchFields.H"
#include "fvmSup.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(neutronics, 0);
    defineRunTimeSelectionTable(neutronics, dictionary);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::neutronics::neutronics
(
    fvMesh& mesh
)
:
    IOdictionary
    (
        IOobject
        (
            "neutronicsProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    reactorState_
    (
        IOobject
        (
            "reactorState",
            mesh_.time().timeName(),
            "uniform",
            mesh_.time(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        )
    ),
    keff_(reactorState_.lookupOrDefault("keff",1.0)),
    pTarget_(reactorState_.lookupOrDefault("pTarget",1.0)),
    powerDensity_
    (
        IOobject
        (
            "powerDensity",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimPower/dimVol, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    secondaryPowerDenisty_
    (
        IOobject
        (
            "secondaryPowerDenisty",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimPower/dimVol, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    disp_
    (
        IOobject
        (
            "disp",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("d_zero", dimLength, vector::zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    initialResidual_(1.0),
    eigenvalueNeutronics_
    (
        IOdictionary::lookupOrDefault("eigenvalueNeutronics", false)
    ),
    liquidFuel_
    (
        mesh.time().controlDict().lookupOrDefault("liquidFuel", false)
    )
{
    Info << "Initial keff = " << keff_ << endl;
}

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::neutronics> Foam::neutronics::New
(
    fvMesh& mesh
)
{
    word modelName;

    IOdictionary dict
    (
        IOobject
        (
            "neutronicsProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dict.lookup("model") >> modelName;

    Info<< "Selecting neutronics model type " << modelName << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelName);

    if (!ctorPtr)
    {
        FatalErrorIn
        (
            "neutronics::New(const volScalarField&, "
            "const volVectorField&, basicThermo&)"
        )   << "Unknown neutronics model " << modelName
            << endl << endl
            << "Valid models types are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }
    return
        autoPtr<neutronics>
        (
            ctorPtr(mesh)
        );
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::neutronics::deformMesh
(
    const meshToMesh& TMToNeutro,
    const volVectorField& dispOrig
)
{
    const volPointInterpolation& neutroMeshPointInterpolation =
        volPointInterpolation::New(mesh_);

    tmp<pointVectorField> neutroPointsDisplacementOld =
        neutroMeshPointInterpolation.interpolate(disp_);

    disp_ *= 0.0;
    TMToNeutro.mapSrcToTgt(dispOrig, plusEqOp<vector>(), disp_);
    disp_.correctBoundaryConditions();

    tmp<pointVectorField> neutroPointsDisplacement =
        neutroMeshPointInterpolation.interpolate(disp_);

    tmp<pointField> displacedPoints =
        mesh_.points()
    +   neutroPointsDisplacement->internalField()
    -   neutroPointsDisplacementOld->internalField();

    mesh_.movePoints(displacedPoints);
}

// ************************************************************************* //
