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

#include "SNNeutronics.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "coordinateSystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(SNNeutronics, 0);

    addToRunTimeSelectionTable
    (
        neutronics,
        SNNeutronics,
        dictionary
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SNNeutronics::SNNeutronics
(
    fvMesh& mesh
)
:
    neutronics(mesh),//SNNeutronics is derived from neutronics
    xs_(mesh),
    quadratureSet_
    (
        IOobject
        (
            "quadratureSet",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    discreteDirections_(8*quadratureSet_.lookupOrDefault<int>("discreteDirectionsOct", 1)),
    directionVersors_(discreteDirections_),
    directionWeights_(discreteDirections_),
    legendreMatrices_(xs_.legendreMoments()),
    facePhi_("facePhi",vector(1, 0, 0) & mesh_.Sf()),
    facePhis_(discreteDirections_),
    angularFlux_(xs_.energyGroups()),
    prec_(xs_.precGroups()),
    precStar_(xs_.precGroups()),
    defaultFlux_
    (
        IOobject
        (
            "defaultFlux",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    defaultPrec_
    (
        IOobject
        (
            "defaultPrec",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimensionSet(0,-3,0,0,0,0,0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    flux_(xs_.energyGroups()),
    oneGroupFlux_
    (
        IOobject
        (
            "oneGroupFlux",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        defaultFlux_
    ),
    neutroSource_
    (
        IOobject
        (
            "neutroSource",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimensionSet(0,-3,-1,0,0,0,0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    delayedNeutroSource_
    (
        IOobject
        (
            "delayedNeutroSource",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimensionSet(0,-3,-1,0,0,0,0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    TFuel_
    (
        IOobject
        (
            "TFuel",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimTemperature, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    TClad_
    (
        IOobject
        (
            "TClad",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimTemperature, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    TCool_
    (
        IOobject
        (
            "TCool",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimTemperature, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    rhoCool_
    (
        IOobject
        (
            "rhoCool",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimTemperature, SMALL),
        zeroGradientFvPatchScalarField::typeName
    ),
    keff0_(0.0),
    keff1_(0.0),
    keff2_(0.0),
    aitkenIterCounter_(0),
    aitkenIterNo_(0)
{
    #include "readQuadratureSet.H"
    #include "calcLegendreMatrices.H"
    #include "createNeutronicsFieldsSN.H"
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::SNNeutronics::~SNNeutronics()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::SNNeutronics::getCouplingFieldRefs
(
    const objectRegistry& src,
    const meshToMesh& neutroToFluid
)
{
    #include "defaultGetCouplingFieldRefs.H"
}

void Foam::SNNeutronics::interpolateCouplingFields
(
    const meshToMesh& neutroToFluid
)
{
    #include "defaultInterpolateCouplingFields.H"
}

void Foam::SNNeutronics::correct
(
    scalar& residual, 
    label couplingIter
) 
{
    #include "solveNeutronicsSN.H"
}

// ************************************************************************* //
