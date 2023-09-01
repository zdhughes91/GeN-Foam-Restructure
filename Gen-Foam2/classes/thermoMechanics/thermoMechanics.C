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

#include "thermoMechanics.H"
#include "zeroGradientFvPatchFields.H"
#include "fvmSup.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(thermoMechanics, 0);
    defineRunTimeSelectionTable(thermoMechanics, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermoMechanics::thermoMechanics
(
    fvMesh& mesh
)
:
    IOdictionary
    (
        IOobject
        (
            "thermoMechanicalProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    meshDisp_
    (
        IOobject
        (
            "meshDisp",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector("", dimLength, vector::zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    initialResidual_(1.0)
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::thermoMechanics> Foam::thermoMechanics::New
(
    fvMesh& mesh
)
{
    word modelName;

    {
        IOdictionary dict
        (
            IOobject
            (
                "thermoMechanicalProperties",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        dict.lookup("model") >> modelName;
    }

    Info<< "Selecting thermoMechanics model type " << modelName << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelName);

    if (!ctorPtr)
    {
        FatalErrorIn
        (
            "thermoMechanics::New(const volScalarField&, "
            "const volVectorField&, basicThermo&)"
        )   << "Unknown thermoMechanics model " << modelName
            << endl << endl
            << "Valid models types are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }
    return
    autoPtr<thermoMechanics>
    (
        ctorPtr(mesh)
    );


}


// ************************************************************************* //
