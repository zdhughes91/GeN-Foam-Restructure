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

#include "virtualMass.H"
#include "FFPair.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(virtualMass, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::virtualMass::virtualMass
(
    const FFPair& pair,
    const dictionary& dict
)
:
    IOdictionary
    (
        IOobject
        (
            typeName,
            pair.mesh().time().timeName(),
            pair.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dict
    ),
    mesh_(pair.mesh()),
    pair_(pair),
    fluid1_(pair.fluid1()),
    fluid2_(pair.fluid2()),
    U1_(fluid1_.U()),
    U2_(fluid2_.U()),
    phi1_(fluid1_.phi()),
    phi2_(fluid2_.phi()),
    VmPtr_
    (
        virtualMassCoefficientModel::New
        (
            pair,
            dict,
            mesh_
        )
    ),
    Vm_
    (
        IOobject
        (
            IOobject::groupName("Vm", pair.name()),
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
    )
{
    VmForces_.set
    (
        U1_.name(),
        fvVectorMatrix(fluid1_.U(), dimMass*dimLength/dimTime/dimTime)
    );
    VmForces_.set
    (
        U2_.name(),
        fvVectorMatrix(fluid2_.U(), dimMass*dimLength/dimTime/dimTime)
    );
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::virtualMass::correct()
{
    //- Store Vm in alphaRhoVm_
    VmPtr_->correctField(Vm_);

    //- Actually update alphaRhoVm
    bool virtualMassUsesMixtureDensity
    (
        pair_.pimple().dict().lookupOrDefault<bool>
        (
            "virtualMassUsesMixtureDensity", 
            false
        )
    );
    volScalarField alphaRhoVm(Vm_*pair_.alphaDispersed());
    if (virtualMassUsesMixtureDensity)
        alphaRhoVm *= 
            (
                fluid1_.rho()*fluid1_.normalized()
            +   fluid2_.rho()*fluid2_.normalized()
            );
    else
        alphaRhoVm *= pair_.rhoContinuous();

    //-
    VmForces_[U1_.name()] = 
        alphaRhoVm*
        (
            fvm::ddt(U1_)
        +   fvm::div(phi1_, U1_)
        +   fvm::SuSp(-fvc::div(phi1_), U1_)
        -   (
                fvc::ddt(U2_)
            +   fvc::div(phi2_, U2_)
            -   fvc::div(phi2_)*U2_
            )
        );

    VmForces_[U2_.name()] = 
        alphaRhoVm*
        (
            fvm::ddt(U2_)
        +   fvm::div(phi2_, U2_)
        +   fvm::SuSp(-fvc::div(phi2_), U2_)
        -   (
                fvc::ddt(U1_)
            +   fvc::div(phi1_, U1_)
            -   fvc::div(phi1_)*U1_
            )
        );
}

// ************************************************************************* //
