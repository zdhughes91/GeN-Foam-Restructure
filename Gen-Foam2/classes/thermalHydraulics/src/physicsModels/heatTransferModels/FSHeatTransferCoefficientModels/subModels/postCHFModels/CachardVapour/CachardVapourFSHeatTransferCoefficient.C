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
#include "CachardVapourFSHeatTransferCoefficient.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace FSHeatTransferCoefficientModels
{
    defineTypeNameAndDebug(CachardVapour, 0);
    addToRunTimeSelectionTable
    (
        FSHeatTransferCoefficientModel, 
        CachardVapour, 
        FSHeatTransferCoefficientModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FSHeatTransferCoefficientModels::CachardVapour::CachardVapour
(
    const FSPair& pair,
    const dictionary& dict,
    const objectRegistry& objReg
)
:
    FSHeatTransferCoefficientModel
    (
        pair,
        dict,
        objReg
    ),
    p_(pair.mesh().lookupObject<volScalarField>("p")),
    alpha_(pair.fluidRef().normalized()),
    kg_(pair.mesh().lookupObject<volScalarField>("kappa.vapour")),
    Dh_(pair.fluidRef().Dh()),
    gravity_(9.81)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

/*
IMPLEMENTATION NOTES


*/

Foam::scalar Foam::FSHeatTransferCoefficientModels::CachardVapour::value
(
    const label& celli
) const
{
    const scalar& kgi(kg_[celli]);
    const scalar& Dhi(Dh_[celli]);
    const scalar DRi = 0.01;
    const scalar& alphai(alpha_[celli]);
    const scalar& pi(p_[celli]);


    // Film Thickness delta 
    scalar deltai(Dhi/2*(pow(1+alphai*((4/Foam::constant::mathematical::pi)*pow(pi/max(DRi,1e-6),2)-1),0.5)-1));
    scalar hCachard = 2*kgi/max(deltai,1e-6);

    return hCachard;
}

// ************************************************************************* //
