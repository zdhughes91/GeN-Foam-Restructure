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
#include "CachardLiquidFSHeatTransferCoefficient.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace FSHeatTransferCoefficientModels
{
    defineTypeNameAndDebug(CachardLiquid, 0);
    addToRunTimeSelectionTable
    (
        FSHeatTransferCoefficientModel, 
        CachardLiquid, 
        FSHeatTransferCoefficientModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FSHeatTransferCoefficientModels::CachardLiquid::CachardLiquid
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
    Twall_(pair.structureRef().Twall()),
    Tsat_(pair.mesh().lookupObject<volScalarField>("T.interface")),
    Tfluid_(pair.fluidRef().thermo().T()),
    alpha_(pair.mesh().lookupObject<volScalarField>("normalized.alpha.vapour")),
    p_(pair.mesh().lookupObject<volScalarField>("p")),
    rhoLiquid_(pair.fluidRef().rho()),
    rhoVapour_(pair.mesh().lookupObject<volScalarField>("thermo:rho.vapour")),
    kg_(pair.mesh().lookupObject<volScalarField>("kappa.vapour")),
    mug_(pair.mesh().lookupObject<volScalarField>("mu.vapour")),
    Dh_(pair.mesh().lookupObject<volScalarField>("Dh.vapour")),
    hGamma_
    (
        IOobject
        (
            "hGammaCachard",
            pair.mesh().time().timeName(),
            pair.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        pair.mesh(),
        dimensionedScalar("",dimensionSet(1,0,-3,-1,0,0,0),0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    sigmaSB_(5.670*1e-8), 
    gravity_(9.81),
    epsilonWall_(dict.get<scalar>("wallEmissivity")),
    epsilonLiq_(dict.get<scalar>("liquidEmissivity"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

/*
IMPLEMENTATION NOTES


*/

Foam::scalar Foam::FSHeatTransferCoefficientModels::CachardLiquid::value
(
    const label& celli
) const
{
    const scalar& Twi(Twall_[celli]);
    const scalar& Tsati(Tsat_[celli]);
    const scalar& Tli(Tfluid_[celli]);
    const scalar& alphai(alpha_[celli]);
    const scalar& pi(p_[celli]);
    const scalar& rhoLiqi(rhoLiquid_[celli]);
    const scalar& rhoVapi(rhoVapour_[celli]);
    const scalar& kgi(kg_[celli]);
    const scalar& mugi(mug_[celli]);
    const scalar DRi = 0.01; // Bundle Rod Diameter
    const scalar& Dhi(Dh_[celli]);

    // Radiation due to Phase 
    scalar dT(Twi-Tli);
    dT = (dT >= 0.0) ? max(dT, 1e-6) : min(dT, -1e-6);
    // Film Thickness for Rod Bundle (see TRACE)
    scalar deltai(Dhi/2*(pow(1+alphai*((4/Foam::constant::mathematical::pi)*pow(pi/max(DRi,1e-6),2)-1),0.5)-1));
    // Non dimensionnal Film Thickness (see TRACE)
    scalar deltastari(deltai*pow(rhoVapi*gravity_*(rhoLiqi-rhoVapi)/pow(max(mugi,1e-6),2),1.0/3));
    scalar hIniti(kgi/max(deltai,1e-6));    
    // Liquid Nusselt Number 
    scalar Nuwli(max(0,1.3*(0.268*pow(deltastari,0.77)-0.34)));
    scalar hwlPhase(hIniti*Nuwli*(Twi-Tsati)/dT);
    
    // HTC due to Radiation
    scalar a(max(epsilonLiq_*sqrt(1.0-alphai),1e-6));
    scalar b(1.0/epsilonWall_-1.0); 
    scalar hwlRad(sigmaSB_*(pow(Twi,2)+pow(Tsati,2))*(Twi+Tsati)/max(1/a+b,1e-6));

    // IOobject for the hwl Phase -> gives the mass rate change to MultiRegimeBoilingTRACECHF. 
    hGamma_[celli] = hwlPhase;

    return hwlPhase+hwlRad;

}

// ************************************************************************* //
