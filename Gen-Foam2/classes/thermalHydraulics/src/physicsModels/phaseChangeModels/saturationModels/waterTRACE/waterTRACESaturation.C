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

#include "FFPair.H"
#include "waterTRACESaturation.H"
#include "addToRunTimeSelectionTable.H"
#include "phaseChangeModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace saturationModels
{
    defineTypeNameAndDebug(waterTRACE, 0);
    addToRunTimeSelectionTable
    (
        saturationModel,
        waterTRACE,
        saturationModels
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::saturationModels::waterTRACE::
waterTRACE
(
    const phaseChangeModel& pcm,
    const dictionary& dict, 
    const objectRegistry& objReg
)
:
    saturationModel
    (
        pcm,
        dict,
        objReg
    ),
    iT_(pcm.pair().iT()),
    p_(pcm.mesh().lookupObject<volScalarField>("p")),
    /*
    // GeN-Foam Model
    A_(1e6*2.590718143628726e-10),
    B_(273.159),
    C_(4.247368421052632),
    */
    // Constant of water 
    Rv_(461.4975)  // J/kg/K

{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::saturationModels::waterTRACE::valuePSat
(
    const label& celli
) const
{
    const scalar& T(iT_[celli]);
    /*
    // ----- Stefan Radman Version - Interpolation NIST  ----- //
    return A_*pow(T-B_, C_);
    */

    // ----- TRACE Version - different formulas according to temperature ----- //
     if (T<370.4251)
    {
        scalar Ts = T;
        if (T<=273.15)
        {
            Ts = 273.15;
        }
        const scalar AAP(24821.0);
        const scalar BAP(338.0);
        const scalar CAP(-5.3512);
        const scalar DAP(20.387);
        scalar ps(AAP*pow(Ts/BAP,CAP)*exp(DAP*(Ts-BAP)/Ts));
        return ps; 
    }
    else if (T<609.62462615967)
    {
        const scalar AB(117.8);
        const scalar BB(1e-5);
        const scalar CB(0.223);
        const scalar DB(255.2);
        scalar ps(BB*pow((T-DB)/AB,1.0/CB));
        return ps;
    }
    else if (T<647.3)
    {
        const scalar AC(7.2166948490268e11);
        const scalar BC(-8529.6481905883);
        const scalar CC(1166669.3278328);
        scalar ps(AC*exp((BC+CC/T)/T));
        return ps;
    }
    else
    {
        const scalar AD(22.12e6);
        const scalar BD(7.6084087);
        const scalar CD(4924.9229);
        scalar ps(AD*exp(BD-CD/T));
        return ps;
    }
}

Foam::scalar Foam::saturationModels::waterTRACE::valuePSatPrime
(
    const label& celli
) const
{
    const scalar& T(iT_[celli]);
    const scalar& ps(p_[celli]);

    /*
    // ----- Stefan Radman Version - Interpolation NIST  ----- //
    return C_*A_*pow(T-B_, C_-1.0);
    */

    if (T<370.4251)
    {
        scalar Ts = T;
        if (T<=273.15)
        {
            Ts = 273.15;
        }
        const scalar AA(3180619.59);
        const scalar BA(2470.2120);
        scalar hs(AA-BA*Ts);
        scalar slope(hs*ps/Rv_/sqr(Ts)); 
        return slope; 
    }
    else if (T<609.62462615967)
    {
        const scalar AB(0.223);
        const scalar BB(255.2);
        scalar slope(ps/(AB*(T-BB)));
        return slope;
    }
    else if (T<647.3)
    {
        const scalar AC(-8529.6481905883);
        const scalar BC(2333338.6556656);
        scalar slope(-ps*(AC+BC/T)/sqr(T));
        return slope;
    }
    else
    {
        const scalar AD(2.0304886238506e-4);
        scalar slope(ps/(AD*sqr(T)));
        return slope;
    }
}

Foam::scalar Foam::saturationModels::waterTRACE::valueLnPSat
(
    const label& celli
) const
{
    return log(valuePSat(celli));
}

Foam::scalar Foam::saturationModels::waterTRACE::valueTSat
(
    const label& celli
) const
{
    const scalar& pi(p_[celli]); 
    /*
    // ----- Stefan Radman Version - Interpolation NIST  ----- //
    scalar TsRAD(pow(pi/A_, 1.0/C_) + B_ );
    //Info << "Tsat Interpolation NIST --" << TsRAD << endl;
    */

    //-------------------------------------------------------------//
    
    // ----- TRACE Version - different formulas according to pressure ----- //

    if (pi<90.56466*1000)
    {
        scalar ps = pi;
        if (pi<610.8)
        {
            ps = 610.8;
        }
        const scalar AA(-2263.0);
        const scalar BA(0.434);
        const scalar CA(100000.0);
        const scalar DA(6.064);
        const scalar AAP(24821.0);
        const scalar BAP(338.0);
        const scalar CAP(-5.3512);
        const scalar DAP(20.387);
        const scalar AAH(3180619.59);
        const scalar BAH(2470.2120);

        scalar Tsapprox(0);
        scalar psapprox (0);
        scalar hsapprox(0);
        scalar Ts(AA/(BA*log(ps/CA)-DA));
        for (int step=0; step<2; step++)
        {
            Tsapprox=Ts;
            psapprox=AAP*pow(Tsapprox/BAP,CAP)*exp(DAP*(Tsapprox-BAP)/Tsapprox);
            hsapprox=AAH - BAH*Tsapprox;
            Ts=Tsapprox/(1-Rv_*Tsapprox/hsapprox*log(ps/psapprox));
        }
        return Ts; 
    }
    else if (pi<13.969971285053*1e6)
    {
        const scalar AB(117.8);
        const scalar BB(1.0*pow(10.0,-5));
        const scalar CB(0.223);
        const scalar DB(255.2);
        scalar Ts(AB*pow(BB*pi,CB)+DB);
        return Ts;
    }
    else if (pi<22.12*1e6)
    {
        const scalar AC(4264.8240952941);
        const scalar BC(13666986.708428);
        const scalar CC(1166669.3278328);
        const scalar DC(27.304833093884);
        scalar Ts((AC+sqrt(-BC+CC*log(pi)))/(DC-log(pi)));
        return Ts;
    }
    else
    {
        const scalar AD(4924.9229);
        const scalar BD(24.520401);
        scalar Ts(AD/(BD-log(pi)));
        return Ts;
    }
}

// ************************************************************************* //
