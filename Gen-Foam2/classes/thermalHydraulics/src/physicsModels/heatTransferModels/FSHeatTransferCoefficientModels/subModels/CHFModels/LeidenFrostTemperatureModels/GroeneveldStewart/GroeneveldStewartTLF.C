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
#include "GroeneveldStewartTLF.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace LeidenFrostTemperatureModels
{
    defineTypeNameAndDebug(GroeneveldStewart, 0);
    addToRunTimeSelectionTable
    (
        TLFModel,
        GroeneveldStewart,
        LeidenFrostTemperatureModels
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LeidenFrostTemperatureModels::GroeneveldStewart::GroeneveldStewart
(
    const FSPair& pair,
    const dictionary& dict,
    const objectRegistry& objReg
)
:
    TLFModel
    (
        pair,
        dict,
        objReg
    ),
    criticalPressure_(dict.get<scalar>("criticalPressure"))
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::LeidenFrostTemperatureModels::GroeneveldStewart::value
(
    const label& celli
) const
{
    const scalar& pli(p_[celli]);
    const scalar& Tsati(Tsat_[celli]);

    scalar Tmin(0);
    if (pli<9*1e6) // Correlation from GroeneveldStewart, valid for pressure P<9 MPa
    {
        Tmin = 557.85+44.1*pli*(1e-6)-3.72*pow(pli*1e-6,2);
    }
    else // Ramp up to critical pressure 
    {
        scalar DeltaTmin(557.85+44.1*9-3.72*pow(9.0,2)-Tsati);
        Tmin = Tsati+(criticalPressure_-pli)/(criticalPressure_-9*1e6)*DeltaTmin;  
    }

    return Tmin;

}

// ************************************************************************* //
