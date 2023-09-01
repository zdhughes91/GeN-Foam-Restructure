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
#include "waterLatentHeat.H"
#include "addToRunTimeSelectionTable.H"
#include "phaseChangeModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace latentHeatModels
{
    defineTypeNameAndDebug(water, 0);
    addToRunTimeSelectionTable
    (
        latentHeatModel,
        water,
        latentHeatModels
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::latentHeatModels::water::
water
(
    const phaseChangeModel& pcm,
    const dictionary& dict, 
    const objectRegistry& objReg
)
:
    latentHeatModel
    (
        pcm,
        dict,
        objReg
    )
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::latentHeatModels::water::value(const label& celli)
const
{
    //- The critical T of water is 647.3 K, but the fit spew out 681.718 as
    //  a limit temperature for the fit. Understandable, since the data range
    //  used for the fit is 0.01-350 *C, and that above that the latent heat 
    //  drops very, very rapidly to 0, which I do not want to deal with as I do
    //  not deal with supercritical transitions and shit like that. This is 
    //  fine. Latent heat here is in SI units (J/kg)

    scalar T(min(max(iT_[celli], 273.16), 680));
    return 1e6*(-2.44+0.82*log(681.718-T));
}

// ************************************************************************* //
