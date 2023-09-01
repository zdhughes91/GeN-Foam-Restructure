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
#include "BrowningPotterSaturation.H"
#include "addToRunTimeSelectionTable.H"
#include "phaseChangeModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace saturationModels
{
    defineTypeNameAndDebug(BrowningPotter, 0);
    addToRunTimeSelectionTable
    (
        saturationModel,
        BrowningPotter,
        saturationModels
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::saturationModels::BrowningPotter::
BrowningPotter
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
    p_(pcm.mesh().lookupObject<volScalarField>("p"))
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::saturationModels::BrowningPotter::valuePSat
(
    const label& celli
) const
{
    return exp(valueLnPSat(celli));
}

Foam::scalar Foam::saturationModels::BrowningPotter::valuePSatPrime
(
    const label& celli
) const
{
    const scalar& T(iT_[celli]);
    return valuePSat(celli)*(-0.4672/T + 12633.37/sqr(T));
}

Foam::scalar Foam::saturationModels::BrowningPotter::valueLnPSat
(
    const label& celli
) const
{
    const scalar& T(iT_[celli]);
    //- The + log(1e6) is to have p in Pa rather than MPa
    return 
        (11.9463 - 12633.37/T - 0.4672*log(T)) + log(1e6); 
}

Foam::scalar Foam::saturationModels::BrowningPotter::valueTSat
(
    const label& celli
) const
{
    return 
        923840.0/
        (
        -   11275 
        +   Foam::sqrt
            (
                127125625 + 1847680*
                (
                    7.8270 
                -   log(p_[celli]/1e6) // p in MPa
                ) 
            )
        );
}

// ************************************************************************* //
