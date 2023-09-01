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
#include "NoKazimiInterfacialArea.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace interfacialAreaModels
{
    defineTypeNameAndDebug(NoKazimi, 0);
    addToRunTimeSelectionTable
    (
        interfacialAreaModel,
        NoKazimi,
        interfacialAreaModels
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfacialAreaModels::NoKazimi::NoKazimi
(
    const FFPair& pair,
    const dictionary& dict,
    const objectRegistry& objReg
)
:
    interfacialAreaModel
    (
        pair,
        dict,
        objReg
    ),
    vapour_
    (
        (pair.fluid1().isGas()) ? pair.fluid1() : pair.fluid2()
    ),
    D_(dict.get<scalar>("pinDiameter")),
    P_(dict.get<scalar>("pinPitch")),
    PD_(P_/D_),
    A_
    (
        4.0*constant::mathematical::pi/
        (
            D_*(2*Foam::sqrt(3.0)*sqr(PD_)-constant::mathematical::pi)
        )
    )
{
    if 
    (
        !(pair.fluid1().isLiquid() and pair.fluid2().isGas())
    and !(pair.fluid2().isLiquid() and pair.fluid1().isGas())
    )
    {
        FatalErrorInFunction
            << "NoKazimi model only work for gas-liquid systems! Set the "
            << "stateOfMatter keyword in the fluid properties dictionaries to"
            << "specify the stateOfMatter of each fluid"
            << exit(FatalError);
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::interfacialAreaModels::NoKazimi::value
(
    const label& celli
) const
{
    const scalar& a(vapour_[celli]);
    const scalar& aN(vapour_.normalized()[celli]);
    return 
        A_*a*
        min
        (
            (1.0-aN)/0.043,
            1.0
        );
}

// ************************************************************************* //
