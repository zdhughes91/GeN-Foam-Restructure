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
#include "AutruffeFFDragCoefficient.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace FFDragCoefficientModels
{
    defineTypeNameAndDebug(Autruffe, 0);
    addToRunTimeSelectionTable
    (
        FFDragCoefficientModel, 
        Autruffe, 
        FFDragCoefficientModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FFDragCoefficientModels::Autruffe::Autruffe
(
    const FFPair& pair,
    const dictionary& dict,
    const objectRegistry& objReg
)
:
    FFDragCoefficientModel
    (
        pair,
        dict,
        objReg
    ),
    vapour_
    (
        (pair.fluid1().isGas()) ? pair.fluid1() : pair.fluid2()
    )
{
    const fluid& fluid1(pair.fluid1());
    const fluid& fluid2(pair.fluid2());
    if 
    (
        !(fluid1.isLiquid() and fluid2.isGas()) and
        !(fluid2.isLiquid() and fluid1.isGas())
    )
    {
        FatalErrorInFunction
            << "The Autruffe model only works for liquid-gas systems. Set "
            << "the stateOfMatter entry in "
            << "phaseProperties." << fluid1.name() << "Properties and/or "
            << "phaseProperties." << fluid2.name() << "Properties) to "
            << "distinguish between gas and liquid"
            << exit(FatalError);
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::FFDragCoefficientModels::Autruffe::value
(
    const label& celli
) const
{
    const scalar& a(vapour_.normalized()[celli]);
    //- The coeff is 2.155 but, the FFDragFactor multiplies this value times
    //  0.5, so the coeff here is doubled
    return 
        4.31*max(pow((1.0-a)*(1.0+75.0*(1.0-a)),0.95),0.005)
        *(vapour_.rho()[celli]/pair_.rhoContinuous()[celli])
        *(pair_.DhDispersed()[celli]/pair_.DhContinuous()[celli]);
}


// ************************************************************************* //
