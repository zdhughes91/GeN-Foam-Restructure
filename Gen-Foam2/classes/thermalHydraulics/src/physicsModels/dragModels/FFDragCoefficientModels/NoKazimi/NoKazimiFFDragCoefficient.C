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
#include "NoKazimiFFDragCoefficient.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace FFDragCoefficientModels
{
    defineTypeNameAndDebug(NoKazimi, 0);
    addToRunTimeSelectionTable
    (
        FFDragCoefficientModel, 
        NoKazimi, 
        FFDragCoefficientModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FFDragCoefficientModels::NoKazimi::NoKazimi
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
            << "The NoKazimi model only works for liquid-gas systems. Set "
            << "the stateOfMatter entry in "
            << "phaseProperties." << fluid1.name() << "Properties and/or "
            << "phaseProperties." << fluid2.name() << "Properties) to "
            << "distinguish between gas and liquid"
            << exit(FatalError);
    }

    //- Compute A_
    scalar P(dict.get<scalar>("pinPitch"));
    scalar D(dict.get<scalar>("pinDiameter"));
    A_ = 
        (4.0*constant::mathematical::pi)/
        (D*(2.0*sqrt(3.0)*sqr(P/D)-constant::mathematical::pi));
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::FFDragCoefficientModels::NoKazimi::value
(
    const label& celli
) const
{
    const scalar& a(vapour_.normalized()[celli]);
    //- The coeff is 0.0025 but, the FFDragFactor multiplies this value times
    //  0.5, so the coeff here is doubled. Also, this model does not use the
    //  dispersed phase hydraulic diameter as char. dimension, rather, it
    //  uses its own char. dimension computed via the A_ thing. For this reason
    //  I need to multply the coeff here by the dispersed phase dimension, to
    //  cancel it out
    return
        pair_.DhDispersed()[celli]*A_*min(a, 0.99)*min((1.0-a)/0.043, 1.0)*
        0.005*(1.0 + 234.3*pow((1.0-sqrt(a)), 1.15));
}


// ************************************************************************* //
