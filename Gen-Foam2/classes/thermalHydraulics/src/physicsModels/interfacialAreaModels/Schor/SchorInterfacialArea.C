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
#include "SchorInterfacialArea.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace interfacialAreaModels
{
    defineTypeNameAndDebug(Schor, 0);
    addToRunTimeSelectionTable
    (
        interfacialAreaModel,
        Schor,
        interfacialAreaModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfacialAreaModels::Schor::Schor
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
        (pair.fluid1().isGas()) ?
        pair.fluid1()
        :
        pair.fluid2()
    ),
    D_(dict.get<scalar>("pinDiameter")),
    P_(dict.get<scalar>("pinPitch")),
    PD_(P_/D_),
    A_(2.0*Foam::sqrt(3.0)*Foam::sqr(PD_)),
    PI_(3.1415927),
    alpha1_(0.55),
    alpha2_(0.65),
    deltaAlpha21_(alpha2_-alpha1_),
    alpha3_(dict.getOrDefault<scalar>("cutoffAlpha", 0.957)),
    iA1_((4.0/D_)*Foam::sqrt(PI_*alpha1_/(3.0*(A_-PI_)))),
    iA2_
    (
        (4.0/D_)*Foam::sqrt(PI_)
        /
        (
            3.0*(A_-PI_)
        )*
        Foam::sqrt
        (
            (1.0-alpha2_)*A_ + PI_*alpha2_
        )
    ),
    iA3_
    (
        dict.getOrDefault<scalar>("minInterfacialAreaAtLargeAlpha", 0)
    )
{
    if 
    (
        !(pair.fluid1().isLiquid() and pair.fluid2().isGas())
    and !(pair.fluid2().isLiquid() and pair.fluid1().isGas())
    )
    {
        FatalErrorInFunction
            << "Shor model only work for gas-liquid systems! Set the "
            << "stateOfMatter keyword in the fluid properties dictionaries to"
            << "specify the stateOfMatter of each fluid"
            << exit(FatalError);
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::interfacialAreaModels::Schor::value
(
    const label& celli
) const
{
    const scalar& alpha(vapour_.normalized()[celli]);
    scalar scaleFactor(vapour_[celli]/max(alpha, 1e-9));
    return 
        scaleFactor*
        (
            (alpha < alpha1_) ?
            (4.0/D_)*Foam::sqrt(PI_*alpha/(3.0*(A_-PI_)))
            :
            (
                (alpha < alpha2_) ?
                ((alpha2_-alpha)*iA1_+(alpha-alpha1_)*iA2_)/deltaAlpha21_
                :
                max
                (
                    (4.0/D_)*Foam::sqrt(PI_)
                    /
                    (
                        3.0*(A_-PI_)
                    )*
                    Foam::sqrt
                    (
                        (1.0-alpha)*A_ + PI_*alpha
                    )*
                    min
                    (
                        (1-alpha)/(1.0-alpha3_),
                        1.0
                    ),
                    iA3_/scaleFactor
                )
            )
        );
}

// ************************************************************************* //
