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
#include "NoKazimiFFHeatTransferCoefficient.H"
#include "mathematicalConstants.H"
#include "physicoChemicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace FFHeatTransferCoefficientModels
{
    defineTypeNameAndDebug(NoKazimi, 0);
    addToRunTimeSelectionTable
    (
        FFHeatTransferCoefficientModel, 
        NoKazimi, 
        FFHeatTransferCoefficientModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FFHeatTransferCoefficientModels::NoKazimi::NoKazimi
(
    const FFPair& pair,
    const dictionary& dict,
    const objectRegistry& objReg
)
:
    FFHeatTransferCoefficientModel
    (
        pair,
        dict,
        objReg
    ),
    alpha_(bulkFluid_.normalized()),
    rho_(bulkFluid_.rho()),
    kappa_(bulkFluid_.kappa()),
    Cp_(bulkFluid_.Cp()),
    Pr_(bulkFluid_.Cp()),
    magU_(bulkFluid_.magU()),
    p_(pair.mesh().lookupObject<volScalarField>("p")),
    iA_(pair.iA()),
    iT_(pair.iT()),
    coeff_
    (
        sqrt
        (
            // W/1000 g/mol->kg/mol
            bulkFluid_.thermo().W()().average().value()/1000.0/ 
            (
                2.0*constant::mathematical::pi*
                constant::physicoChemical::R.value()
            )
        )
    )
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::FFHeatTransferCoefficientModels::NoKazimi::value
(
    const label& celli
) const
{
    scalar alpha(max(alpha_[celli], 1e-9));
    const scalar& rho(rho_[celli]);
    const scalar& L(pair_.L()[celli]);
    return 
        min
        (
            iA_[celli]*kappa_[celli]/(alpha*(1.0-alpha))
        +   27.73*bulkFluid_[celli]*rho*magU_[celli]*Cp_[celli]*Pr_[celli],
            coeff_*sqr(rho*L)/(p_[celli]*sqrt(iT_[celli]))
        );
}


// ************************************************************************* //
