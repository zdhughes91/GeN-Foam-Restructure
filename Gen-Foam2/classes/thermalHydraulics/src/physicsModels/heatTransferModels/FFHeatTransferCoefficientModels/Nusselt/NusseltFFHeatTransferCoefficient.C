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
#include "NusseltFFHeatTransferCoefficient.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace FFHeatTransferCoefficientModels
{
    defineTypeNameAndDebug(Nusselt, 0);
    addToRunTimeSelectionTable
    (
        FFHeatTransferCoefficientModel, 
        Nusselt, 
        FFHeatTransferCoefficientModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FFHeatTransferCoefficientModels::Nusselt::Nusselt
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
    Re_(pair.Re()),
    kappa_(bulkFluid_.kappa()),
    Pr_(pair.PrContinuous()),//(bulkFluid_.Pr()),(bulkFluid_.Pr()),
    Dh_(pair.DhDispersed()),//(bulkFluid_.Dh()),(bulkFluid_.Dh()),
    A_(dict.get<scalar>("const")),
    B_(dict.get<scalar>("coeff")),
    C_(dict.get<scalar>("expRe")),
    D_(dict.get<scalar>("expPr")),
    usePeclet_(C_ == D_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::FFHeatTransferCoefficientModels::Nusselt::value
(
    const label& celli
) const
{
    //- I am creating a scalar on return to (hopefully) force Return Value
    //  Optimizations (RVOs, C++ performance stuff)
    scalar Dhi(max(Dh_[celli], 1e-4));
    if (B_ != 0)
    {
        if (usePeclet_)
            return
                scalar
                ( 
                    (kappa_[celli]/Dhi)*
                    (A_ + B_*pow(Re_[celli]*Pr_[celli], C_))
                );
        else
            return 
                scalar
                (
                    (kappa_[celli]/Dhi)*
                    (A_ + B_*pow(Re_[celli], C_)*pow(Pr_[celli], D_))
                );
    }
    else
        return scalar((kappa_[celli]/Dhi)*A_);

}


// ************************************************************************* //
