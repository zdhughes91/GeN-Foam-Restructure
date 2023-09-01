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

#if defined __has_include
#  if __has_include(<commDataLayer.H>) 
#    include <commDataLayer.H>
#    define isCommDataLayerIncluded
#  endif
#endif

#ifdef isCommDataLayerIncluded


#include "FSPair.H"
#include "NusseltWallAndHfromFMUFSHeatTransferCoefficient.H"
#include "addToRunTimeSelectionTable.H"
#include "commDataLayer.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace FSHeatTransferCoefficientModels
{
    defineTypeNameAndDebug(NusseltWallAndHfromFMU, 0);
    addToRunTimeSelectionTable
    (
        FSHeatTransferCoefficientModel,
        NusseltWallAndHfromFMU,
        FSHeatTransferCoefficientModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FSHeatTransferCoefficientModels::NusseltWallAndHfromFMU::NusseltWallAndHfromFMU
(
    const FSPair& pair,
    const dictionary& dict,
    const objectRegistry& objReg
)
:
    FSHeatTransferCoefficientModel
    (
        pair,
        dict,
        objReg
    ),
    Re_(pair.Re()),
    kappa_(pair.fluidRef().kappa()),
    Pr_(pair.fluidRef().Pr()),
    Dh_(pair.fluidRef().Dh()),
    A_(dict.get<scalar>("const")),
    B_(dict.get<scalar>("coeff")),
    C_(dict.get<scalar>("expRe")),
    D_(dict.get<scalar>("expPr")),
    H_wall(dict.get<scalar>("addH")),
    usePeclet_(C_ == D_),
    HNameFromFMU_("HCoupled")
{
    // Storing object in dataLayer
    word HKeyFromFMU("HNameFromFMU");
    if (dict.found(HKeyFromFMU))
    {
        HNameFromFMU_ = dict.get<word>(HKeyFromFMU);
        // Communicating with the FMU
        const Time& runTime = this->db().time();
        commDataLayer& data = commDataLayer::New(runTime); 
        // Store in data layer and set its initial value to 0 
        // in the dictionary      
        data.storeObj(
            scalar(0.0),
            HNameFromFMU_,
            commDataLayer::causality::in
            );
        Info << "Using FMUs for the Nusselt in " << dict.dictName() << endl;
    }

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::FSHeatTransferCoefficientModels::NusseltWallAndHfromFMU::value
(
    const label& celli
) const
{

    // Communicating with the FMU
    const Time& runTime = this->db().time();
    commDataLayer& data = commDataLayer::New(runTime);
    const scalar H_fromFMU =
        data.getObj<scalar>(HNameFromFMU_,commDataLayer::causality::in);
                
    //- I am creating a scalar on return to (hopefully) force Return Value
    //  Optimizations (RVOs, C++ performance stuff)
    scalar H_fluid;
    if (B_ != 0)
    {
        if (usePeclet_)
            H_fluid = scalar
            (
                (kappa_[celli]/Dh_[celli])*
                (A_ + B_*pow(Re_[celli]*Pr_[celli], C_))
            );
        else
            H_fluid = scalar
            (
                (kappa_[celli]/Dh_[celli])*
                (A_ + B_*pow(Re_[celli], C_)*pow(Pr_[celli], D_))
            );
    }
    else
        H_fluid = scalar((kappa_[celli]/Dh_[celli])*A_);

    //return(H_fluid*H_wall / (H_fluid + H_wall));
    return
        (1/
            (
                (1/max(SMALL,H_fluid))
                +(1/max(SMALL,H_wall))
                +(1/max(SMALL,H_fromFMU))
            )
        );
}

#endif
// ************************************************************************* //
