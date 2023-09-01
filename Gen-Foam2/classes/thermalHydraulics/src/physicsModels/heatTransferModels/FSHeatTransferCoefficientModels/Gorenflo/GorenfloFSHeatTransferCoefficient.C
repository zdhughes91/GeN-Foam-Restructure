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
#include "GorenfloFSHeatTransferCoefficient.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace FSHeatTransferCoefficientModels
{
    defineTypeNameAndDebug(Gorenflo, 0);
    addToRunTimeSelectionTable
    (
        FSHeatTransferCoefficientModel, 
        Gorenflo, 
        FSHeatTransferCoefficientModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FSHeatTransferCoefficientModels::Gorenflo::Gorenflo
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
    useExplicitHeatFlux_
    (
        this->lookupOrDefault<bool>("useExplicitHeatFlux", false)
    ),
    Twall_
    (
        (!useExplicitHeatFlux_) ?
        &pair.structureRef().Twall() :
        nullptr
    ),
    Tsat_
    (
        (!useExplicitHeatFlux_) ?
        &pair.mesh().lookupObject<volScalarField>("T.interface") :
        nullptr
    ),
    q_
    (
        (useExplicitHeatFlux_) ? 
        &pair.mesh().lookupObject<volScalarField>("heatFlux.structure") :
        nullptr
    ),
    p_(pair.mesh().lookupObject<volScalarField>("p")),
    pCrit_(2.209e7), //- Specific to Water
    h0_(5600),
    q0_(20000),
    R0_(4e-7),
    R_(dict.getOrDefault<scalar>("absoluteSurfaceRoughness", R0_)),
    A_(pow(R_/R0_, 0.133))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

/*
IMPLEMENTATION NOTES

The original pool boiling heat transfer coefficient by Gorenflo is the following 
form:

htcPB = h0_*F*((q/q0_)^n)*(R_/R0_)^0.133

with F being a function of the reduced pressure, i.e. p/pCtri_. Needless to 
say, this form is not particularly suitable for a numerical implementation in 
its current form as htcPB is required to compute q in the first place. In 
general, the pool boiling heat flux is given by:

q = htc*(Twall-Tsat)

with htc being the total heat transfer coefficient between wall and fluid,
Twall being the wall temperature and Tsat the saturation temperature. By 
assuming that htc ~ htcPB (i.e. that most of the heat flux is due to the pool
boiling mechanism), one can tha equate:

htcPB = h0_*F*((htcPB*(Twall-Tsat)/q0_)^n)*(R_/R0_)^0.133

which can be re-arranged to isolate htcPB on the left hand side as:

htcPB = (h0_*F*((Twall-Tsat)^n)*(R_/R0_)^0.133/(q0_^n))^(1.0/(1.0-n))

Which is how the correlation is implemented here.
The assumption of htc ~ htcPB can be rather crude but it is widely employed
to yield a feasable numerical implementation in other computer codes too (e.g.
TRACE)
*/

Foam::scalar Foam::FSHeatTransferCoefficientModels::Gorenflo::value
(
    const label& celli
) const
{
    if (!useExplicitHeatFlux_)
    {
        scalar deltaT((*Twall_)[celli]-(*Tsat_)[celli]);
        if (deltaT > 0.0)
        {
            scalar pR(p_[celli]/pCrit_);
            scalar F(1.73*pow(pR, 0.27) + (6.1+0.68/(1.0-pR))*sqr(pR));
            scalar n(0.9-0.3*pow(pR, 0.15));
            return pow(h0_*A_*F*pow(deltaT/q0_,n), (1.0/(1.0-n)));
        }
        else
            return 0.0;
    }
    else
    {
        const scalar& q((*q_)[celli]);
        if (q > 0.0)
        {
            //- Reduced pressure
            scalar pR(p_[celli]/pCrit_);
            scalar F(1.73*pow(pR, 0.27) + (6.1+0.68/(1.0-pR))*sqr(pR));
            scalar n(0.9-0.3*pow(pR, 0.15));
            return h0_*F*pow(q/q0_, n)*A_;
        }
        else
            return 0.0;
    }
}

// ************************************************************************* //
