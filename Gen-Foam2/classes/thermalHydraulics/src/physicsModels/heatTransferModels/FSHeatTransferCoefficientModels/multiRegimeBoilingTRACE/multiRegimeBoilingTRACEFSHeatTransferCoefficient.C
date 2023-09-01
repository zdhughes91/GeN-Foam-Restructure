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
#include "FFPair.H"
#include "multiRegimeBoilingTRACEFSHeatTransferCoefficient.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace FSHeatTransferCoefficientModels
{
    defineTypeNameAndDebug(multiRegimeBoilingTRACE, 0);
    addToRunTimeSelectionTable
    (
        FSHeatTransferCoefficientModel, 
        multiRegimeBoilingTRACE, 
        FSHeatTransferCoefficientModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FSHeatTransferCoefficientModels::multiRegimeBoilingTRACE::
multiRegimeBoilingTRACE
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
    Tw_(pair.structureRef().Twall()),
    Tf_(pair.fluidRef().thermo().T()),
    FFPairPtr_(nullptr),
    htcFCPtr_
    (
        FSHeatTransferCoefficientModel::New
        (
            pair,
            this->subDict("forcedConvectionModel"),
            pair.mesh()
        )
    ),
    htcPBPtr_
    (
        FSHeatTransferCoefficientModel::New
        (
            pair,
            this->subDict("poolBoilingModel"),
            pair.mesh()
        )
    ),
    FPtr_
    (
        flowEnhancementFactorModel::New
        (
            pair,
            this->subDict("flowEnhancementFactorModel"),
            pair.mesh()
        )
    ),
    SPtr_
    (
        suppressionFactorModel::New
        (
            pair,
            this->subDict("suppressionFactorModel"),
            pair.mesh()
        )
    ),
    SCBFPtr_(nullptr),
    TONBPtr_(nullptr),
    dmdtWPtr_(nullptr)
{
    if (this->found("filmCondensationModel"))
    {
        htcCndPtr_.reset
        (
            FSHeatTransferCoefficientModel::New
            (
                pair,
                this->subDict("filmCondensationModel"),
                pair.mesh()
            )
        );
    }
    if (this->found("nucleateBoilingOnsetModel"))
    {
        TONBPtr_.reset
        (
            TONBModel::New
            (
                pair,
                this->subDict("nucleateBoilingOnsetModel"),
                pair.mesh()
            )
        );
    }
    if (this->found("subCooledBoilingFractionModel"))
    {
        SCBFPtr_.reset
        (
            subCooledBoilingFractionModel::New
            (
                pair,
                this->subDict("subCooledBoilingFractionModel"),
                pair.mesh()
            )
        );
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar 
Foam::FSHeatTransferCoefficientModels::multiRegimeBoilingTRACE::value
(
    const label& celli
) const
{
    if (FFPairPtr_ == nullptr)
    {
        HashTable<const FFPair*> FFPairs(pair_.mesh().lookupClass<FFPair>());
        FFPairPtr_ = FFPairs[FFPairs.toc()[0]];
    }
    if (dmdtWPtr_ == nullptr)
    {
        //- Ptr to wall mass transfer term (i.e. due to subcooled boiling)
        dmdtWPtr_ = 
            &
            (
                pair_.mesh().lookupObjectRef<volScalarField>
                (
                    "dmdtW."+FFPairPtr_->name()
                )
            );
    }

    //- Refs
    const scalar& alphai(1.0-this->pair_.fluidRef().normalized()[celli]);
    const scalar& Twi(Tw_[celli]);
    const scalar& Tfi(Tf_[celli]);
    const scalar& Tsati(FFPairPtr_->iT()[celli]);

    //- Two-phase forced-convection heat transfer coefficient
    //  and explicit heat flux, used later throughout the map
    scalar htc2pFCi(htcFCPtr_->value(celli)*FPtr_->value(celli));
    scalar qFCi(htc2pFCi*(Twi-Tfi));

    /*-----------------------------------------------------------------------*/
    /*- HTC MAP TAKEN FROM TRACE, BUT THIS SUPPORTS USER-SELECTABLE MODELS --*/
    /*-----------------------------------------------------------------------*/
    if (Twi < Tsati)    //- Wall below saturation, either condensing or
                        //  nothing special happens
    {
        if (htcCndPtr_.valid()) //- If a film condensation model has been
                                //  specified
        {
            if (alphai <= 0.8) //- Nothing special
                return htc2pFCi;
            else    //- Film condensation
            {
                //- Film condensation heat transfer coefficient
                scalar htcCndi(htcCndPtr_->value(celli));
                
                //- Linear interpolation if 0.8 < alphai < 0.9
                if (alphai < 0.9)
                {
                    scalar f(0.9-alphai/0.1);
                    return f*htc2pFCi+(1.0-f)*htcCndi;
                }
                else //- Film condensation only
                    return htcCndi;
            }   
        }
        else
            return htc2pFCi;
    }
    else
    {
        //- If a TONB model has not been provided, no sub-cooled boiling
        //  and set it equal to Tsati
        scalar TONBi
        (
            (TONBPtr_.valid()) ?
            TONBPtr_->value(celli, htc2pFCi) :
            Tsati
        );
        if (Twi < TONBi)    //- Wall above saturation but below onset of
                            //  nucleate boiling, nothing special
        {
            return htc2pFCi;
        }
        else //- Wall above onset of nucleate boiling
        {
            scalar TCHFi(1e69); //- Needs dedicated model for its setting
            if (Twi < TCHFi)    //- Below CHF, i.e. either nucleate boiling
                                //  or subcooled boiling
            {
                //- Calc hPB at boiling onset, i.e. the hPB when the wall 
                //  temperature is TONBi. This is done by caching Twi, 
                //  setting it to TONB,  using the htcPB model to calc the 
                //  hPB with Tw=TONB, then re-setting the Twi
                scalar Twi0(Twi);
                const_cast<scalar&>(Twi) = TONBi;
                scalar hPBONBi(htcPBPtr_->value(celli));
                const_cast<scalar&>(Twi) = Twi0;
                scalar hPBi(htcPBPtr_->value(celli));

                //- Explicit heat fluxes
                scalar qBIi(hPBONBi*(TONBi-Tsati));
                scalar qPBi(hPBi*(Twi-Tsati));
                scalar qNBi = pow(pow(qFCi, 3)+pow(qPBi-qBIi, 3), 1.0/3.0);

                //- Please note that only the dmdtW is set in this scope as
                //  the subcooled boiling heat transfer coefficient and the 
                //  nucleate boiling heat transfer coefficient are computed 
                //  in the same way, outside and after this scope
                scalar f(1.0);
                if (SCBFPtr_.valid())
                {
                    f = SCBFPtr_->value(celli, qNBi);
                }
                scalar qSCDmdti(f*(qNBi-qFCi));

                //- If fluid1 is liquid and 2 is vapour then L > 0 and 
                //  this sub-cooled boiling term is also > 0. If fluid2
                //  is liquid and fluid1 is vapour then L < 0 and 
                //  everything still  works out in terms of the dmdtW 
                //  sign, as I recall that it is positive for phase
                //  changes from fluid1 to fluid2 and negative 
                //  vice-versa
                (*dmdtWPtr_)[celli] = 
                    pair_.structureRef().iAact()[celli]*qSCDmdti/
                    mag(FFPairPtr_->L()[celli]);

                scalar dT(Twi-Tfi);
                dT = (dT >= 0.0) ? max(dT, 1e-3) : min(-dT, -1e-3);
                return qNBi/dT;
            }
            else    //- Post-CHF heat transfer, currently missing 
                    //  implementation
            {
                return 1.0;
            }
        }
    }
}

// ************************************************************************* //
