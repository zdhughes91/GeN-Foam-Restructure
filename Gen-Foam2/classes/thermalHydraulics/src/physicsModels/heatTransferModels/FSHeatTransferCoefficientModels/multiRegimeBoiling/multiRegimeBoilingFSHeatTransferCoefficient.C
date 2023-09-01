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
#include "multiRegimeBoilingFSHeatTransferCoefficient.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace FSHeatTransferCoefficientModels
{
    defineTypeNameAndDebug(multiRegimeBoiling, 0);
    addToRunTimeSelectionTable
    (
        FSHeatTransferCoefficientModel, 
        multiRegimeBoiling, 
        FSHeatTransferCoefficientModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FSHeatTransferCoefficientModels::multiRegimeBoiling::
multiRegimeBoiling
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
    dmdtWPtr_(nullptr),
    exp_(this->get<scalar>("superpositionExponent")),
    oneByExp_(1.0/exp_),
    qMode_(this->get<bool>("heatFluxSuperposition")),
    htc2pFCi_(0.0),
    TONBi_(0.0),
    htcNBi_(0.0)
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
    if (this->found("suppressionFactorModel"))
    {
        SPtr_.reset
        (
            suppressionFactorModel::New
            (
                pair,
                this->subDict("suppressionFactorModel"),
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

// * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * * //

void Foam::FSHeatTransferCoefficientModels::multiRegimeBoiling::setDmdtW
(
    const label& celli,
    const scalar& qNBi,
    const scalar& qFCi
) const
{
    //- No sub-cooled boiling possible if a model for the temperature of 
    //  onest of nucleate boiling (ONB) is not provided
    if (TONBPtr_.valid())
    {
        //- Fraction of the "real" (i.e. after superposition) pool boiling
        //  heat flux that results in net vapour generation
        scalar f(1.0);
        if (SCBFPtr_.valid())
        {
            f = SCBFPtr_->value(celli, qNBi);
        }

        //- Heat flux that results in vapour generation
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
            FFPairPtr_->L()[celli];
    }
}

Foam::scalar Foam::FSHeatTransferCoefficientModels::multiRegimeBoiling::qNB
(
    const label& celli,
    const scalar& Twi1
) const
{
    const scalar& Tfi(Tf_[celli]);
    const scalar& Tsati(FFPairPtr_->iT()[celli]);
    
    //- Cache and overwrite current wall temperature in celli if requested
    const scalar& Twi(Tw_[celli]);
    scalar Twi0(Twi);
    bool cacheValues(Twi1 == Twi);
    if (!cacheValues)
        const_cast<scalar&>(Twi) = Twi1;

    //- Enhanced two-phase forced convection heat flux, pool boiling htc and
    //  its heat flux
    scalar qFCi(htc2pFCi_*(Twi-Tfi));
    scalar htcPBi(htcPBPtr_->value(celli));
    scalar qPBi(htcPBi*(Twi-Tsati));

    //- Overall nucleate boiling htc and heat flux;
    scalar qNBi(0.0);

    //- If computing the nucleate boiling heat flux via the TRACE approach 
    //  (i.e. superposing the heat fluxes, not the htcs)
    if (qMode_)
    {
        if (SPtr_.valid())  //- If a suppression factor model is the one that
                            //  is used to turn-off pool boiling as the flow
                            //  becomes increasingly annular
            qNBi = 
                pow
                (
                    pow(qFCi, exp_) + pow(SPtr_->value(celli)*qPBi, exp_), 
                    oneByExp_
                );
        else if (TONBPtr_.valid())  //- If doing things the TRACE way, i.e.
                                    //  avoiding discontinuities at the onest
                                    //  of boiling via a qBI term
        {
            //- Calc pool boiling heat flux at the onset of boiling (ONB, 
            //  needless to say, cannot do that if a TONB model has not been
            //  specified)
            scalar Twi00(Twi);
            const_cast<scalar&>(Twi) = TONBi_;
            scalar hPBONBi(htcPBPtr_->value(celli));
            const_cast<scalar&>(Twi) = Twi00;
            scalar qBIi(hPBONBi*(TONBi_-Tsati));
            qNBi = pow(pow(qFCi, exp_) + pow(qPBi-qBIi, exp_), oneByExp_);
        }
        else //- If previous methods are not applicable
            qNBi = pow(pow(qFCi, exp_) + pow(qPBi, exp_), oneByExp_);
        if (cacheValues)
        {
            scalar dT(Twi-Tfi);
            dT = (dT >= 0.0) ? max(dT, 1e-3) : min(-dT, -1e-3);
            htcNBi_ = qNBi/dT;
        }
    }
    //- If computing the nucleate boiling heat flux via other approaches that
    //  superpose the htcs directly
    else
    {
        scalar htcNBi(0.0);
        if (SPtr_.valid())
            htcNBi =
            (
                pow
                (
                    pow(htc2pFCi_, exp_)+pow(SPtr_->value(celli)*htcPBi, exp_), 
                    oneByExp_
                )
            );
        else
            htcNBi =
            (
                pow(pow(htc2pFCi_, exp_) + pow(htcPBi, exp_), oneByExp_)
            );
        qNBi = htcNBi_*(Twi-Tfi);
        if (cacheValues)
            htcNBi_ = htcNBi;
    }

    //- Calculate subcooled boiling term
    if (cacheValues) // and Tfi < Tsati)
        setDmdtW(celli, qNBi, qFCi);
    else //- Reset wall temperature in celli and return qNB
        const_cast<scalar&>(Twi) = Twi0;
        
    return qNBi;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar 
Foam::FSHeatTransferCoefficientModels::multiRegimeBoiling::value
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
    const scalar& Tsati(FFPairPtr_->iT()[celli]);

    //- Two-phase forced-convection heat transfer coefficient
    //  and explicit heat flux, used later throughout the map
    htc2pFCi_ = (htcFCPtr_->value(celli)*FPtr_->value(celli));

    /*-----------------------------------------------------------------------*/
    /* HTC MAP FOLLOWS TRACE PHILOSOPHY, BUT SUPPORTS USER-SELECTABLE MODELS */
    /*-----------------------------------------------------------------------*/
    if (Twi < Tsati)    //- Wall below saturation, either condensing or
                        //  nothing special happens
    {
        if (htcCndPtr_.valid()) //- If a film condensation model has been
                                //  specified
        {
            if (alphai <= 0.8) //- Nothing special
                return htc2pFCi_;
            else    //- Film condensation
            {
                //- Film condensation heat transfer coefficient
                scalar htcCndi(htcCndPtr_->value(celli));
                
                //- Linear interpolation if 0.8 < alphai < 0.9
                if (alphai < 0.9)
                {
                    scalar f(0.9-alphai/0.1);
                    return f*htc2pFCi_+(1.0-f)*htcCndi;
                }
                else //- Film condensation only
                    return htcCndi;
            }   
        }
        else
            return htc2pFCi_;
    }
    else
    {
        //- If a TONB model has not been provided, no sub-cooled boiling
        //  and set it equal to Tsati
        TONBi_ =
        (
            (TONBPtr_.valid()) ?
            TONBPtr_->value(celli, htc2pFCi_) :
            Tsati
        );
        if (Twi < TONBi_)    //- Wall above saturation but below onset of
                            //  nucleate boiling, nothing special
        {
            return htc2pFCi_;
        }
        else //- Wall above saturation and onset of nucleate boiling
        {
            scalar TCHFi(1e69); //- Needs dedicated model for its setting
            if (Twi < TCHFi)    //- Below CHF, i.e. either nucleate boiling
                                //  or subcooled boiling
            {
                //- Calling qNBi also calculates the htcNBi_ so that is why
                //  it is called here
                scalar qNBi(qNB(celli, Twi));
                (void) qNBi; //- Suppress unused variable warining
                return htcNBi_;
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
