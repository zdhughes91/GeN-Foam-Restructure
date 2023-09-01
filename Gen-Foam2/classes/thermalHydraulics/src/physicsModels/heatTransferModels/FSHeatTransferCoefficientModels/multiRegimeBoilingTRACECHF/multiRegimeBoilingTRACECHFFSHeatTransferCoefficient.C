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
#include "multiRegimeBoilingTRACECHFFSHeatTransferCoefficient.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace FSHeatTransferCoefficientModels
{
    defineTypeNameAndDebug(multiRegimeBoilingTRACECHF, 0);
    addToRunTimeSelectionTable
    (
        FSHeatTransferCoefficientModel, 
        multiRegimeBoilingTRACECHF, 
        FSHeatTransferCoefficientModels
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FSHeatTransferCoefficientModels::multiRegimeBoilingTRACECHF::
multiRegimeBoilingTRACECHF
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
    pair_(pair),
    p_(pair.mesh().lookupObject<volScalarField>("p")),
    Tw_(pair.structureRef().Twall()),
    Tf_(pair.fluidRef().thermo().T()),    
    // CHF Bool - 0 = preCHF, 1 = postCHF, 2 = transition region
    CHFBool_
    (
        IOobject
        (
            "CHFBoolMultiRegimeBoilingTRACECHF",
            pair.mesh().time().timeName(),
            pair.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        pair.mesh(),
        dimensionedScalar("",dimensionSet(0,0,0,0,0,0,0),0.0),
        zeroGradientFvPatchScalarField::typeName

    ),
    wfTB_
    (
        IOobject
        (
            "wfTransitionBoiling",
            pair.mesh().time().timeName(),
            pair.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        pair.mesh(),
        dimensionedScalar("",dimensionSet(0,0,0,0,0,0,0),0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    TLF_
    (
        IOobject
        (
            "LeidenfrostTemperature",
            pair.mesh().time().timeName(),
            pair.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        pair.mesh(),
        dimensionedScalar("",dimensionSet(0,0,0,1,0,0,0),0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
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
    TLFPtr_(nullptr),
    htcAFPtr_(nullptr),
    R_(dict.getOrDefault<scalar>("absoluteSurfaceRoughness", 4e-7)),
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
    if (this->found("leidenfrostModel"))
    {
        TLFPtr_.reset
        (
            TLFModel::New
            (
                pair,
                this->subDict("leidenfrostModel"),
                pair.mesh()
            )
        );
    }

    if (this->found("annularFlowModel"))
    {
        htcAFPtr_.reset
        (
            FSHeatTransferCoefficientModel::New
            (
                pair,
                this->subDict("annularFlowModel"),
                pair.mesh()
            )
        );
    }
    if (this->found("criticalHeatFluxModel"))
    {
        qCHFPtr_.reset
        (
            CHFModel::New
            (
                pair,
                this->subDict("criticalHeatFluxModel"),
                pair.mesh()
            )
        );
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar 
Foam::FSHeatTransferCoefficientModels::multiRegimeBoilingTRACECHF::value
(
    const label& celli
) const
{
    // - Pointers breaking encapsulation - needed
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
    // ---------------------------------- //
    //--- Refs --//
    const scalar& alphai(1.0-this->pair_.fluidRef().normalized()[celli]);
    const scalar& Twi(Tw_[celli]);
    const scalar& Tfi(Tf_[celli]);
    const scalar& Tsati(FFPairPtr_->iT()[celli]);
    const scalar& pi(pair_.mesh().lookupObject<volScalarField>("p")[celli]);
    
    //CHF Bool : 0 if preCHF, 1 if postCHF, 2 if transition
    CHFBool_[celli]=0;
    // weight for the transition boiling -> see TRACE
    wfTB_[celli] = 0;
    // Leidenfrost Temperature 
    TLF_[celli] = 0;

    // Heated area per volume 
    const scalar Areai(pair_.structureRef().iAact()[celli]);
    // Latent heat value 
    const scalar Li(mag(FFPairPtr_->L()[celli]));

    //- Two-phase forced-convection heat transfer coefficient
    //  and explicit heat flux, used later throughout the map
    scalar htc2pFCi(htcFCPtr_->value(celli)*FPtr_->value(celli));
    scalar qFCi(htc2pFCi*(Twi-Tfi));

    if (Twi < Tsati)    //- Wall below saturation, either condensing or
                        //  nothing special happens
    {
        if (htcCndPtr_.valid()) //- If a film condensation model has been specified
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
        //Info << "TONB --- " << TONBi << " ---" << endl;
        if (Twi < TONBi)    //- Wall above saturation but below onset of
                            //  nucleate boiling, nothing special
        {
            return htc2pFCi;
        }
        else //- Wall above onset of nucleate boiling
// -------------------------- Pre CHF Situation -------------------------- // 
        {
            scalar qCHFi // Critical Heat Flux -> Only constant coded yet 
            (
                (qCHFPtr_.valid()) ?
                qCHFPtr_->value(celli) :
                1e15
            );

// ------------------------------------------------------------------------ //
            // DETERMINATION OF TCHF //
// ------------------------------------------------------------------------ //
            //Computation of TCHF accurately by using only the Pool Boiling Model 
            // Will not be always valid if the model for PB is changed 

            // Gorenflo Model used             
            const scalar pCrit(2.209e7); //- Specific to Water
            const scalar h0(5600);
            const scalar q0(20000);
            const scalar R0(4e-7);
            scalar A(pow(R_/R0, 0.133));

            scalar pR(pi/pCrit);
            pR = min(pR,0.99999);
            scalar F(1.73*pow(pR, 0.27) + (6.1+0.68/(1.0-pR))*sqr(pR));
            scalar n(0.9-0.3*pow(pR, 0.15));

            scalar TCHF(Tsati+pow(qCHFi,1.0-n)*pow(q0,n)/h0/F/A);

            // Print of Temperature for checking results
            //Info << "Twall --- " << Twi << " ---" << endl;
            //Info << "TCHF  --- " << TCHF << " ---" << endl;
            //Info << "TONB --- " << TONBi << " ---" << endl;
            //Info << "Tsat --- " << Tsati << " ---" << endl; 

// ------------------------------------------------------------------------ // 
// ------------------------------------------------------------------------ // 

            scalar dT(Twi-Tfi);
            dT = (dT >= 0.0) ? max(dT, 1e-3) : min(dT, -1e-3);
            if (Twi<TCHF) // This part was already coded by Stefan
            {
                scalar hPBi(htcPBPtr_->value(celli));
                scalar qPBi(hPBi*(Twi-Tsati));

                //- Calc hPB at boiling onset, i.e. the hPB when the wall 
                //  temperature is TONBi. This is done by caching Twi, 
                //  setting it to TONB,  using the htcPB model to calc the 
                //  hPB with Tw=TONB, then re-setting the Twi
                scalar Twi0(Twi);
                const_cast<scalar&>(Twi) = TONBi;
                scalar hPBONBi(htcPBPtr_->value(celli));
                const_cast<scalar&>(Twi) = Twi0;


                //- Explicit heat fluxes
                scalar qBIi(hPBONBi*(TONBi-Tsati));
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
                    Areai*qSCDmdti/Li;

                // Print of Heat Flux for checking results
                //Info << "qNBi --- " << qNBi << " ---" << endl;
                //Info << "qFCi --- " << qFCi << " ---" << endl;
                //Info << "qPBi --- " << qPBi << " ---" << endl;                
                //Info << "qBIi --- " << qBIi << " ---" << endl;
                return qNBi/dT;
            }

// -------------------------- Post CHF Situation -------------------------- // 
            else    
            {   
                CHFBool_[celli] = 1;
                //- If a TLF model has not been provided, the Leidenfrost temperature is 
                // set to Tw -> no hysteresis zone. 
                scalar TLFi
                (
                    (TLFPtr_.valid()) ?
                    TLFPtr_->value(celli) :
                    Twi
                );
                TLF_[celli]=TLFi; 
                
                // For now, only annular flow regime for FILM BOILING (FB) -> alpha < 0.6 needed

                if (Twi<TLFi) // Transition boiling region - hysteresis zone
                {
                    CHFBool_[celli] = 2; // Transition

                    // Film Boiling HTC needed at the TLF temperature for Film Boiling (see TRACE), here Inverted 
                    // Annular flow is considered
                    scalar Twi0(Twi);
                    const_cast<scalar&>(Twi) = TLFi;
                    scalar hFBiMIN
                    (
                        (htcAFPtr_.valid()) ?
                        htcAFPtr_->value(celli) :
                        0
                    );
                    // Production of vapour at TLF 
                    scalar hGammaiMIN(pair_.mesh().lookupObject<volScalarField>("hGammaCachard")[celli]);
                    const_cast<scalar&>(Twi) = Twi0;


                    scalar dTCHF(TCHF-TLFi);
                    dTCHF = (dTCHF >= 0.0) ? max(dTCHF, 1e-6) : min(dTCHF, -1e-6);
                    //scalar wfTB(sqrt(1-alphai)*sqr((Twi-TLFi)/dTCHF));
                    // Weight for the transition boiling - also needed for MultiRegimeVapour
                    scalar wfTB(sqr((Twi-TLFi)/dTCHF));
                    wfTB_[celli] = wfTB;

                    dT=Twi-Tfi;
                    dT = (dT >= 0.0) ? max(dT, 1e-6) : min(dT, -1e-6);
                    scalar hwlTB(qCHFi*wfTB/dT);

                    // Same as Stefan did :
                    (*dmdtWPtr_)[celli] = 
                        (qCHFi*Areai*wfTB+Areai*(1-wfTB)*hGammaiMIN*(Twi-Tfi))/max(Li,1e-6);
                    return (hwlTB+(1-wfTB)*hFBiMIN);
            
                }
                else // Transition is ended -> Film Boiling HTC - Here only Inverted Annular Flow 
                { 
                    scalar hFBi
                    (
                        (htcAFPtr_.valid()) ?
                        htcAFPtr_->value(celli) :
                        0
                    );
                    // net production of Vapour in FB IAF
                    const volScalarField& hGamma_(pair_.mesh().lookupObject<volScalarField>("hGammaCachard"));
                    scalar hGammai(hGamma_[celli]);
        
                    (*dmdtWPtr_)[celli] = 
                        Areai*hGammai*(Twi-Tfi)/max(Li,1e-6);
                    return hFBi;
                }
            }
        } 
    }
}

// ************************************************************************* //
