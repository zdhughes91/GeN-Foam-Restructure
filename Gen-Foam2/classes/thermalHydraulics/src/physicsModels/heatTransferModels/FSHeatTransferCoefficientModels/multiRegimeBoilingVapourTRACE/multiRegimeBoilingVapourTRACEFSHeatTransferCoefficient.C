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
#include "multiRegimeBoilingVapourTRACEFSHeatTransferCoefficient.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace FSHeatTransferCoefficientModels
{
    defineTypeNameAndDebug(multiRegimeBoilingVapourTRACE, 0);
    addToRunTimeSelectionTable
    (
        FSHeatTransferCoefficientModel, 
        multiRegimeBoilingVapourTRACE, 
        FSHeatTransferCoefficientModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FSHeatTransferCoefficientModels::multiRegimeBoilingVapourTRACE::
multiRegimeBoilingVapourTRACE
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
    p_(pair.mesh().lookupObject<volScalarField>("p")),
    Tw_(pair.structureRef().Twall()),
    CHFBool_
    (
                pair_.mesh().lookupObject<volScalarField>
                (
                    "CHFBoolMultiRegimeBoilingTRACECHF"
                )
    ),
    wfTB_
    (
                pair_.mesh().lookupObject<volScalarField>
                (
                    "wfTransitionBoiling"
                )
    ),
    TLF_
    (
                pair_.mesh().lookupObject<volScalarField>
                (
                    "LeidenfrostTemperature"
                )
    ),
    FFPairPtr_(nullptr),
    htcAFPtr_(nullptr)
{
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
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar 
Foam::FSHeatTransferCoefficientModels::multiRegimeBoilingVapourTRACE::value
(
    const label& celli
) const
{
    if (FFPairPtr_ == nullptr)
    {
        HashTable<const FFPair*> FFPairs(pair_.mesh().lookupClass<FFPair>());
        FFPairPtr_ = FFPairs[FFPairs.toc()[0]];
    }

    //- Refs
    const scalar& Twi(Tw_[celli]);
    const scalar& POSTCHFi(CHFBool_[celli]);    // 0 = PreCHF, 1 = Transition, 2 = PostCHF
    const scalar& wfTBi(wfTB_[celli]);          // Weight for transition
    const scalar& TLFi(TLF_[celli]);            // Leidenfrost Temperature -> Used for transition

    if (POSTCHFi>=1)    //- In POSTCHF MODEL
    {
        // If alpha < 0.6, inverted annular model used, else need to be coded = later.
        /*if (alphai<0.6)
        {
            // Annular flow model 
        }
        else 
        {
            // Need to be coded - Invert slug + Dispersed flow 
            return 0.0  
        }*/
        
        if (POSTCHFi==1) // Post CHF Region only for inverted annular
        {   
            scalar hFBi(htcAFPtr_->value(celli));
            return hFBi; 
        }
        else // Transition region 
        { 
            // HTC at TLF for Film Boiling (see TRACE) - here only inverted annular 
            scalar Twi0(Twi);
            const_cast<scalar&>(Twi) = TLFi;
            scalar hFBiMIN(htcAFPtr_->value(celli));
            const_cast<scalar&>(Twi) = Twi0;
            return (1-wfTBi)*hFBiMIN;
        }
    }
    else // Pre-CHF -> no heat flux from wall to vapour - return 0
    {
        return 0.0;
    }
}
// ************************************************************************* //
