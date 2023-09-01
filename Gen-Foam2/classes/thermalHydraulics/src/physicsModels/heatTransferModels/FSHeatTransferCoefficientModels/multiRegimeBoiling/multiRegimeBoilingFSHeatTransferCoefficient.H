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

Class
    Foam::FSHeatTransferCoefficientModels::multiRegimeBoiling

Description
    Heat transfer coefficient that covers many regimes, from single phase
    convection, two-phase convection, sub-cooled and nucleate boiling, to
    film boiling (the latter, hopefully).

    Author: 
    Stefan Radman (sradman@pm.me / stefan.radman@epfl.ch)

SourceFiles
    multiRegimeBoiling.C

\*---------------------------------------------------------------------------*/

#ifndef multiRegimeBoilingFSHeatTransferCoefficient_H
#define multiRegimeBoilingFSHeatTransferCoefficient_H

#include "FSHeatTransferCoefficientModel.H"
#include "subCooledBoilingFractionModel.H"
#include "flowEnhancementFactorModel.H"
#include "suppressionFactorModel.H"
#include "TONBModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class FFPair;

namespace FSHeatTransferCoefficientModels
{

/*---------------------------------------------------------------------------*\
                           Class superposition Declaration
\*---------------------------------------------------------------------------*/

class multiRegimeBoiling
:
    public FSHeatTransferCoefficientModel
{
protected:

    //- Ref to wall temperature
    const volScalarField& Tw_;

    //- Ref to fluid temperature
    const volScalarField& Tf_;

    //- Ptr to FFPair (breaks encapsulation but so does the objectRegistry)
    mutable const FFPair* FFPairPtr_;

    //- Forced convection heat transfer coefficient
    autoPtr<FSHeatTransferCoefficientModel> htcFCPtr_;

    //- Pool boiling heat transfer coefficient
    autoPtr<FSHeatTransferCoefficientModel> htcPBPtr_;

    //- Film condensation heat transfer coefficient
    autoPtr<FSHeatTransferCoefficientModel> htcCndPtr_;

    //- Flow enhcement factor model: increases the convective heat trasfer
    //  coefficient to account for bubble turbulence, increase in flow velocity
    //  due to the vapour generation
    autoPtr<flowEnhancementFactorModel> FPtr_;

    //- Suppression factor model: "turns-off" the pool boiling component
    //  as the boiling transition to film boiling
    autoPtr<suppressionFactorModel> SPtr_;

    //- Fraction of sub-cooled heat transfer that results in vapour generation
    autoPtr<subCooledBoilingFractionModel> SCBFPtr_;

    //- Temperature of the onset of nucleate boiling model
    autoPtr<TONBModel> TONBPtr_;

    //- Ptr to mass transfer term at the wall (of the phaseChange class) so
    //  that it can be set here to account e.g. for sub-cooled boiling.
    //  Please, never forget that encapsulation leads to total madness if 
    //  followed religiously. Still I encourage people to find more
    //  elegant and simpler solutions, if any
    mutable volScalarField* dmdtWPtr_;

    //- Superposition exponent
    scalar exp_;
    scalar oneByExp_;

    //- Superposition is performed on heat fluxes rather than heat transfer 
    //  coefficients
    bool qMode_;

    //- Cache (yikes!)
    mutable scalar htc2pFCi_;
    mutable scalar TONBi_;
    mutable scalar htcNBi_;

//- Protected methods

    //- Set the subcooled boiling mass transfer term for celli 
    void setDmdtW
    (
        const label& celli, 
        const scalar& qNBi, 
        const scalar& qFCi
    ) const;

    //- Return the nucleate boiling heat flux in celli for a given value of the
    //  wall temperature (will be used to numerically calculate the derivative 
    //  of qNB, which is required to find the critical heat flux temperature 
    //  given the critical flux value)
    scalar qNB
    (
        const label& celli, 
        const scalar& Twi1
    ) const; 

public:

    TypeName("multiRegimeBoiling");

    //- Constructors

        //- Construct for a fluid-structure pair
        multiRegimeBoiling
        (
            const FSPair& pair,
            const dictionary& dict,
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~multiRegimeBoiling(){}


    //- Member Functions

        //- Return htc for cell celli
        virtual scalar value(const label& celli) const;

        
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FSHeatTransferCoefficientModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
