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
    Foam::FSHeatTransferCoefficientModels::CachardLiquid

Description
    Heat transfer coefficient for Inverted Annular Flow - for the liquid phase. 
    This code cooperate with the CachardVapourFSHeatTransferCoefficient. 
    Warning, this model has been coded for rod bundle, especially for delta model,
    the thickness of the film of vapor. 

    Author:
    Gauthier LAZARE - gauthier.lazare@gmail.com

SourceFiles
    CachardLiquidFSHeatTransferCoefficient.C

\*---------------------------------------------------------------------------*/

#ifndef CachardLiquidFSHeatTransferCoefficient_H
#define CachardLiquidFSHeatTransferCoefficient_H

#include "FSHeatTransferCoefficientModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace FSHeatTransferCoefficientModels
{

/*---------------------------------------------------------------------------*\
                           Class CachardLiquid Declaration
\*---------------------------------------------------------------------------*/

class CachardLiquid
:
    public FSHeatTransferCoefficientModel
{
private:

    //- Ref to structure wall temperature
    const volScalarField& Twall_;
 
    //- Ref to saturation temperature 
    const volScalarField& Tsat_;

    //- Ref to liquid temperature 
    const volScalarField& Tfluid_;

    //- Ref to normalized void fraction
    const volScalarField& alpha_;

    //- Ref to pressure
    const volScalarField& p_;

    //- Ref to normalized void fraction
    const volScalarField& rhoLiquid_;

    //-Ref to DR rod diameter
    const volScalarField&  rhoVapour_;

    //-Ref to kg kappa vapour
    const volScalarField&  kg_;

    //-Ref to mug mu vapour
    const volScalarField&  mug_;

    //-Ref to Dh vapour
    const volScalarField&  Dh_;
    
    //- Ref to the mass transfer rate - breaking encapsulation to reach MultiRegimeBoiling
    mutable volScalarField hGamma_;

    //- Reference parameter
    //- Stefan Boltzman Constant in W/m2/K
    scalar sigmaSB_;
    //- Gravity 
    scalar gravity_;

    //- Input parameters
    //- Wall emissivity 
    scalar epsilonWall_;
    //- Liquid emissivity
    scalar epsilonLiq_;

    

public:

    TypeName("CachardLiquid");

    //- Constructors

        //- Construct for a fluid-structure pair
        CachardLiquid
        (
            const FSPair& pair,
            const dictionary& dict,
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~CachardLiquid(){}


    //- Member Functions

        //- Heat transfer coefficient
        virtual scalar value(const label& celli) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FSHeatTransferCoefficientModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //