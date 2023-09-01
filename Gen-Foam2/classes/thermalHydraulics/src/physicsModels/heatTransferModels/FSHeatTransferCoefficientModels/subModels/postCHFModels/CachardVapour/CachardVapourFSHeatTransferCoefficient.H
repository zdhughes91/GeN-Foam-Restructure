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
    Foam::FSHeatTransferCoefficientModels::CachardVapour

Description
    Heat transfer coefficient for Inverted Annular Flow - for the vapour phase.
    This code cooperate with the CachardLiquidFSHeatTransferCoefficient.  
    Warning, this model has been coded for rod bundle, especially for the delta model,
    the thickness of the film of vapor. 

    Author: 
    Gauthier LAZARE - gauthier.lazare@gmail.com

SourceFiles
    CachardVapourFSHeatTransferCoefficient.C

\*---------------------------------------------------------------------------*/

#ifndef CachardVapourFSHeatTransferCoefficient_H
#define CachardVapourFSHeatTransferCoefficient_H

#include "FSHeatTransferCoefficientModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace FSHeatTransferCoefficientModels
{

/*---------------------------------------------------------------------------*\
                           Class CachardLiquid Declaration
\*---------------------------------------------------------------------------*/

class CachardVapour
:
    public FSHeatTransferCoefficientModel
{
private:

    //- Ref to pressure
    const volScalarField& p_;

    //- Ref to normalized void fraction
    const volScalarField& alpha_;

    //- Ref to the thermal conductivity of the vapour 
    const volScalarField& kg_;

    //- Ref to the hydraulic diameter
    const volScalarField& Dh_;

    //- Reference parameter
    //- Gravity constant
    const scalar gravity_;

    //- Values send to CachardLiquid
    /*//- hInit = k_g/delta -> initial heat transfer coefficient without turbulence
    mutable volScalarField hInit_;

    //- Nusselt Number 
    mutable volScalarField Nuwl_;*/
    

public:

    TypeName("CachardVapour");

    //- Constructors

        //- Construct for a fluid-structure pair
        CachardVapour
        (
            const FSPair& pair,
            const dictionary& dict,
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~CachardVapour(){}


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