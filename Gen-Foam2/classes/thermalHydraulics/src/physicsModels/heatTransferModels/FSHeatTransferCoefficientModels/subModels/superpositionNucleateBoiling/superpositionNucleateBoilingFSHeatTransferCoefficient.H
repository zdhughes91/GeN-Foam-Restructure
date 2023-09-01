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
    Foam::FSHeatTransferCoefficientModels::superpositionNucleateBoiling

Description
    Heat transfer coefficient for both single-phase convective scenarios that
    evolve into nucleate boiling ones, but no sub-cooled boiling region in
    between.

    Author: 
    Stefan Radman (sradman@pm.me / stefan.radman@epfl.ch)

SourceFiles
    superpositionNucleateBoiling.C

\*---------------------------------------------------------------------------*/

#ifndef superpositionNucleateBoilingFSHeatTransferCoefficient_H
#define superpositionNucleateBoilingFSHeatTransferCoefficient_H

#include "FSHeatTransferCoefficientModel.H"
#include "flowEnhancementFactorModel.H"
#include "suppressionFactorModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace FSHeatTransferCoefficientModels
{

/*---------------------------------------------------------------------------*\
                           Class superposition Declaration
\*---------------------------------------------------------------------------*/

class superpositionNucleateBoiling
:
    public FSHeatTransferCoefficientModel
{
private:

    //- Forced convection heat transfer coefficient
    autoPtr<FSHeatTransferCoefficientModel> htcFCPtr_;

    //- Pool boiling heat transfer coefficient
    autoPtr<FSHeatTransferCoefficientModel> htcPBPtr_;

    //- 
    autoPtr<flowEnhancementFactorModel> FPtr_;

    //- 
    autoPtr<suppressionFactorModel> SPtr_;

public:

    TypeName("superpositionNucleateBoiling");

    //- Constructors

        //- Construct for a fluid-structure pair
        superpositionNucleateBoiling
        (
            const FSPair& pair,
            const dictionary& dict,
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~superpositionNucleateBoiling(){}


    //- Member Functions

        //- Drag coefficient
        virtual scalar value(const label& celli) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FSHeatTransferCoefficientModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //