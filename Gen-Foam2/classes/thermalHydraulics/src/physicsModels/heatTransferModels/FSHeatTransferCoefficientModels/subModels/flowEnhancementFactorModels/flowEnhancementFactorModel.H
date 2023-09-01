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
    Foam::flowEnhancementFactorModel

Description

    Author: 
    Stefan Radman (sradman@protonmail.com; stefan.radman@epfl.ch)

SourceFiles
    flowEnhancementFactorModel.C
    newFlowEnhancementFactorModel.C

\*---------------------------------------------------------------------------*/

#ifndef flowEnhancementFactorModel_H
#define flowEnhancementFactorModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fluid.H"
#include "structure.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

//- Put them in a new FSHeatTransferCoefficientSubModels namespace? Eh, I don't 
//  know how handy this is

//-
class FSPair;

/*---------------------------------------------------------------------------*\
                           Class flowEnhancementFactorModel Declaration
\*---------------------------------------------------------------------------*/

class flowEnhancementFactorModel
:
    public IOdictionary
{
protected:

    //- Ref to pair
    const FSPair& pair_;

    //- Max allowable value
    const scalar maxValue_;

public:

    TypeName("flowEnhancementFactorModel");

    //- Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            flowEnhancementFactorModel,
            flowEnhancementFactorModels,
            (
                const FSPair& pair,
                const dictionary& dict,
                const objectRegistry& objReg
            ),
            (pair, dict, objReg)
        );

    //- Constructors

        flowEnhancementFactorModel
        (
            const FSPair& pair,
            const dictionary& dict,
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~flowEnhancementFactorModel(){}

    //- Selectors

        static autoPtr<flowEnhancementFactorModel> New
        (
            const FSPair& pair,
            const dictionary& dict,
            const objectRegistry& objReg
        );

    //- Member Functions

        virtual scalar value(const label& celli) const = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
