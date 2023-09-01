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
    Foam::subCooledBoilingFractionModel

Description
    Model for the prediction of the sub-cooled bulk liquid temperature at 
    which bubble detachment occurs (only for sub-cooled boiling)

    Author: 
    Stefan Radman (sradman@pm.me / stefan.radman@epfl.ch)

SourceFiles
    subCooledBoilingFractionModel.C
    newsubCooledBoilingFractionModel.C

\*---------------------------------------------------------------------------*/

#ifndef subCooledBoilingFractionModel_H
#define subCooledBoilingFractionModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fluid.H"
#include "structure.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

class FSPair;
class FFPair;

/*---------------------------------------------------------------------------*\
                           Class subCooledBoilingFractionModel Declaration
\*---------------------------------------------------------------------------*/

class subCooledBoilingFractionModel
:
    public IOdictionary
{
protected:

    //- Ref to pair
    const FSPair& pair_;

    //- Ptr to FFPair (fUcK eNcApSuLaTi1o0N!11!1)
    mutable const FFPair* FFPairPtr_;

public:

    TypeName("subCooledBoilingFractionModel");

    //- Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            subCooledBoilingFractionModel,
            subCooledBoilingFractionModels,
            (
                const FSPair& pair,
                const dictionary& dict,
                const objectRegistry& objReg
            ),
            (pair, dict, objReg)
        );

    //- Constructors

        subCooledBoilingFractionModel
        (
            const FSPair& pair,
            const dictionary& dict,
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~subCooledBoilingFractionModel(){}

    //- Selectors

        static autoPtr<subCooledBoilingFractionModel> New
        (
            const FSPair& pair,
            const dictionary& dict,
            const objectRegistry& objReg
        );

    //- Member Functions

        //-
        void setPtr() const;

        //- Return value of the fraction of the pool boiling heat flux
        //  that is resulting in a net vapour generation
        virtual scalar value
        (
            const label& celli,
            const scalar& qi
        ) const = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
