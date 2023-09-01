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
    Foam::regimeMapModel

Description
    Run-time selectable class to handle flow regime map models. 

    Author: 
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

SourceFiles
    regimeMapModel.C

\*---------------------------------------------------------------------------*/

#ifndef regimeMapModel_H
#define regimeMapModel_H

#include "dictionary.H"
#include "dimensionedScalar.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "HashPtrTable.H"
#include "runTimeSelectionTables.H"
#include "myOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class regimeMapModel Declaration
\*---------------------------------------------------------------------------*/

class regimeMapModel
    : public IOdictionary
{
protected :

    //- Typedefs
    typedef HashTable
    <
        label,
        word,
        word::hash
    >   labelTable;

    //- Ref to mesh
    const fvMesh& mesh_;

    //- Table of regimes indexed by regime name
    //regimeTable regimes_;

    //-
    labelTable regimeNameToLabel_;
    wordList regimeLabelToName_;

    //- For eaech mesh cell, this consists of a list of tuples of parameters
    //  for each regime that exists in said cell, wherein the
    //  first element in each tuple is the regime label and the second the
    //  regime coefficient. It is handled via DynamicLists as they are
    //  significantly faster to allocate, append to, and clear than regular
    //  lists (which is something that is done quite often for this object)
    DynamicList<DynamicList<Tuple2<label,scalar>>> regimeLabelCoeffs_;

public:

    //- Runtime type information
    TypeName("regimeMapModel");

    //- Runtime selection table declaration
    declareRunTimeSelectionTable
    (
        autoPtr,
        regimeMapModel,
        regimeMapModels,
        (
            const fvMesh& mesh,
            const dictionary& dict
        ),
        (mesh, dict)
    );

    //- Constructors

        regimeMapModel
        (
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Destructor
        virtual ~regimeMapModel(){}

    //- Selectors
    
        //- The selector does not have to be the same as the constructor
        static autoPtr<regimeMapModel> New
        (
            const fvMesh& mesh,
            const dictionary& dict
        );

    //- Members

        //- Updated the marker field of each regime in regimeTable
        //  regimeTable according to the desired map model
        virtual void correct() = 0;

        //- Fill in the models list based on the different models specified
        //  for each flow regime
        template<class modelType, class firstArgType>
        void constructModels
        (
            List<autoPtr<modelType>>& models,
            const firstArgType& firstArg,
            const dictionary& dict,
            const objectRegistry& objReg
        ) const;

        //- If different models in the models list use different correlations
        //  to predict a value for a certain quantity, this function returns
        //  the interpolated value (for the mesh cell in exam) between the 
        //  different models, as each model corrseponds to a certain flow 
        //  regime
        template<class valueType, class modelType>
        valueType interpolateValue
        (
            const List<autoPtr<modelType>>& models, 
            const label& celli
        ) const;

        //- Access

            int size() const
            {
                return regimeNameToLabel_.size();
            }

            const labelTable& regimeNameToLabel() const
            {
                return regimeNameToLabel_;
            }

            const wordList& regimeLabelToName() const
            {
                return regimeLabelToName_;
            }

            const DynamicList<DynamicList<Tuple2<label,scalar>>>& 
            regimeLabelCoeffs() const
            {
                return regimeLabelCoeffs_;
            }

    //- Static member functions

        static const regimeMapModel& lookupRegimeMap
        (
            const objectRegistry& objReg, 
            const word& name
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "regimeMapModelTemplates.C"
#endif

#endif

// ************************************************************************* //
