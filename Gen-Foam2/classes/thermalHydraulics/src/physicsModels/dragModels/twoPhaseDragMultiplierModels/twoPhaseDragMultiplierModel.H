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
    Foam::twoPhaseDragMultiplierModel

Description

    Author: 
    Stefan Radman (sradman@protonmail.com / stefan.radman@epfl.ch)

SourceFiles
    twoPhaseDragMultiplierModel.C
    newTwoPhaseDragMultiplierModel.C

\*---------------------------------------------------------------------------*/

#ifndef twoPhaseDragMultiplierModel_H
#define twoPhaseDragMultiplierModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "fluid.H"
#include "FSPair.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class twoPhaseDragMultiplierModel Declaration
\*---------------------------------------------------------------------------*/

class twoPhaseDragMultiplierModel
:
    public IOdictionary
{
protected:

    //- Ref to mesh
    const fvMesh& mesh_;

    //- Ptr to fluid whose properties are used as basis for the calculation of
    //  the two phase pressure drop
    const fluid* mFluidPtr_;

    //- Ptr to other fluid in the system
    const fluid* oFluidPtr_;

    //- Ptr to multiplier fluid (single-phase) friction factor tensor
    const volTensorField* mKd_;

    //- Maximum allowable value of phi2
    scalar maxPhi2_;

public:

    //- Runtime type information
    TypeName("twoPhaseDragMultiplierModel");

    //- Declare runtime construction
        
        declareRunTimeSelectionTable
        (
            autoPtr,
            twoPhaseDragMultiplierModel,
            twoPhaseDragMultiplierModels,
            (
                const fvMesh& mesh,
                const dictionary& dict,
                const objectRegistry& objReg
            ),
            (mesh, dict, objReg)
        );

    //- Constructors

        twoPhaseDragMultiplierModel
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~twoPhaseDragMultiplierModel(){}

    //- Selectors

        static autoPtr<twoPhaseDragMultiplierModel> New
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const objectRegistry& objReg
        );

    //- Member Functions

        //- True if flow is essentially single phase in celli
        bool onePhase(const label& celli) const;

        //- Return multiplier
        virtual scalar phi2(const label& celli) const;

        //- Return total drag factor tensor divided by
        //  multiplier fluid velocity.
        virtual tensor value
        (
            const label& celli
        ) const;

        //-
        void correctField(volTensorField& KdTotByU);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
