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
    Foam::saturationModel

Description

SourceFiles
    saturationModel.C
    newSaturationModel.C

\*---------------------------------------------------------------------------*/

#ifndef saturationModel_H
#define saturationModel_H

#include <map>
#include "fluid.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class FFPair;
class phaseChangeModel;

/*---------------------------------------------------------------------------*\
                         Class saturationModel Declaration
\*---------------------------------------------------------------------------*/

class saturationModel
:
    public IOdictionary
{
protected:

    //- Ref to mesh
    const fvMesh& mesh_;

    //- Ref to liquid
    const fluid& liquid_;

    //- Ref to vapour
    const fluid& vapour_;

    //- Pointer to a function (of the scalar value...(const label&) const)
    scalar (Foam::saturationModel::*funcPtr_)(const label&) const;

    //- Map that associates a function ptr to a word
    std::map<word, scalar (Foam::saturationModel::*)(const label&) const> map_;

public:

    TypeName("saturationModel");

    //- Declare runtime construction
    declareRunTimeSelectionTable
    (
        autoPtr,
        saturationModel,
        saturationModels,
        (
            const phaseChangeModel& pcm,
            const dictionary& dict, 
            const objectRegistry& objReg
        ),
        (pcm, dict, objReg)
    );

    //- Constructors

        saturationModel
        (
            const phaseChangeModel& pcm,
            const dictionary& dict, 
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~saturationModel(){}

    //- Selectors

        static autoPtr<saturationModel> New
        (
            const phaseChangeModel& pcm,
            const dictionary& dict, 
            const objectRegistry& objReg
        );

    //- Member Functions

        virtual scalar valueTSat(const label& celli) const = 0;

        virtual scalar valuePSat(const label& celli) const = 0;

        virtual scalar valueLnPSat(const label& celli) const = 0;

        virtual scalar valuePSatPrime(const label& celli) const = 0;

        //-
        virtual void correctField
        (
            volScalarField& field, 
            const word& func
        );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
