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
    Foam::contactPartitionModel

Description
    Class to handle fluid-structure interfacial contact contact fractions

    Author: 
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

SourceFiles
    contactPartitionModel.C
    newContactPartitionModel.C

\*---------------------------------------------------------------------------*/

#ifndef contactPartitionModel_H
#define contactPartitionModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "fluid.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class FSPair;

/*---------------------------------------------------------------------------*\
            Class contactPartitionModel Declaration
\*---------------------------------------------------------------------------*/

class contactPartitionModel
:
    public IOdictionary
{
protected:

    const fvMesh& mesh_;

    const FSPair& pair_;

    const fluid& fluid_;

public:

    //- Runtime type information
    TypeName("contactPartitionModel");

    //- Declare runtime construction
        
        declareRunTimeSelectionTable
        (
            autoPtr,
            contactPartitionModel,
            contactPartitionModels,
            (
                const FSPair& pair,
                const dictionary& dict,
                const objectRegistry& objReg
            ),
            (pair, dict, objReg)
        );

    // Constructors

        contactPartitionModel
        (
            const FSPair& pair,
            const dictionary& dict,
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~contactPartitionModel();

    // Selectors

        static autoPtr<contactPartitionModel> New
        (
            const FSPair& pair,
            const dictionary& dict,
            const objectRegistry& objReg
        );

    // Member Functions

        //- Return fraction of structure interfacial area in
        //  contact with the fluid
        virtual scalar value(const label& celli) const = 0;

        //-
        virtual void correctField(volScalarField& f);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
