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
    Foam::functionObjects::patchScalarFieldValue

Group
    grpFieldFunctionObjects

Description
    Writes the values of an fvPatchScalarField 

Usage
    Example of function object specification:
    \verbatim
    patchScalarFieldValue1
    {
        type                patchScalarFieldValue;
        libs                ("libfieldFunctionObjects.so");
        ...
        writeToFile             true;
        log                     false;
        field               fieldName
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property      | Description              | Required   | Default value
        type          | type name: patchScalarFieldValue  | yes        |
        writeToFile   | write extents data to file | no       | yes
        log           | write extents data to standard output | no | yes
        field         | field name               | yes        | no
    \endtable

    Output data is written to the file \<timeDir\>/patchScalarFieldValue.dat

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile

SourceFiles
    patchScalarFieldValue.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_patchScalarFieldValue_H
#define functionObjects_patchScalarFieldValue_H

#include "fvMeshFunctionObject.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "writeFile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                        Class patchScalarFieldValue Declaration
\*---------------------------------------------------------------------------*/

class patchScalarFieldValue
:
    public fvMeshFunctionObject,
    public writeFile
{
public: 


protected:

    // Protected data

        //- Name of patch
        word patchName_;

        //- ID of patch
        label patchID_;

        //- Name of field
        word fieldName_;

        //- Field ptr
        const volScalarField* fieldPtr_;

    // Protected Member Functions

        //- Output file header information
        virtual void writeFileHeader(Ostream& os);

        //- No copy construct
        patchScalarFieldValue(const patchScalarFieldValue&) = delete;

        //- No copy assignment
        void operator=(const patchScalarFieldValue&) = delete;


public:

    //- Runtime type information
    TypeName("patchScalarFieldValue");


    // Constructors

        //- Construct from Time and dictionary
        patchScalarFieldValue
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~patchScalarFieldValue() = default;


    // Member Functions

        //- Read the field extents data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the patchScalarFieldValue
        virtual bool write();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
