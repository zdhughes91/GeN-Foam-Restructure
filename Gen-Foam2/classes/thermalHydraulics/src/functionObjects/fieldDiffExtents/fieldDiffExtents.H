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
    Foam::functionObjects::fieldDiffExtents

Group
    grpFieldFunctionObjects

Description
    Calculates the spatial minimum and maximum extents of a difference between
    two fields.

    The extents are derived from the bound box limits after identifying the
    locations where one field exceeds the maskField

Usage
    Example of function object specification:
    \verbatim
    fieldDiffExtents1
    {
        type        fieldDiffExtents;
        libs        ("libFFSSrc.so");
        ...
        writeToFile yes;
        log         yes;
        fields      (T.structure);
        maskFields  (T.interface);
        patches     ();
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property      | Description                 | Required  | Default value
        type          | type name: fieldDiffExtents | yes       |
        writeToFile   | write extents data to file  | no        | yes
        log           | write extents data to standard output|no| yes
        internalField | Process the internal field  | no        | yes
        threshold     | Field value to identify extents boundary| yes 
        referencePosition | Reference position      | no        | (0 0 0)
        fields        | list of fields to process   | yes       |
        maskFields    | list of mask fields         | yes       |
        patches       | list of patches to process  | no      | \<all patches\>
    \endtable

    Output data is written to the file \<timeDir\>/fieldDiffExtents.dat

Note
    For non-scalar fields, the magnitude of the field is employed and compared
    to the mask field magnitude

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile

SourceFiles
    fieldDiffExtents.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_fieldDiffExtents_H
#define functionObjects_fieldDiffExtents_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "volFields.H"
#include "processorPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                        Class fieldDiffExtents Declaration
\*---------------------------------------------------------------------------*/

class fieldDiffExtents
:
    public fvMeshFunctionObject,
    public writeFile
{
protected:

    // Protected data

        //- Flag to write the internal field extents
        bool internalField_;

        //- Reference position; default = (0 0 0)
        point C0_;

        //- Fields to assess
        wordList fieldNames_;

        //- Mask fields
        wordList maskFieldNames_;

        //- Patches to assess
        labelHashSet patchIDs_;


    // Protected Member Functions

        //- Output file header information
        virtual void writeFileHeader(Ostream& os);

        //- Return the field mask
        template<class Type>
        tmp<volScalarField> calcMask
        (
            const GeometricField<Type, fvPatchField, volMesh>& field,
            const GeometricField<Type, fvPatchField, volMesh>& maskField
        ) const;

        //- Main calculation
        template<class Type>
        void calcfieldDiffExtents
        (
            const word& fieldName,
            const word& maskFieldName
        );

        //- No copy construct
        fieldDiffExtents(const fieldDiffExtents&) = delete;

        //- No copy assignment
        void operator=(const fieldDiffExtents&) = delete;


public:

    //- Runtime type information
    TypeName("fieldDiffExtents");


    // Constructors

        //- Construct from Time and dictionary
        fieldDiffExtents
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~fieldDiffExtents() = default;


    // Member Functions

        //- Read the field extents data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the fieldDiffExtents
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
tmp<volScalarField> fieldDiffExtents::calcMask
(
    const GeometricField<scalar, fvPatchField, volMesh>& field,
    const GeometricField<scalar, fvPatchField, volMesh>& maskField
) const;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fieldDiffExtentsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
