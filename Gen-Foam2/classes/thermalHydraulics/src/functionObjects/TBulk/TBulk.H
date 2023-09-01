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
    Foam::functionObjects::TBulk

Group
    grpFieldFunctionObjects

Description
    Calculates the total mass flux through a set of faceSets, faceZones or
    patches.

Usage
    Example of function object specification:
    \verbatim
    TBulk1
    {
        type                TBulk;
        libs                ("libfieldFunctionObjects.so");
        ...
        writeToFile         yes;
        log                 yes;
        regionType          patch;
        regionName          outlet;
        alphaRhoPhiName     "alphaRhoPhi.liquid"
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property      | Description              | Required   | Default value
        type          | type name: TBulk      | yes        |
        writeToFile   | write extents data to file | no       | yes
        log           | write extents data to standard output | no | yes
        thermoName    | name of fluid thermophysicalProperties dict | yes
        alphaRhoPhiName | name of mass face flow field | yes |
        regionType    | either cellZone or patch | yes |
        regionName    | name of the regionType to process | yes |
    \endtable

    Output data is written to the file \<timeDir\>/TBulk.dat

Note
    For non-scalar fields, the magnitude of the field is employed and compared
    to the threshold value.

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile

SourceFiles
    TBulk.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_TBulk_H
#define functionObjects_TBulk_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "volFieldSelection.H"
#include "rhoThermo.H"
#include "IOobjectList.H"
#include "faceSet.H"
#include "surfaceFields.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                        Class TBulk Declaration
\*---------------------------------------------------------------------------*/

class TBulk
:
    public fvMeshFunctionObject,
    public writeFile
{
public: 

    //-
    enum regionType
    {
        patch,
        faceSet,
        faceZone
    };

    static const Enum<regionType> regionTypeNames_;

protected:

    // Protected data

        //-
        word regionName_;

        //- Patch related
        label patchID_;
    
        //- faceZone or faceSet related
        labelList faces_;

        //- Region type
        regionType regionType_;

        //- Field/object names
        word thermoName_;
        word alphaRhoPhiName_;

        //- Field/object ptrs
        const rhoThermo* thermoPtr_;
        const surfaceScalarField* alphaRhoPhiPtr_;

    // Protected Member Functions

        //- Output file header information
        virtual void writeFileHeader(Ostream& os);

        //- No copy construct
        TBulk(const TBulk&) = delete;

        //- No copy assignment
        void operator=(const TBulk&) = delete;


public:

    //- Runtime type information
    TypeName("TBulk");


    // Constructors

        //- Construct from Time and dictionary
        TBulk
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~TBulk() = default;


    // Member Functions

        //- Read the field extents data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the TBulk
        virtual bool write();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
