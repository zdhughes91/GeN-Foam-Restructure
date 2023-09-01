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
    Foam::functionObjects::massFlow

Group
    grpFieldFunctionObjects

Description
    Calculates the total mass flux through a set of faceSets, faceZones or
    patches.

Usage
    Example of function object specification:
    \verbatim
    massFlow1
    {
        type                massFlow;
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
        type          | type name: massFlow      | yes        |
        writeToFile   | write extents data to file | no       | yes
        log           | write extents data to standard output | no | yes
        alphaRhoPhiName | name of mass face flow field | yes |
        regionType    | either cellZone or patch | yes |
        regionName    | name of the regionType to process | yes |
        scaleFactor   | rescale the output by this factor | no | 1.0
    \endtable

    Output data is written to the file \<timeDir\>/massFlow.dat

Note
    For non-scalar fields, the magnitude of the field is employed and compared
    to the threshold value.

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile

SourceFiles
    massFlow.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_massFlow_H
#define functionObjects_massFlow_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "volFieldSelection.H"
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
                        Class massFlow Declaration
\*---------------------------------------------------------------------------*/

class massFlow
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

        //- Field names
        word alphaRhoPhiName_;

        //- Field ptrs
        const surfaceScalarField* alphaRhoPhiPtr_;

        //- Optional scale factor
        scalar scaleFactor_;

    // Protected Member Functions

        //- Output file header information
        virtual void writeFileHeader(Ostream& os);

        //- No copy construct
        massFlow(const massFlow&) = delete;

        //- No copy assignment
        void operator=(const massFlow&) = delete;


public:

    //- Runtime type information
    TypeName("massFlow");


    // Constructors

        //- Construct from Time and dictionary
        massFlow
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~massFlow() = default;


    // Member Functions

        //- Read the field extents data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the massFlow
        virtual bool write();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
