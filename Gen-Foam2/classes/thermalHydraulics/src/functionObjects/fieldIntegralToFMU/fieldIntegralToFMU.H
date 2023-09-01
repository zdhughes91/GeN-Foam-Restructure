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
    Foam::externalIOObject::fieldIntegralToFMU

Description
    Output an integral over a cellZone to an FMU

Usage
    Minimal example by using \c system/externalCouplingDict:
    \verbatim
    fieldIntegralToFMU
    {
        type fieldIntegralToFMU;
        nameFMU hxPowerOut;
        fieldName heatFlux.structure;
        region     fluidRegion;
        cellZone hx;
    }  
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                   | Type | Req'd | Dflt
      type         | Type name: fieldIntegralToFMU       | word |  yes  | -
      nameFMU     | name of the output scalar    | word |  yes  | -
      fieldName     | name of the measurement field | word |  yes  | -
      region      | name of the region         | word |  yes  | -
      cellZone      | position of the cellZone        | word |  yes  | -
    \endtable

SourceFiles
    fieldIntegralToFMU.C

\*---------------------------------------------------------------------------*/

#ifndef externalIOObject_fieldIntegralToFMU_H
#define externalIOObject_fieldIntegralToFMU_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "volFieldSelection.H"
#include "rhoThermo.H"
#include "IOobjectList.H"
#include "faceSet.H"
#include "surfaceFields.H"
#include "fvc.H"

#include "externalIOObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace externalIOObject //functionObjects
{

/*---------------------------------------------------------------------------*\
                        Class fieldIntegralToFMU Declaration
\*---------------------------------------------------------------------------*/

class fieldIntegralToFMU
:
    public externalIOObject
    // public fvMeshFunctionObject  //This would also work
    // but one would have to put the functionObject in
    // the controlDict instead of the externalCouplingDict
    // which goes agaist the logif of FMU4FOAM  
{
public: 

protected:

    // Protected data
        
        //- Name of fields
        word fieldName_;

        //- Mesh
        const fvMesh& mesh_;

        //- Name of cellZone
        word cellZone_;

        //- Name of scalar to be exhanged with the FMU
        word nameFMU_;

        //- Field ptrs
        const volScalarField* fieldPtr_;


    // Protected Member Functions

        //- No copy construct
        fieldIntegralToFMU(const fieldIntegralToFMU&) = delete;

        //- No copy assignment
        void operator=(const fieldIntegralToFMU&) = delete;


public:

    //- Runtime type information
    TypeName("fieldIntegralToFMU");


    // Constructors

        //- Construct from Time and dictionary
        fieldIntegralToFMU
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~fieldIntegralToFMU() = default;


    // Member Functions

        //- Read the field extents data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the fieldIntegralToFMU
        virtual bool write();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace externalIOObject
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
