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

\*---------------------------------------------------------------------------*/

#include "patchScalarFieldValue.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(patchScalarFieldValue, 0);
    addToRunTimeSelectionTable
    (
        functionObject, 
        patchScalarFieldValue, 
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::patchScalarFieldValue::writeFileHeader(Ostream& os)
{
    /*
    if (writtenHeader_)
    {
        os << endl;
    }
    else
    {
        word headerText
        (
            "Value of volScalarField "+fieldName_+" on patch "+patchName_
        );
        writeHeader(os, headerText);
    }
    
    word time("Time = "+mesh_.time().timeName());

    os << time << endl;
    */

    writtenHeader_ = true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::patchScalarFieldValue::patchScalarFieldValue
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name, typeName, dict),
    fieldPtr_(nullptr)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::patchScalarFieldValue::read(const dictionary& dict)
{
    if (fvMeshFunctionObject::read(dict) && writeFile::read(dict))
    {
        fieldName_ = dict.get<word>("field");

        patchName_ = dict.get<word>("patch");

        const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
        
        patchID_ = pbm.findIndex(patchName_);
        
        return true;
    }

    return false;
}


bool Foam::functionObjects::patchScalarFieldValue::execute()
{
    return true;
}


bool Foam::functionObjects::patchScalarFieldValue::write()
{
    //writeFileHeader(file());

    file() << endl;

    if (fieldPtr_ == nullptr)
    {
        fieldPtr_ = &(mesh_.lookupObject<volScalarField>(fieldName_));
    }
    const volScalarField& field(*fieldPtr_);
    const scalarField& bField(field.boundaryField()[patchID_]);

    const scalar t = mesh_.time().timeOutputValue();

    file() << "( " << t << " ( "; 
    forAll(bField, i)
    {
        file() << bField[i] << " ";
    }
    file() << "))";
    
    return true;
}


// ************************************************************************* //
