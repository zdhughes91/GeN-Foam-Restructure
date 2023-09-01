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

#include "stopIfMaxFieldDiff.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(stopIfMaxFieldDiff, 0);
    addToRunTimeSelectionTable(functionObject, stopIfMaxFieldDiff, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::stopIfMaxFieldDiff::stopIfMaxFieldDiff
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    field1Ptr_(nullptr),
    field2Ptr_(nullptr)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::stopIfMaxFieldDiff::read(const dictionary& dict)
{
    if (fvMeshFunctionObject::read(dict))
    {

        field1Name_ = dict.get<word>("field1");
        field2Name_ = dict.get<word>("field2");

        return true;
    }

    return false;
}


bool Foam::functionObjects::stopIfMaxFieldDiff::execute()
{
    return true;
}


bool Foam::functionObjects::stopIfMaxFieldDiff::write()
{   
    //- Set pointers
    if (field1Ptr_ == nullptr)
    {
        field1Ptr_ = &mesh_.lookupObject<volScalarField>(field1Name_);
    }
    if (field2Ptr_ == nullptr)
    {
        field2Ptr_ = &mesh_.lookupObject<volScalarField>(field2Name_);
    }
    const volScalarField& field1(*field1Ptr_);
    const volScalarField& field2(*field2Ptr_);

    //- Check
    if (max(field1-field2).value() > 0.0)
    {
        Log << "    Terminated by " << type() << " " << name() << " at time = "
            << mesh_.time().timeName() << " s" << endl;

        Time& time(const_cast<Time&>(mesh_.time()));
        time.writeAndEnd();
    }

    Log << endl;

    return true;
}


// ************************************************************************* //
