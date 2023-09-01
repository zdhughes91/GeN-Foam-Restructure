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
#if defined __has_include
#  if __has_include(<commDataLayer.H>) 
#    include <commDataLayer.H>
#    define isCommDataLayerIncluded
#  endif
#endif

#ifdef isCommDataLayerIncluded

#include "fieldIntegralToFMU.H"
#include "addToRunTimeSelectionTable.H"
#include "commDataLayer.H"
#include "externalIOObject.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
/*
// public fvMeshFunctionObject  //This would also work
// but one would have to put the functionObject in
// the controlDict instead of the externalCouplingDict
// which goes agaist the logif of FMU4FOAM
namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(fieldIntegralToFMU, 0);
    addToRunTimeSelectionTable(functionObject, fieldIntegralToFMU, dictionary);
}
}
*/

namespace Foam
{
namespace externalIOObject
{
    defineTypeNameAndDebug(fieldIntegralToFMU, 0);
    addToRunTimeSelectionTable(externalIOObject, fieldIntegralToFMU, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::externalIOObject::fieldIntegralToFMU::fieldIntegralToFMU
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    externalIOObject(name, runTime, dict),
    fieldName_(dict.get<word>("fieldName")),
    mesh_
    (
        refCast<const fvMesh>
        (
            time_.lookupObject<objectRegistry>
            (
                dict.getOrDefault("region", polyMesh::defaultRegion)
            )
        )
    ),
    cellZone_(dict.get<word>("cellZone")),
    nameFMU_(dict.get<word>("nameFMU")),    
    fieldPtr_(nullptr)
{
    read(dict);
    execute();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool Foam::externalIOObject::fieldIntegralToFMU::read(const dictionary& dict)
{

    commDataLayer& data = commDataLayer::New(time_);

    data.storeObj(0.0,nameFMU_,commDataLayer::causality::out);
    
    return false;
}

bool Foam::externalIOObject::fieldIntegralToFMU::execute()
{
    commDataLayer& data = commDataLayer::New(time_);
    
    scalar& result = data.getObj<scalar>(nameFMU_,commDataLayer::causality::out);

    const volScalarField& field = mesh_.lookupObject<volScalarField>(fieldName_);

    label cellZoneID = mesh_.cellZones().findZoneID(cellZone_);
    const cellZone& tgtCellZone = mesh_.cellZones()[cellZoneID];

    scalarField fieldZone(field,tgtCellZone);
    scalarField  volZone(mesh_.V(),tgtCellZone);

    result = gSum(fieldZone * volZone);

    return false;
}

bool Foam::externalIOObject::fieldIntegralToFMU::write()
{ 
    return false;
}

#endif
// ************************************************************************* //
