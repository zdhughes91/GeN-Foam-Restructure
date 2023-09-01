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

#include "TBulk.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(TBulk, 0);
    addToRunTimeSelectionTable(functionObject, TBulk, dictionary);
}
}

const Foam::Enum
<
    Foam::functionObjects::TBulk::regionType
>
Foam::functionObjects::TBulk::regionTypeNames_
(
    {
        { 
            regionType::patch, 
            "patch" 
        },
        { 
            regionType::faceSet, 
            "faceSet" 
        },
        { 
            regionType::faceZone, 
            "faceZone" 
        }
    }
);


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::TBulk::writeFileHeader(Ostream& os)
{
    if (writtenHeader_)
    {
        writeBreak(os);
    }
    else
    {
        writeHeader(os, "Field extents");
    }

    writeCommented(os, "Time");

    os  << endl;

    writtenHeader_ = true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::TBulk::TBulk
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name, typeName, dict),
    regionName_(""),
    patchID_(0),
    faces_(0),
    thermoPtr_(nullptr),
    alphaRhoPhiPtr_(nullptr)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::TBulk::read(const dictionary& dict)
{
    if (fvMeshFunctionObject::read(dict) && writeFile::read(dict))
    {
        regionType_ = regionType
        (
            regionTypeNames_.get
            (
                dict.get<word>
                (
                    "regionType"
                )
            )
        );
        regionName_ = dict.get<word>("regionName");
        thermoName_ = dict.get<word>("thermoName");
        alphaRhoPhiName_ = dict.get<word>("alphaRhoPhiName");

        if (regionType_ == regionType::patch)
        {
            const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
            patchID_ = pbm.findPatchID(regionName_);
        }
        else if (regionType_ == regionType::faceZone)
        {
            const faceZoneMesh& faceZones(mesh_.faceZones());
            const labelList& faceis(faceZones[regionName_]);
            forAll(faceis, i)
            {
                faces_.append(faceis[i]);
            }
        }
        else if (regionType_ == regionType::faceSet)
        {
            IOobjectList objects
            (
                mesh_,
                mesh_.time().findInstance
                (
                    polyMesh::meshSubDir/"sets",
                    word::null,
                    IOobject::READ_IF_PRESENT,
                    mesh_.facesInstance()
                ),
                polyMesh::meshSubDir/"sets"
            );
            IOobjectList faceSets(objects.lookupClass(faceSet::typeName));
            if (faceSets.found(regionName_))
            {
                Foam::faceSet set(*faceSets[regionName_]);
                forAllIter(faceSet, set, iter)
                {
                    label facei(*iter);
                    faces_.append(facei);
                }
            }
        }

        return true;
    }

    return false;
}


bool Foam::functionObjects::TBulk::execute()
{
    return true;
}


bool Foam::functionObjects::TBulk::write()
{
    writeFileHeader(file());

    Log << type() << " " << name() <<  " write:" << nl;

    //- Set pointers
    if (thermoPtr_ == nullptr)
    {
        thermoPtr_ = &mesh_.lookupObject<rhoThermo>(thermoName_);
    }
    if (alphaRhoPhiPtr_ == nullptr)
    {
        alphaRhoPhiPtr_ = 
            &mesh_.lookupObject<surfaceScalarField>(alphaRhoPhiName_);
    }

    const rhoThermo& thermo(*thermoPtr_);
    const volScalarField& T(thermo.T());
    tmp<volScalarField> Cp(thermo.Cp());
    const surfaceScalarField& alphaRhoPhi(*alphaRhoPhiPtr_);

    scalar hDot(0.0);
    scalar hDotByT(0.0);

    if (regionType_ == regionType::patch)
    {
        const fvPatchScalarField& Cpp = Cp().boundaryField()[patchID_];
        const fvPatchScalarField& Tp = T.boundaryField()[patchID_];
        const fvsPatchField<scalar>& alphaRhoPhip 
            = alphaRhoPhi.boundaryField()[patchID_];
        const fvPatch& patch(mesh_.boundary()[patchID_]);
        const scalarField& magSf(patch.magSf());
        forAll(magSf, i)
        {
            scalar magAlphaRhoCpPhipi(mag(alphaRhoPhip[i])*Cpp[i]);
            hDotByT += magAlphaRhoCpPhipi;
            hDot += magAlphaRhoCpPhipi*Tp[i];
        }
    }
    else
    {
        surfaceScalarField Tf(fvc::interpolate(T));
        surfaceScalarField Cpf(fvc::interpolate(Cp));
        forAll(faces_, i)
        {
            const label& facei(faces_[i]);
            scalar magAlphaRhoCpPhii(mag(alphaRhoPhi[facei])*Cpf[facei]);
            hDotByT += magAlphaRhoCpPhii;
            hDot += magAlphaRhoCpPhii*Tf[facei];
        }
    }

    reduce(hDotByT, sumOp<scalar>());
    reduce(hDot, sumOp<scalar>());

    scalar Tb(hDot/max(hDotByT, 1e-9));

    Log << "    " << regionTypeNames_[regionType_] << " " << regionName_ 
        << " TBulk = " << Tb << " K" << endl;
    file() << Tb;
    this->setResult(regionName_+"_TBulk", Tb);

    Log << endl;

    return true;
}


// ************************************************************************* //
