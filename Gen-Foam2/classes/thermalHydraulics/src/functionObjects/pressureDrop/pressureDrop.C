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

#include "pressureDrop.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(pressureDrop, 0);
    addToRunTimeSelectionTable(functionObject, pressureDrop, dictionary);
}
}

const Foam::Enum
<
    Foam::functionObjects::pressureDrop::regionType
>
Foam::functionObjects::pressureDrop::regionTypeNames_
(
    {
        { 
            regionType::patch, 
            "patch" 
        },
        { 
            regionType::faceSet_, 
            "faceSet" 
        }
    }
);


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::pressureDrop::writeFileHeader(Ostream& os)
{
    if (writtenHeader_)
    {
        os << endl;
    }
    else
    {
        word headerText
        (
            "Pressure drop between "+region1Name_+" and "+region2Name_
        );
        writeHeader(os, headerText);
    }
    
    word time("Time = "+mesh_.time().timeName());

    os << time << endl;

    writtenHeader_ = true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::pressureDrop::pressureDrop
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name, typeName, dict),
    faces1_(0),
    faces2_(0),
    S1_(0),
    S2_(0),
    pPtr_(nullptr)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::pressureDrop::read(const dictionary& dict)
{
    if (fvMeshFunctionObject::read(dict) && writeFile::read(dict))
    {
        region1Type_ = regionType
        (
            regionTypeNames_.get
            (
                dict.get<word>
                (
                    "region1Type"
                )
            )
        );
        region2Type_ = regionType
        (
            regionTypeNames_.get
            (
                dict.get<word>
                (
                    "region2Type"
                )
            )
        );
        region1Name_ = dict.get<word>("region1");
        region2Name_ = dict.get<word>("region2");

        pName_ = dict.lookupOrDefault<word>("pressure", "p_rgh");

        if (region1Type_ == regionType::patch)
        {
            const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
            patchID1_ = pbm.findIndex(region1Name_);
            if (S1_ == 0.0)
            {
                const fvPatch& patch(mesh_.boundary()[patchID1_]);
                const scalarField& magSf(patch.magSf());
                forAll(magSf, i)
                {
                    S1_ += magSf[i];
                }
                reduce(S1_, sumOp<scalar>());
            }
        }
        else if (region1Type_ == regionType::faceSet_)
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
            faceSet set(*faceSets[region1Name_]);
            forAllIter(faceSet, set, iter)
            {
                label facei(*iter);
                faces1_.append(facei);
            }
            if (S1_ == 0.0)
            {
                const scalarField& magSf(mesh_.magSf());
                forAll(faces1_, i)
                {
                    S1_ += magSf[faces1_[i]];
                }
                reduce(S1_, sumOp<scalar>());
            }
        }
        
        if (region2Type_ == regionType::patch)
        {
            const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
            patchID2_ = pbm.findIndex(region2Name_);
            if (S2_ == 0.0)
            {
                const fvPatch& patch(mesh_.boundary()[patchID2_]);
                const scalarField& magSf(patch.magSf());
                forAll(magSf, i)
                {
                    S2_ += magSf[i];
                }
                reduce(S2_, sumOp<scalar>());
            }
        }
        else if (region2Type_ == regionType::faceSet_)
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
            faceSet set(*faceSets[region2Name_]);
            forAllIter(faceSet, set, iter)
            {
                label facei(*iter);
                faces2_.append(facei);
            }
            if (S2_ == 0.0)
            {
                const scalarField& magSf(mesh_.magSf());
                forAll(faces2_, i)
                {
                    S2_ += magSf[faces2_[i]];
                }
                reduce(S2_, sumOp<scalar>());
            }
        }

        return true;
    }

    return false;
}


bool Foam::functionObjects::pressureDrop::execute()
{
    return true;
}


bool Foam::functionObjects::pressureDrop::write()
{
    writeFileHeader(file());

    Log << type() << " " << name() <<  " write:" << nl;

    if (pPtr_ == nullptr)
    {
        pPtr_ = &mesh_.lookupObject<volScalarField>(pName_);
    }
    const volScalarField& p(*pPtr_);
    surfaceScalarField pi(fvc::interpolate(p));

    //- Pressures at regions 1, 2
    scalar p1(0);
    scalar p2(0);

    switch (region1Type_)
    {
        case regionType::patch :
        {               
            const fvPatchScalarField& pp = p.boundaryField()[patchID1_];
            const fvPatch& patch(mesh_.boundary()[patchID1_]);
            const scalarField& magSf(patch.magSf());
            forAll(magSf, i)
            {
                p1 += pp[i]*magSf[i];
            }
            reduce(p1, sumOp<scalar>());
            p1 /= S1_;
            break;
        }
        case regionType::faceSet_ :
        {
            const scalarField& magSf(mesh_.magSf());
            forAll(faces1_, i)
            {
                const label& facei(faces1_[i]);
                p1 += pi[facei]*magSf[facei];
            }
            reduce(p1, sumOp<scalar>());
            p1 /= S1_;
            break;
        }
    }

    switch (region2Type_)
    {
        case regionType::patch :
        {               
            const fvPatchScalarField& pp = p.boundaryField()[patchID2_];
            const fvPatch& patch(mesh_.boundary()[patchID2_]);
            const scalarField& magSf(patch.magSf());
            forAll(magSf, i)
            {
                p2 += pp[i]*magSf[i];
            }
            reduce(p2, sumOp<scalar>());
            p2 /= S2_;
            break;
        }
        case regionType::faceSet_ :
        {
            const scalarField& magSf(mesh_.magSf());
            forAll(faces2_, i)
            {
                const label& facei(faces2_[i]);
                p2 += pi[facei]*magSf[facei];
            }
            reduce(p2, sumOp<scalar>());
            p2 /= S2_;
            break;
        }
    }

    scalar deltaP(p1-p2);

    Log << "    " << region1Name_ << " "<< region2Name_ << " deltaP = " <<
        deltaP << " Pa" << endl;
    file() << "deltaP = " << deltaP << " Pa";

    Log << endl;

    return true;
}


// ************************************************************************* //
