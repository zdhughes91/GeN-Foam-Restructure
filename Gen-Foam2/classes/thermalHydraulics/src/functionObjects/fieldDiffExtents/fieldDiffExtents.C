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

#include "fieldDiffExtents.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(fieldDiffExtents, 0);
    addToRunTimeSelectionTable(functionObject, fieldDiffExtents, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::fieldDiffExtents::writeFileHeader(Ostream& os)
{
    if (writtenHeader_)
    {
        writeBreak(os);
    }
    else
    {
        writeHeader(os, "Field extents");
        writeHeaderValue(os, "Reference position", C0_);
    }

    writeCommented(os, "Time");

    forAll(fieldNames_, i)
    {
        word fieldName(fieldNames_[i]);
        if (internalField_)
        {
            writeTabbed(os, fieldName + "_internal");
        }
        for (const label patchi : patchIDs_)
        {
            const word& patchName = mesh_.boundaryMesh()[patchi].name();
            writeTabbed(os, fieldName + "_" + patchName);
        }
    }

    forAll(maskFieldNames_, i)
    {
        word fieldName(maskFieldNames_[i]);
        if (internalField_)
        {
            writeTabbed(os, fieldName + "_internal");
        }
        for (const label patchi : patchIDs_)
        {
            const word& patchName = mesh_.boundaryMesh()[patchi].name();
            writeTabbed(os, fieldName + "_" + patchName);
        }
    }

    os  << endl;

    writtenHeader_ = true;
}


template<>
Foam::tmp<Foam::volScalarField> Foam::functionObjects::fieldDiffExtents::calcMask
(
    const GeometricField<scalar, fvPatchField, volMesh>& field,
    const GeometricField<scalar, fvPatchField, volMesh>& maskField
) const
{
    dimensionedScalar oneField("", field.dimensions(), 1.0);
    dimensionedScalar oneMaskField("", maskField.dimensions(), 1.0);

    return pos((field/oneField) - (maskField/oneMaskField));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fieldDiffExtents::fieldDiffExtents
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name, typeName, dict),
    internalField_(true),
    C0_(Zero),
    fieldNames_(0),
    maskFieldNames_(0),
    patchIDs_()
{
    read(dict);

    // Note: delay creating the output file header to handle field names
    // specified using regular expressions
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::fieldDiffExtents::read(const dictionary& dict)
{
    if (fvMeshFunctionObject::read(dict) && writeFile::read(dict))
    {
        dict.readIfPresent<bool>("internalField", internalField_);

        dict.readIfPresent<vector>("referencePosition", C0_);

        patchIDs_.clear();
        const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

        wordReList patchNames;
        if (dict.readIfPresent("patches", patchNames))
        {
            for (const wordRe& name : patchNames)
            {
                patchIDs_.insert(pbm.findIndices(name));
            }
        }
        else
        {
            // Add all non-processor and non-empty patches
            forAll(pbm, patchi)
            {
                const polyPatch& pp = pbm[patchi];
                if (!isA<processorPolyPatch>(pp) && !isA<emptyPolyPatch>(pp))
                {
                    patchIDs_.insert(patchi);
                }
            }
        }

        if (!internalField_ && patchIDs_.empty())
        {
            IOWarningInFunction(dict)
                << "No internal field or patches selected - no field extent "
                << "information will be generated" << endl;
        }

        fieldNames_ = dict.get<wordList>("fields");
        maskFieldNames_ = dict.get<wordList>("maskFields");

        if (fieldNames_.size() != maskFieldNames_.size())
        {
            FatalErrorInFunction
            << "The fields and maskFields lists must have a 1-1 correspondence, "
            << "i.e. the same size" << exit(FatalError);
        }

        return true;
    }

    return false;
}


bool Foam::functionObjects::fieldDiffExtents::execute()
{
    return true;
}


bool Foam::functionObjects::fieldDiffExtents::write()
{
    writeFileHeader(file());

    Log << type() << " " << name() <<  " write:" << nl;

    forAll(fieldNames_, i)
    {
        word fieldName(fieldNames_[i]);
        word maskFieldName(maskFieldNames_[i]);
        calcfieldDiffExtents<scalar>(fieldName, maskFieldName);
        calcfieldDiffExtents<vector>(fieldName, maskFieldName);
        calcfieldDiffExtents<sphericalTensor>(fieldName, maskFieldName);
        calcfieldDiffExtents<symmTensor>(fieldName, maskFieldName);
        calcfieldDiffExtents<tensor>(fieldName, maskFieldName);
    }

    Log << endl;

    return true;
}


// ************************************************************************* //
