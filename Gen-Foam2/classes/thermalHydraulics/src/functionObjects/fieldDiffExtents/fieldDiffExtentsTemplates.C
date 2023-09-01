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
#include "volFields.H"
#include "boundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::volScalarField> Foam::functionObjects::fieldDiffExtents::calcMask
(
    const GeometricField<Type, fvPatchField, volMesh>& field,
    const GeometricField<Type, fvPatchField, volMesh>& maskField
) const
{
    dimensionedScalar oneField("", field.dimensions(), 1.0);
    dimensionedScalar oneMaskField("", maskField.dimensions(), 1.0);

    return pos(mag(field/oneField) - mag(maskField/oneMaskField));
}


template<class Type>
void Foam::functionObjects::fieldDiffExtents::calcfieldDiffExtents
(
    const word& fieldName,
    const word& maskFieldName
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;

    const VolFieldType* fieldPtr =
        obr_.findObject<VolFieldType>(fieldName);

    if (!fieldPtr)
    {
        return;
    }

    const VolFieldType* maskFieldPtr =
        obr_.findObject<VolFieldType>(maskFieldName);

    if (!maskFieldPtr)
    {
        return;
    }

    auto extents = [this](const scalarField& mask, const vectorField& C)
    {
        boundBox extents(boundBox::invertedBox);
        forAll(mask, i)
        {
            if (mask[i] > 0.5)
            {
                extents.add(C[i] - C0_);
            }
        };

        extents.reduce();

        if (extents.empty())
        {
            extents.add(point::zero);
        }

        return extents;
    };

    Log << "field: " << fieldName << nl;

    file() << mesh_.time().timeName() << endl;

    tmp<volScalarField> tmask = calcMask<Type>(*fieldPtr, *maskFieldPtr);
    const volScalarField& mask = tmask();

    // Internal field
    if (internalField_)
    {
        boundBox bb(extents(mask, mesh_.C()));
        Log << "    internal field: " << bb << nl;
        file() << bb;

        this->setResult(fieldName + "_internal_min" , bb.min());
        this->setResult(fieldName + "_internal_max", bb.max());
    }

    // Patches
    for (const label patchi : patchIDs_)
    {
        const fvPatchScalarField& maskp = mask.boundaryField()[patchi];
        boundBox bb(extents(maskp, maskp.patch().Cf()));
        const word& patchName = maskp.patch().name();
        Log << "    patch " << patchName << ": " << bb << nl;
        file() << bb;
        this->setResult(fieldName + "_" + patchName + "_min", bb.min());
        this->setResult(fieldName + "_" + patchName + "_max", bb.max());
    }

    Log << endl;
    file() << endl;
}


// ************************************************************************* //
