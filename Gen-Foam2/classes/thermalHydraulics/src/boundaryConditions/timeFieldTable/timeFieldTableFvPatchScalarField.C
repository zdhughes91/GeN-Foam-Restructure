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

#include "timeFieldTableFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "one.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeFieldTableFvPatchScalarField::
timeFieldTableFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF),
    currentIndex_(0),
    table_(),
    tStart_(0)
{}


Foam::timeFieldTableFvPatchScalarField::
timeFieldTableFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict, false),
    currentIndex_(0),
    table_
    (
        dict.get
        <
            List
            <
                Tuple2
                <
                    scalar, 
                    scalarField
                >
            >
        >("table")
    ),
    tStart_(table_[0].first())
{    
    updateField();
}


Foam::timeFieldTableFvPatchScalarField::
timeFieldTableFvPatchScalarField
(
    const timeFieldTableFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
    currentIndex_(ptf.currentIndex_),
    table_(ptf.table_),
    tStart_(ptf.tStart_)
{}


Foam::timeFieldTableFvPatchScalarField::
timeFieldTableFvPatchScalarField
(
    const timeFieldTableFvPatchScalarField& ptf
)
:
    fixedValueFvPatchField<scalar>(ptf),
    currentIndex_(ptf.currentIndex_),
    table_(ptf.table_),
    tStart_(ptf.tStart_)
{}


Foam::timeFieldTableFvPatchScalarField::
timeFieldTableFvPatchScalarField
(
    const timeFieldTableFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(ptf, iF),
    currentIndex_(ptf.currentIndex_),
    table_(ptf.table_),
    tStart_(ptf.tStart_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::timeFieldTableFvPatchScalarField::updateField()
{
    //- Simulation time
    const scalar t = db().time().timeOutputValue();

    if (t > tStart_)
    {
        scalar t0(0);
        scalar t1(0);
        if (currentIndex_ < (table_.size()-1))
        {
            t0 = table_[currentIndex_].first();
            t1 = table_[currentIndex_+1].first();
            while (true)
            {
                if ((t >= t0 and t <= t1))
                    break;
                currentIndex_++;
                t0 = table_[currentIndex_].first();
                t1 = table_[currentIndex_+1].first(); 
            }
        }

        if (currentIndex_ < (table_.size()-1))
        {
            //- Read table values at correct index
            scalarField& f0(table_[currentIndex_].second());
            scalarField& f1(table_[currentIndex_+1].second());

            //- Interpolation coefficient
            scalar c((t1-t)/(t1-t0));

            //- Linear interpolation between the fields at the provided times
            //  if the time falls in between two time bins 
            this->operator==
            (
                c*f0+(1.0-c)*f1
            );
        }
        else //- i.e. t > tLast
        {
            this->operator==(table_[table_.size()-1].second());
        }
    }
    else //- i.e. t <= tStart
    {
        this->operator==(table_[0].second());
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::timeFieldTableFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    updateField();
}


void Foam::timeFieldTableFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeEntry
    <
        List
        <
            Tuple2
            <
                scalar, 
                scalarField
            >
        >
    >("table", table_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchScalarField,
       timeFieldTableFvPatchScalarField
   );
}


// ************************************************************************* //
