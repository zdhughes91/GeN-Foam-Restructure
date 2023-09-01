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

#include "powerOffCriterionFieldValue.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace powerOffCriterionModels
{
    defineTypeNameAndDebug(fieldValue, 0);
    addToRunTimeSelectionTable
    (
        powerOffCriterionModel, 
        fieldValue, 
        powerOffCriterionModels
    );
}
}

const Foam::Enum
<
    Foam::powerOffCriterionModels::fieldValue::fieldOp
>
Foam::powerOffCriterionModels::fieldValue::fieldOpNames_
(
    {
        { 
            fieldOp::max,
            "max" 
        },
        { 
            fieldOp::min,
            "min" 
        }
    }
);

const Foam::Enum
<
    Foam::powerOffCriterionModels::fieldValue::criterion
>
Foam::powerOffCriterionModels::fieldValue::criterionNames_
(
    {
        { 
            criterion::above, 
            "valueAboveThreshold" 
        },
        { 
            criterion::below, 
            "valueBelowThreshold" 
        }
    }
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::powerOffCriterionModels::fieldValue::fieldValue
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    powerOffCriterionModel
    (
        mesh,
        dict
    ),
    fieldName_(this->get<word>("fieldName")),
    fieldOp_
    (
        fieldOpNames_.get
        (
            this->get<word>
            (
                "fieldOperation"
            )
        )
    ),
    criterion_
    (
        criterionNames_.get
        (
            this->get<word>
            (
                "criterion"
            )
        )
    ),
    threshold_(this->get<scalar>("threshold")),
    timeDelay_(this->lookupOrDefault<scalar>("timeDelay", 0.0)),
    time0_(0.0),
    fieldPtr_(nullptr),
    InfoFlag_(false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::powerOffCriterionModels::fieldValue::~fieldValue()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::powerOffCriterionModels::fieldValue::powerOffCriterion()
{   
    if (fieldPtr_ == nullptr)
    {
        fieldPtr_ = &mesh_.lookupObject<volScalarField>(fieldName_);
    }

    scalar value(0.0);
    bool flag(false);

    switch(fieldOp_)
    {
        case fieldOp::max :
        {
            value = Foam::max(*fieldPtr_).value();
            break;
        }
        case fieldOp::min :
        {
            value = Foam::min(*fieldPtr_).value();
            break;
        }
    }

    switch(criterion_)
    {
        case criterion::above :
        {
            flag = (value >= threshold_);
            break;
        }
        case criterion::below :
        {
            flag = (value <= threshold_);
            break;
        }
    }

    if (flag)
    {
        if (timeDelay_ != 0.0)
        {
            scalar time(mesh_.time().timeOutputValue());
            if (time0_ == 0.0) time0_ = time;
            flag = (time >= (time0_ + timeDelay_));
        }
        else if (!InfoFlag_) time0_ = mesh_.time().timeOutputValue();
    }
    if (flag and !InfoFlag_) InfoFlag_ = true;
    if (InfoFlag_) 
            Info << "Power off at t = " << time0_+timeDelay_ << " s" << endl;
    return flag;
}


// ************************************************************************* //
