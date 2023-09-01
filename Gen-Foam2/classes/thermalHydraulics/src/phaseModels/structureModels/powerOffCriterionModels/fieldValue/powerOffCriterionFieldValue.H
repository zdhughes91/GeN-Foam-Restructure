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
    Foam::powerOffCriterionModels::fieldValue

Description
    Turns power off if the field value maximum exceeds a provided value.

    Author: 
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

Usage 
    ("" is the placeholder for the variable to be inserted)
    \verbatim
    powerOffControlModel
    {
        type        "fieldValue";
        value       "";   
    }
    \endverbatim

SourceFiles
    fieldValue.C

\*---------------------------------------------------------------------------*/

#ifndef fieldValue_H
#define fieldValue_H

#include "powerOffCriterionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace powerOffCriterionModels
{

/*---------------------------------------------------------------------------*\
                           Class fieldValue Declaration
\*---------------------------------------------------------------------------*/

class fieldValue
:
    public powerOffCriterionModel
{
public:

    enum fieldOp
    {
        max,
        min
    };

    static const Enum<fieldOp> fieldOpNames_;

    enum criterion
    {
        above,
        below
    };

    static const Enum<criterion> criterionNames_;

private:

    const word fieldName_;

    const fieldOp fieldOp_;

    const criterion criterion_;

    const scalar threshold_;

    const scalar timeDelay_;

    scalar time0_;

    const volScalarField* fieldPtr_;

    bool InfoFlag_;

public:

    //- Runtime type information
    TypeName("fieldValue");

    //- Constructors

        fieldValue
        (
            const fvMesh& mesh,
            const dictionary& dict
        );

    //- Destructor
    virtual ~fieldValue();

    //- Member Functions

        //- Interfacial area
        virtual bool powerOffCriterion();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace powerOffControlModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //