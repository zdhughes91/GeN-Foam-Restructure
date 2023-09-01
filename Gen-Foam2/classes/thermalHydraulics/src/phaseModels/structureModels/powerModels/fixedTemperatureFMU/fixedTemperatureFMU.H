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
    Foam::powerModels::fixedTemperatureFMU

Description
    Describes a structure with a temperature taken from an FMU

    Authors: 
    Carlo Fiorina (carlo.fiorina@epfl.ch / carlo.fiorina@outlook.com)

Usage
    \verbatim
        type            NusseltWallAndHfromFMU;
        HNameFromFMU    HCoupled;
        const   0.0;
        coeff   0.023;
        expRe   0.8;
        expPr   0.4;
        addH    30000;
    \endverbatim

SourceFiles
    fixedTemperatureFMU.C

\*---------------------------------------------------------------------------*/

#ifndef fixedTemperatureFMU_H
#define fixedTemperatureFMU_H

#include "powerModel.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace powerModels
{

/*---------------------------------------------------------------------------*\
                           Class fixedTemperatureFMU Declaration
\*---------------------------------------------------------------------------*/

class fixedTemperatureFMU
:
    public powerModel
{
protected :

    //- Surface temperature
    mutable volScalarField T_;

public:

    //- Runtime type information
    TypeName("fixedTemperatureFMU");

    //- Constructors

        fixedTemperatureFMU
        (
            structure& structure,
            const dictionary& dicts
        );

        //- Destructor
        virtual ~fixedTemperatureFMU();

    //- Member Functions

        //- In two-phase solvers :
        //  HT = H1*frac1*T1 + H2*frac2*T2
        //  H =  H1*frac1 + H2*frac2
        //- In a mono-phase solver :
        //  HT = H*T
        //  H = H
        virtual void correct
        (
            const volScalarField& HT,
            const volScalarField& H
        ){}

        //- Set provided field with structure surface temperature
        virtual void correctT(volScalarField& T) const;

        //- Turn off power
        virtual void powerOff();

        //- Temperature update based on given profile
        void temperatureUpdate() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace powerModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //