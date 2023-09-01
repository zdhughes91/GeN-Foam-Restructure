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
    Foam::FSDragCoefficientModels::NoKazimi

Description
    NoKazimi correlation for fluid-structure pressure drop
    
    https://dspace.mit.edu/handle/1721.1/60581

    Author:
    Stefan Radman (sradman@protonmail.com / stefan.radman@epfl.ch)

    Example of mandatory dictionary entries:
    \verbatim
    pinPitch        0.006;
    pinDiameter     0.0076;
    wireDiameter    0.0016;
    wireLeadLen     0.2;
    \endverbatim

SourceFiles
    NoKazimi.C

\*---------------------------------------------------------------------------*/

#ifndef NoKazimiFSDragCoefficient_H
#define NoKazimiFSDragCoefficient_H

#include "FSDragCoefficientModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Forward declarations
class FSPair;

namespace FSDragCoefficientModels
{

/*---------------------------------------------------------------------------*\
                           Class NoKazimi Declaration
\*---------------------------------------------------------------------------*/

class NoKazimi
:
    public FSDragCoefficientModel
{
protected:

    //- Pre-computed coeffs
    scalar A_;
    scalar B_;
    scalar C_;

public:

    TypeName("NoKazimi");

    //- Constructors

        //- Construct for a fluid-structure pair
        NoKazimi
        (
            const FSPair& pair,
            const dictionary& dict,
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~NoKazimi(){}

    //- Member Functions

        //- Drag coefficient
        virtual scalar value(const label& celli) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FSDragCoefficientModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //