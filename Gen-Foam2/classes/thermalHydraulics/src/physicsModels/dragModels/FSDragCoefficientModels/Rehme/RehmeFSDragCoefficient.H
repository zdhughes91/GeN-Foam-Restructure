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
    Foam::FSDragCoefficientModels::Rehme

Description
    Rheme correlation for pressure drop

    @ARTICLE{Rehme197315,
    author={Rehme, Klaus},
    title={PRESSURE DROP CORRELATIONS FOR FUEL ELEMENT SPACERS.},
    journal={Nuclear Technology},
    year={1973},
    volume={17},
    number={1},
    pages={15-23},
    doi={10.13182/NT73-A31250},
    }

    Example of mandatory dictionary entries:
    \verbatim
    numberOfPins        37;
    pinDiameter         0.006000;
    pinPitch            0.007600;
    wireDiameter        0.0016;
    wireLeadLen         0.1;        //- I.e. pitch of the wrapper wire
    wetWrapPerimeter    0.175237;   //- Wetted wrapper perimeter (wrapper only)
    \endverbatim

    Author: 
    Stefan Radman (sradman@protonmail.com / stefan.radman@epfl.ch)

SourceFiles
    Rehme.C

\*---------------------------------------------------------------------------*/

#ifndef RehmeFSDragCoefficient_H
#define RehmeFSDragCoefficient_H

#include "FSDragCoefficientModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Forward declarations
class FSPair;

namespace FSDragCoefficientModels
{

/*---------------------------------------------------------------------------*\
                           Class Rehme Declaration
\*---------------------------------------------------------------------------*/

class Rehme
:
    public FSDragCoefficientModel
{
protected:

    //- Pre-computed coeffs
    scalar A_;
    scalar B1_;
    scalar B2_;

public:

    TypeName("Rehme");

    //- Constructors

        //- Construct for a fluid-structure pair
        Rehme
        (
            const FSPair& pair,
            const dictionary& dict,
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~Rehme(){}

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