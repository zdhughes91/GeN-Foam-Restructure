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
    Foam::twoPhaseDragMultiplierModels::LockhartMartinelli

Description
    Class to handle fluid-structure friction coefficient (drag) multipliers.

    Author: 
    Stefan Radman (sradman@protonmail.com / stefan.radman@epfl.ch)

SourceFiles
    LockhartMartinelli.C
    newLockhartMartinelli.C

\*---------------------------------------------------------------------------*/

#ifndef LockhartMartinelli_H
#define LockhartMartinelli_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "twoPhaseDragMultiplierModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace twoPhaseDragMultiplierModels
{

/*---------------------------------------------------------------------------*\
                    Class LockhartMartinelli Declaration
\*---------------------------------------------------------------------------*/

class LockhartMartinelli
:
    public twoPhaseDragMultiplierModel
{
private:

    //- C in phi2 = 1 + C/X + 1/X^2
    scalar C_;

public:

    //- Runtime type information
    TypeName("LockhartMartinelli");

    // Constructors

        LockhartMartinelli
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~LockhartMartinelli(){}

    // Member Functions

        //- Return multiplier
        virtual scalar phi2(const label& celli) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace twoPhaseDragMultiplierModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
