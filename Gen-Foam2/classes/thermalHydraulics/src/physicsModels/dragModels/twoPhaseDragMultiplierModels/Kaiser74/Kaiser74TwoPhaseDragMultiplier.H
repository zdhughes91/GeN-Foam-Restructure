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
    Foam::twoPhaseDragMultiplierModels::Kaiser74

Description
    Model by Kaiser et al. for the multiplier, specifically for sodium
    two-phase flows:

    phi = 8.2/(X^0.55)

    X2 is the well-known Lockhart-Martinelli
    parameter defined as:

    X2 = (mu_I/mu_O)^0.2 * (flowQuality_I/flowQuality_O)^1.8 * rho_O/rho_I

    Please note that X in use in phi is X = sqrt(X2) and that phi2 = phi^2

    Kaiser, A., Peppler, W., Vöröss, L., 1974. "Type of flow, pressure drop, 
    and critical heat flux of a two-phase sodium flow". Nucl. Eng. Des. 30, 
    pp. 305–315.

    Author: 
    Stefan Radman (sradman@protonmail.com / stefan.radman@epfl.ch)

SourceFiles
    Kaiser74.C

\*---------------------------------------------------------------------------*/

#ifndef Kaiser74TwoPhaseDragMultiplier_H
#define Kaiser74TwoPhaseDragMultiplier_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "twoPhaseDragMultiplierModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace twoPhaseDragMultiplierModels
{

/*---------------------------------------------------------------------------*\
                    Class Kaiser74 Declaration
\*---------------------------------------------------------------------------*/

class Kaiser74
:
    public twoPhaseDragMultiplierModel
{
public:

    //- Runtime type information
    TypeName("Kaiser74");

    // Constructors

        Kaiser74
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~Kaiser74(){}

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
