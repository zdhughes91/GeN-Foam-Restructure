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
    Foam::powerModels::fixedPower

Description
    Describes a generic structure with a constant or time depdendent internal 
    power production.
    Solves an energy equation for the structure, updating the average
    structure temperature (assumed equal to the structure surface temperature
    for fluid heat transfer purposes).

    Authors: 
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com),
    Carlo Fiorina (carlo.fiorina@epfl.ch / carlo.fiorina@outlook.com)

SourceFiles
    fixedPower.C

\*---------------------------------------------------------------------------*/

#ifndef fixedPower_H
#define fixedPower_H

#include "powerModel.H"
#include "Function1.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Forward declaration
class structure;

namespace powerModels
{

/*---------------------------------------------------------------------------*\
                           Class fixedPower Declaration
\*---------------------------------------------------------------------------*/

class fixedPower
:
    public powerModel
{
protected :

    //- Volumetric power
    volScalarField powerDensity_;

    //- Optional table to describe any possible dependence of the power
    //  on time
    PtrList<Function1<scalar>> timeProfile_;

    //- Optional offset time
    List<scalar> t0_;

    //- True if a timeProfilePtr_ exists
    List<bool> timeDependent_;

    //- Structure temperature
    volScalarField T_;

    //- Volumetric heat capacity of the structure
    volScalarField alphaRhoCp_;

public:

    //- Runtime type information
    TypeName("fixedPower");

    //- Constructors

        fixedPower
        (
            structure& structure,
            const dictionary& dicts
        );

        //- Destructor
        virtual ~fixedPower();

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
        );

        //- Set provided field with structure surface temperature
        virtual void correctT(volScalarField& T) const;

        //- Power off
        virtual void powerOff();

        //- Power update based on given profile
        void powerUpdate();

        //- Access

            //-
            //bool timeDependent() const
            //{
            //    return timeDependent_;
            //}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace powerModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //