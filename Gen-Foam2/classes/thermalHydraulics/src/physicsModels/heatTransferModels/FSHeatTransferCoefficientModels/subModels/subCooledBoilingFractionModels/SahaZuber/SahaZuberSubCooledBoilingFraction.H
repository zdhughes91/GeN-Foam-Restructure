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
    Foam::subCooledBoilingFractionModels::SahaZuber

Description
    Calculate the temperature of onset of nucleate boiling according to the
    model by N. SahaZuber, G.R. Warrier, and V.K. Dhir, "Onset of Nucleate Boiling 
    and Active Nucleation Site Density During Subcooled Flow Boiling," 
    J. Heat Transfer, 124, 717-728, 2002..
    This model should only to be used in the following parameter range:
    *   p =         0.1 - 13.8 MPa
    *   massFlux =  95 - 2760 kg/m2/s   (i.e. ~ 0.1 - 2.7 m/s for water at ~ 
                                        1000 kg/m3)
    *   heatFlux =  0.28 - 1.89 MW/m^2

    Author: 
    Stefan Radman (sradman@protonmail.com / stefan.radman@epfl.ch)

SourceFiles
    SahaZubersubCooledBoilingFraction.C

\*---------------------------------------------------------------------------*/

#ifndef SahaZubersubCooledBoilingFraction_H
#define SahaZubersubCooledBoilingFraction_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "subCooledBoilingFractionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace subCooledBoilingFractionModels
{

/*---------------------------------------------------------------------------*\
            Class SahaZuber Declaration
\*---------------------------------------------------------------------------*/

class SahaZuber
:
    public subCooledBoilingFractionModel
{
protected:

public:

    //- Runtime type information
    TypeName("SahaZuber");

    // Constructors

        SahaZuber
        (
            const FSPair& pair,
            const dictionary& dict,
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~SahaZuber(){}

    // Member Functions

        //- Return subCooledBoilingFraction
        virtual scalar value
        (
            const label& celli,
            const scalar& qi
        ) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace subCooledBoilingFractionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
