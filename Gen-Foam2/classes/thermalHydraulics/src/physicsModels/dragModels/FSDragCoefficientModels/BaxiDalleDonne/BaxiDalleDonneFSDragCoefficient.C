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

#include "FSPair.H"
#include "BaxiDalleDonneFSDragCoefficient.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace FSDragCoefficientModels
{
    defineTypeNameAndDebug(BaxiDalleDonne, 0);
    addToRunTimeSelectionTable
    (
        FSDragCoefficientModel, 
        BaxiDalleDonne, 
        FSDragCoefficientModels
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FSDragCoefficientModels::BaxiDalleDonne::BaxiDalleDonne
(
    const FSPair& pair,
    const dictionary& dict,
    const objectRegistry& objReg
)
:
    FSDragCoefficientModel
    (
        pair,
        dict,
        objReg
    ),
    A_(0),
    B_(0),
    C_(0)
{
    scalar Dp(dict.get<scalar>("pinDiameter"));
    scalar Dw(dict.get<scalar>("wireDiameter"));
    scalar H(dict.get<scalar>("wireLeadLen"));
    scalar Pt(Dp+1.0444*Dw);
    A_ = (80.0/sqrt(100*H))*pow(Pt/Dp, 1.5); //- H provided in m but needs to
                                             //  be in cm
    B_ = 1.034/pow(Pt/Dp, 0.124);
    C_ = 29.7*pow(Pt/Dp, 6.9)/pow(H/(Dp+Dw), 2.239);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::FSDragCoefficientModels::BaxiDalleDonne::value
(
    const label& celli
) const
{
    const scalar& Rei(Re(celli));
    scalar fl(A_/Rei);  //- There should be a *pow(Twall/Tbulk, 1.5), but I am
                        //  lazy now (plus, what about cases where I do not 
                        //  solve for energy? I could handle this with some
                        //  flags and ifs)
    scalar ft(0.316*(B_+C_*pow(Rei,0.086))/pow(Rei, 0.25));
    if (Rei > 400)
    {
        if (Rei < 5000)
        {
            scalar psi
            (
                min
                (
                    max
                    (
                        ((Rei-400.0)/4600.0), 
                        0.0
                    ), 
                    1.0
                )
            );
            return sqrt(psi)*ft+sqrt(1.0-psi)*fl;
        }
        else 
            return ft;
    }
    return fl;
}

// ************************************************************************* //
