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
    Foam::fluidDiameterModels::WallisFilm

Description
    Correlation to be used in annular flow scenarios in which the phase of
    interest to which this model is applied flows in a film along the 
    structures (e.g. pins) of interest. It assumes that the dimensionless
    film thickness d+ (=d/Dh with Dh being the structure hydraulic diameter)
    is proportional to the Wallis interfacial drag coefficient so that:
    
    0.005*75*alphaN = 0.005*300*d/Dh

    with alphaN being the normalized phase fraction of the phase of interest 
    and the film thickness is thus computed as:

    d = 0.25*alphaN*Dh

    See TRACE theory manual (https://www.nrc.gov/docs/ML1200/ML120060218.pdf,
    pp. 135-136, eqns. 4-37, 4-38) for further info  

    Author:
    Stefan Radman (stefan.radman@alumni.epfl.ch; sradman@protonmail.com)

SourceFiles
    WallisFilmFluidDiameter.C

\*---------------------------------------------------------------------------*/

#ifndef WallisFilm_H
#define WallisFilm_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fluidDiameterModel.H"

namespace Foam
{
namespace fluidDiameterModels
{

/*---------------------------------------------------------------------------*\
            Class WallisFilm Declaration
\*---------------------------------------------------------------------------*/

class WallisFilm
:
    public fluidDiameterModel
{
protected:

    //- Minimum allowable value of the phase fraction to avoid null fluid 
    //  diameter
    scalar residualAlpha_;
    
public:

    //- Runtime type information
    TypeName("WallisFilm");

    // Constructors

        WallisFilm
        (
            const FSPair& pair,
            const dictionary& dict,
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~WallisFilm(){}

    // Member Functions

        //- Return fraction of structure interfacial contact in
        //  contact with the dispersed fluid
        virtual scalar value(const label& celli) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidDiameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
