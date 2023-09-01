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
    Foam::FSDragCoefficientModels::modifiedEngel

Description
    modifiedEngel correlation for pressure drop

    @article{BUBELIS20083299,
    title = "Review and proposal for best fit of wire-wrapped fuel bundle 
        friction factor and pressure drop predictions using various existing 
        correlations",
    journal = "Nuclear Engineering and Design",
    volume = "238",
    number = "12",
    pages = "3299 - 3320",
    year = "2008",
    issn = "0029-5493",
    doi = "https://doi.org/10.1016/j.nucengdes.2008.06.024",
    url = "http://www.sciencedirect.com/science/article/pii/S0029549308003762",
    author = "E. Bubelis and M. Schikorr",
    abstract = "The aim of this paper is to provide an overview of the existing
        wire-wrapped fuel bundle friction factor/pressure drop correlations and
        to qualitatively evaluate which of the existing friction factor 
        correlations are the best in retracing the results of a large set of 
        the experimental data available on wire-wrapped fuel assemblies tested 
        under different coolant conditions."
    }

    Author: 
    Stefan Radman (sradman@protonmail.com / stefan.radman@epfl.ch)

SourceFiles
    modifiedEngel.C

\*---------------------------------------------------------------------------*/

#ifndef modifiedEngelFSDragCoefficient_H
#define modifiedEngelFSDragCoefficient_H

#include "FSDragCoefficientModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Forward declarations
class FSPair;

namespace FSDragCoefficientModels
{

/*---------------------------------------------------------------------------*\
                           Class modifiedEngel Declaration
\*---------------------------------------------------------------------------*/

class modifiedEngel
:
    public FSDragCoefficientModel
{
public:

    TypeName("modifiedEngel");

    //- Constructors

        //- Construct for a fluid-structure pair
        modifiedEngel
        (
            const FSPair& pair,
            const dictionary& dict,
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~modifiedEngel(){}

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