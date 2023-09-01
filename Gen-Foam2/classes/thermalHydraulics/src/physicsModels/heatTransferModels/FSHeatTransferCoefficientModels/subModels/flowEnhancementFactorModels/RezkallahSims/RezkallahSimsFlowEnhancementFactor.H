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
    Foam::flowEnhancementFactorModels::RezkallahSims

Description
    Two-phase flow enhancement factor for the convective heat transfer
    coefficient based on the work of K.S. Rezkallah and G.E. Sims, 
    "An Examination of Correlations of Mean Heat Transfer Coefficients in 
    Two-Phase Two-Component Flow in Vertical Tubes," AIChE Symposium
    Series, 83 (257), 109-114, 1987.
    The TRACE code typically uses this over the Chen model due to greater
    numerical stability of this approach, particularly for low-pressure,
    low mass flow scenarios (see 
    https://www.nrc.gov/docs/ML1200/ML120060218.pdf, pp. 263-267)

    Author: 
    Stefan Radman (sradman@protonmail.com / stefan.radman@epfl.ch)

SourceFiles
    RezkallahSimsFlowEnhcementFactor.C

\*---------------------------------------------------------------------------*/

#ifndef RezkallahSimsFlowEnhancementFactor_H
#define RezkallahSimsFlowEnhancementFactor_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "flowEnhancementFactorModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace flowEnhancementFactorModels
{

/*---------------------------------------------------------------------------*\
            Class RezkallahSims Declaration
\*---------------------------------------------------------------------------*/

class RezkallahSims
:
    public flowEnhancementFactorModel
{
protected:

    //-
    const scalar exp_;

public:

    //- Runtime type information
    TypeName("RezkallahSims");

    // Constructors

        RezkallahSims
        (
            const FSPair& pair,
            const dictionary& dict,
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~RezkallahSims(){}

    // Member Functions

        //- Return fraction of structure interfacial contact in
        //  contact with the dispersed fluid
        virtual scalar value(const label& celli) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace flowEnhancementFactorModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
