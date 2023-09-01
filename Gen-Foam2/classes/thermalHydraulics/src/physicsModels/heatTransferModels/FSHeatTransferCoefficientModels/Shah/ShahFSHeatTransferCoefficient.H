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
    Foam::FSHeatTransferCoefficientModels::Shah

Description
    Heat transfer coefficient for pool boiling scenarios calculated via the
    correlation devised by Shah, primarly for liquid metals such as sodium or
    Na-K alloys

SourceFiles
    ShahFSHeatTransferCoefficient.C

\*---------------------------------------------------------------------------*/

#ifndef ShahFSHeatTransferCoefficient_H
#define ShahFSHeatTransferCoefficient_H

#include "FSHeatTransferCoefficientModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace FSHeatTransferCoefficientModels
{

/*---------------------------------------------------------------------------*\
                           Class Shah Declaration
\*---------------------------------------------------------------------------*/

class Shah
:
    public FSHeatTransferCoefficientModel
{
private:

    bool useExplicitHeatFlux_;

    //- Ptr to structure wall temperature (only used if !useExplicitHeatFlux_)
    const volScalarField* Twall_;
 
    //- Ptr to fluid temperature (only used if !useExplicitHeatFlux_)
    const volScalarField* Tf_;

    //- Ptr to structure heat flux (only used if useExplicitHeatFlux_)
    const volScalarField* q_;

    //- Ref to pressure
    const volScalarField& p_;

    //- Critical pressure in Pa
    scalar pCrit_;

    //- Model coefficients
    scalar C0_;
    scalar C1_;
    scalar m0_;
    scalar m1_;
    scalar n_;
    scalar exp_;

    //- Interpolation limits (I've chones them arbitrarily but and abrupt
    //  change in the model coefficients is for sure not helpful). The 
    //  abrupt threhsold provided by Shah is at a reduced pressure (i.e.
    //  p_/pCrit_) of 1e-3
    scalar pR0_;
    scalar deltaPR_;

public:

    TypeName("Shah");

    //- Constructors

        //- Construct for a fluid-structure pair
        Shah
        (
            const FSPair& pair,
            const dictionary& dict,
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~Shah(){}


    //- Member Functions

        //- Drag coefficient
        virtual scalar value(const label& celli) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FSHeatTransferCoefficientModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //