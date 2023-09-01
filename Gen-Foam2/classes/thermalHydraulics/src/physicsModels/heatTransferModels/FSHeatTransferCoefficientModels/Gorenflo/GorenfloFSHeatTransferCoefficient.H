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
    Foam::FSHeatTransferCoefficientModels::Gorenflo

Description
    Heat transfer coefficient for pool boiling scenarios calculated via the
    correlation devised by  D. Gorenflo, "Pool Boiling," VDI-Heat Atlas, 
    Sect. Ha, VDI-Verlag, Dusseldorf, 1993. Obtained from secondary reference: 
    J.G. Collier and J. R. Thome, Convective Boiling and Condensation, 
    3rd Edition, p. 155-158, Oxford University Press, Oxford, 1994.

SourceFiles
    GorenfloFSHeatTransferCoefficient.C

\*---------------------------------------------------------------------------*/

#ifndef GorenfloFSHeatTransferCoefficient_H
#define GorenfloFSHeatTransferCoefficient_H

#include "FSHeatTransferCoefficientModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace FSHeatTransferCoefficientModels
{

/*---------------------------------------------------------------------------*\
                           Class Gorenflo Declaration
\*---------------------------------------------------------------------------*/

class Gorenflo
:
    public FSHeatTransferCoefficientModel
{
private:

    bool useExplicitHeatFlux_;

    //- Ptr to structure wall temperature (only used if !useExplicitHeatFlux_)
    const volScalarField* Twall_;
 
    //- Ptr to saturation temperature (only used if !useExplicitHeatFlux_)
    const volScalarField* Tsat_;

    //- Ptr to structure heat flux (only used if useExplicitHeatFlux_)
    const volScalarField* q_;

    //- Ref to pressure
    const volScalarField& p_;

    //- Reference parameter
    //- Critical pressure in Pa
    scalar pCrit_;
    //- Heat transfer coefficient
    scalar h0_;
    //- Heat flux
    scalar q0_;
    //- Absolute surface roghness
    scalar R0_;
    
    //- Input paramenters
    //- Absolute surface roghness
    scalar R_;

    //- Pre-computed coeffs
    scalar A_;

public:

    TypeName("Gorenflo");

    //- Constructors

        //- Construct for a fluid-structure pair
        Gorenflo
        (
            const FSPair& pair,
            const dictionary& dict,
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~Gorenflo(){}


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