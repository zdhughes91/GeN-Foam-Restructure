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
    Foam::FFHeatTransferCoefficientModels::NoKazimi

Description
    NoKazimi model for the liquid side of a liquid-vapour heat transfer
    coefficient. Can be only used if a phaseChangeModel is specified as this
    models relies on latent heat.

    Author: 
    Stefan Radman (sradman@protonmail.com / stefan.radman@epfl.ch)

SourceFiles
    NoKazimi.C

\*---------------------------------------------------------------------------*/

#ifndef NoKazimiFFHeatTransferCoefficient_H
#define NoKazimiFFHeatTransferCoefficient_H

#include "FFHeatTransferCoefficientModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Forward declarations
class FFPair;

namespace FFHeatTransferCoefficientModels
{

/*---------------------------------------------------------------------------*\
                           Class NoKazimi Declaration
\*---------------------------------------------------------------------------*/

class NoKazimi
:
    public FFHeatTransferCoefficientModel
{
private:

    const volScalarField& alpha_;
    
    const volScalarField& rho_;

    const volScalarField& kappa_;

    const volScalarField& Cp_;

    const volScalarField& Pr_;
    
    const volScalarField& magU_;

    //- Ref to pressure
    const volScalarField& p_;

    //- Ref to fluid-fluid interfacial area
    const volScalarField& iA_;

    //- Ref to fluid-fluid interfacial temperature
    const volScalarField& iT_;

    //- A pre-computed coefficients
    scalar coeff_;

public:

    TypeName("NoKazimi");

    //- Constructors

        //- Construct for a fluid-structure pair
        NoKazimi
        (
            const FFPair& pair,
            const dictionary& dict,
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~NoKazimi(){}


    //- Member Functions

        //- Drag coefficient
        virtual scalar value(const label& celli) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FFHeatTransferCoefficientModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //