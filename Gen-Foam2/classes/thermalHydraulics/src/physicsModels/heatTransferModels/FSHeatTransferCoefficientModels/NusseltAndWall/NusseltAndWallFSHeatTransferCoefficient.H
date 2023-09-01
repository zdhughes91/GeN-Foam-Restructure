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
    Foam::FSHeatTransferCoefficientModels::NusseltAndWall

Description
    Heat transfer coefficient for forced convective flows that is computed
    starting from the Nusselt, which is turn obtained via correlation in the
    form

    Nu = A_+B_*(Re^C_)*(Pr^D_)

    with Re, Pr being the Reynolds (superficial in two-phase) and Prandlt
    of the fluid, and A, B, C, D selectable consts

    Additional heat transfer coeff to take into account the wall resistance,
    such that

    1/H = (1/(Nu * kappa / Dh)) + (1/H_wall)

    Example of mandatory dictionary entries:
    \verbatim
        const   4;
        coeff   0.01;
        expRe   0.8;
        expPr   0.8;
        addH    10000;
    \endverbatim

    Authors:
    Stefan Radman (sradman@pm.me / stefan.radman@epfl.ch)
    Thomas Guilbaud (thomas.guilbaud@epfl.ch)

SourceFiles
    NusseltAndWallFSHeatTransferCoefficient.C

https://foam-for-nuclear.gitlab.io/builds/foam-for-nuclear/GeN-Foam/GeN-Foam/classes/neutronics/albedoSP3/albedoSP3FvPatchField.H
\*---------------------------------------------------------------------------*/

#ifndef NusseltAndWallFSHeatTransferCoefficient_H
#define NusseltAndWallFSHeatTransferCoefficient_H

#include "FSHeatTransferCoefficientModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace FSHeatTransferCoefficientModels
{

/*---------------------------------------------------------------------------*\
                           Class Nusselt Declaration
\*---------------------------------------------------------------------------*/

class NusseltAndWall
:
    public FSHeatTransferCoefficientModel
{
private:

    //- Ref to Reynolds of pair
    const volScalarField& Re_;

    //- Fluid thermal conductivity
    const volScalarField& kappa_;

    //- Ref to Prandtl of fluid
    const volScalarField& Pr_;

    //- Ref to characteristic dimension of pair
    const volScalarField& Dh_;

    //- Coefficients in Nu = A_ + B_*Re^(C_)*Pr^(D_)
    const scalar A_;
    const scalar B_;
    const scalar C_;
    const scalar D_;

    //- Additional heat transfer for the wall resistance
    const scalar H_wall;

    //- If C_ == D_, no need to do pow(Re, C_)*pow(Pr, D_), you can just do
    //  pow(Re*Pr, C_) which is faster. Re*Pr = Pe (Peclet)
    const bool usePeclet_;

public:

    TypeName("NusseltAndWall");

    //- Constructors

        //- Construct for a fluid-structure pair
        NusseltAndWall
        (
            const FSPair& pair,
            const dictionary& dict,
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~NusseltAndWall(){}


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
