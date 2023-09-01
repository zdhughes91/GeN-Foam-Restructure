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
    Foam::FFDragFactor

Description
    Base class to represent a FFDragFactor with basing name access
    functionality.

    Author: 
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

SourceFiles
    FFDragFactor.C

\*---------------------------------------------------------------------------*/

#ifndef FFDragFactor_H
#define FFDragFactor_H

#include "volFields.H"
#include "zeroGradientFvPatchFields.H"
#include "HashPtrTable.H"
#include "FFDragCoefficientModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class FFPair;

/*---------------------------------------------------------------------------*\
                           Class FFDragFactor Declaration
\*---------------------------------------------------------------------------*/

class FFDragFactor
    : public IOdictionary
{
protected :

    //- Ref to FFPair
    const FFPair& pair_;

    //- Registry where FFDragCoefficients are stored to avoid confliciting
    //  FFDragCoefficient IOobject names if storing them all in the mesh
    //  registry (is it even a problem? Either way, it's cleaner)
    objectRegistry coeffReg_;

    //- Continutous phase density
    const volScalarField& rhoc_;

    //- Characteristic dimension of dispersed phase
    const volScalarField& Dhd_;

    //- Slip velocity
    const volScalarField& magUr_;
    
    //- Drag coefficient for an isotropiuc drag factor
    autoPtr<FFDragCoefficientModel> fdPtr_;

public:

    TypeName("FFDragFactor");

    //- Constructors

        FFDragFactor
        (
            const FFPair& pair,
            const dictionary& dict
        );

        //- Destructor
        virtual ~FFDragFactor() {}

    //- Member Functions

        //- Drag coefficient used in momentum equations:
        //    ddt(alpha_i*rho_i*U_i) + ... = ... -Kd&U_i
        //  for fluid-structure drag.
        void correctField(volScalarField& Kd) const;

        //- Access

            /*
            const labelList& cells() const
            {
                return cells_;
            }*/
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
