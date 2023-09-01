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
    Foam::virtualMass

Description
    Class to handle virtual mass

    Author: 
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

SourceFiles
    virtualMass.C
    newInterfacialAreaModel.C

\*---------------------------------------------------------------------------*/

#ifndef virtualMass_H
#define virtualMass_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fluid.H"
#include "fvMatrix.H"
#include "virtualMassCoefficientModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class FFPair;

/*---------------------------------------------------------------------------*\
            Class virtualMass Declaration
\*---------------------------------------------------------------------------*/

class virtualMass
:
    public IOdictionary
{
protected:

    const fvMesh& mesh_;

    const FFPair& pair_;

    const fluid& fluid1_;

    const fluid& fluid2_;

    const volVectorField& U1_;
    
    const volVectorField& U2_;
    
    const surfaceScalarField& phi1_;
    
    const surfaceScalarField& phi2_;

    autoPtr<virtualMassCoefficientModel> VmPtr_;

    volScalarField Vm_;

    typedef HashTable<fvVectorMatrix, word, word::hash> VmForceTable;

    VmForceTable VmForces_;

public:

    //- Runtime type information
    TypeName("virtualMass");

    // Constructors

        virtualMass
        (
            const FFPair& pair,
            const dictionary& dict
        );

        //- Destructor
        virtual ~virtualMass(){}

    // Member Functions

        //-
        virtual void correct();

        //- Access

            const volScalarField& Vm() const
            {
                return Vm_;
            }

            const fvVectorMatrix& force(const word& UName) const
            {
                return VmForces_[UName];
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
