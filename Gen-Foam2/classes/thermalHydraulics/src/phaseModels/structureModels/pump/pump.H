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
    Foam::pump

Description
    Allow defining a momentum source, including options for a time 
    dependent source, and a source coming from an FMU

    Authors: 
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com), 
    Carlo Fiorina for the FMU part

SourceFiles
    pump.C

\*---------------------------------------------------------------------------*/

#ifndef pump_H
#define pump_H

#include "volFields.H"
#include "meshToMesh.H"
#include "zeroGradientFvPatchFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Forward declarations
class structure;

/*---------------------------------------------------------------------------*\
                           Class pump Declaration
\*---------------------------------------------------------------------------*/

class pump
: 
    public IOdictionary
{
protected :

    //- Ref to mesh
    const fvMesh& mesh_;

    //- List of cells in which this pump exists
    const labelList& cellList_;

    //- Const value of the momentum source
    vector pumpValue_;
    
    //- Optional table to describe any possible dependence of the momentum
    //  source on time
    autoPtr<Function1<scalar>> timeProfilePtr_;

    //- Optional offset time
    scalar t0_;

    //- True if a timeProfilePtr_ exists
    bool timeProfile_;

    //- Name of the momentum source for FMU coupling
    word pumpMultiplierNameFromFMU_;

    //- True if a momentumSourceNameForFMU exists
    bool fromFMU_;

public:

    TypeName("pump");

    //- Constructors

        pump
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const labelList& cellList
        );

        //- Destructor
        virtual ~pump() {}

    //- Member Functions

        //- Update momentum source field
        void correct(volVectorField& pumpField);    

        //- Access

            //-
            bool timeDependent() const
            {
                return (timeProfile_ || fromFMU_);
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
