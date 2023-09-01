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
    Foam::phaseBase

Description
    Base abstract class for the fluid and structure classes.

    Author: 
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

SourceFiles
    phaseBase.C

\*---------------------------------------------------------------------------*/


#ifndef phaseBase_H
#define phaseBase_H

#include "dictionary.H"
#include "dimensionedScalar.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class phaseBase Declaration
\*---------------------------------------------------------------------------*/

class phaseBase
:
    public volScalarField
{
protected:

    //- Phase dictionary
    const dictionary dict_;

    //- Reference to mesh
    const fvMesh& mesh_;

    //- Name of phase
    const word name_;

    //- Return the residual phase-fraction for given phase
    //  Used to stabilize the phase momentum as the phase-fraction -> 0
    const dimensionedScalar residualAlpha_;
    

public:

    // Constructors

        phaseBase
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const word& name,
            const word& defaultBoundaryType,
            bool readIfPresFlag = true,
            bool writeFlag = true
        );

        //- Destructor
        virtual ~phaseBase();

    // Member Functions

        // Access

            const fvMesh& mesh() const
            {
                return mesh_;
            }
            const word& name() const
            {
                return name_;
            }
            const dictionary& dict() const
            {
                return dict_;
            }
            const dimensionedScalar& residualAlpha() const
            {
                return residualAlpha_;
            } 
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
