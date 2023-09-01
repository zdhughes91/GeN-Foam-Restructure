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
    Foam::powerModels::heatedPin

Description
    Model for representing a heated pin with constant material properties that
    is coupled to the fluid(s) via a convective boundary condition. The
    equation is solved via the finite volume method

    Authors: 
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com), 
    Carlo Fiorina (carlo.fiorina@outlook.com / carlo.fiorina@epfl.ch)
    
SourceFiles
    heatedPin.C

\*---------------------------------------------------------------------------*/

#ifndef heatedPin_H
#define heatedPin_H

#include "powerModel.H"
#include "IOFieldFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace powerModels
{

/*---------------------------------------------------------------------------*\
                           Class heatedPin Declaration
\*---------------------------------------------------------------------------*/

class heatedPin
:
    public powerModel
{
protected:

    //- Field (over the global mesh) of scalarFields (over a 1-D mesh of size
    //  subMeshSize_) representing the radial fuel temperature profile
    //  across fuel and cladding
    IOFieldField<Field, scalar> Trad_;

    //- Power density of the pin
    volScalarField powerDensity_;

    //- Fields representing inner and outer pin temperatures 
    volScalarField Ti_;
    volScalarField To_;

    //- Average temperature
    volScalarField Tav_;

    //- Scalars that are input in the IOdictionary for passing min/max
    //  temperatures
    scalar Tmax_;
    scalar Tmin_;

    //- Number of cells in the 1-D radial mesh
    labelList meshSize_;

    //- Radial position of the mesh nodes
    scalarListList r_;

    //- Inner pin radius (can be 0)
    scalarList ri_;

    //- Outer pin radius
    scalarList ro_;

    //- Delta radius
    scalarList dr_;

    //- Pin volumetric heat capacity
    scalarList rhoCp_;

    //- Pin thermal conductivity
    scalarList k_;

    //- True if pin has hole
    boolList hollow_;

    //- List that maps from cell index to region index
    labelList cellToRegion_;

    //- Just to avoid typing constant::mathematical::pi every time...
    scalar pi_;

    //- Mesh ring areas
    scalarListList dA_;

public:

    //- Runtime type information
    TypeName("heatedPin");

    // Constructors

        heatedPin
        (
            structure& structure,
            const dictionary& dicts
        );

        //- Destructor
        virtual ~heatedPin();

    // Member Functions

        //- Set interfacial area consistently with the pin geometry. In
        //  particular, if the overall total volume-fraction of the pins
        //  in a given mesh cell is alpha, it can be shown that the interfacial
        //  area density iA of the pin surface is 2*alpha/outerRadius
        virtual void setInterfacialArea();

        void updateLocalAvgGlobalMinMaxT
        (
            const label& starti,
            const label& endi,
            const scalarList& r,
            const scalar& dr,
            const scalarField& Trad,
            scalar& Tav,
            scalar& Tmin,
            scalar& Tmax
        );

        //- Update temperature profiles Tc, Tf, in a single domain cell with
        //  current system conditions
        void updateLocalTemperatureProfile
        (
            const label& celli, 
            const scalar& ht, 
            const scalar& h
        );

        //- In two-phase solvers :
        //  HT = H1*frac1*Tfluid1 + H2*frac2*Tfluid2
        //  H =  H1*frac1 + H2*frac2
        //- In a mono-phase solver :
        //  HT = H*Tfluid
        //  H = H
        virtual void correct
        (
            const volScalarField& HT,
            const volScalarField& H
        );

        //- Return surface temperature (i.e. outer cladding temperature)
        virtual void correctT(volScalarField& T) const;

        //- Turn off pin power
        virtual void powerOff();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace thermoSubstructures
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //