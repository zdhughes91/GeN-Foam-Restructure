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
    Foam::powerModels::nuclearFuelPinTest

Description
    Same as nuclearFuelPin. This is a dummy class used in the ESFR
    regression test to verify capability to handle multiple
    "nuclear" power models

    Authors: 
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com), 
    Carlo Fiorina (carlo.fiorina@outlook.com / carlo.fiorina@epfl.ch)
    
SourceFiles
    nuclearFuelPinTest.C

\*---------------------------------------------------------------------------*/

#ifndef nuclearFuelPinTest_H
#define nuclearFuelPinTest_H

#include "powerModel.H"
#include "IOFieldFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace powerModels
{

/*---------------------------------------------------------------------------*\
                           Class nuclearFuelPinTest Declaration
\*---------------------------------------------------------------------------*/

class nuclearFuelPinTest
:
    public powerModel
{
protected:

    //- Field (over the global mesh) of scalarFields (over a 1-D mesh of size
    //  subMeshSize_) representing the radial fuel temperature profile
    //  across fuel and cladding
    IOFieldField<Field, scalar> Trad_;

    //- Power density of the fuel
    //volScalarField powerDensity_;

    //- Fields representing pin centerline, outer fuel, inner cladding and 
    //  outer cladding temperatures
    volScalarField Tfi_;
    volScalarField Tfo_;
    volScalarField Tci_;
    volScalarField Tco_;

    //- Average fuel and cladding temperatures
    volScalarField Tfav_;
    volScalarField Tcav_;

    //- Scalars that are input in the IOdictionary for passing min/max fuel
    //  clad temperatures
    scalar Tfmax_;
    scalar Tfmin_;
    scalar Tcmax_;
    scalar Tcmin_;

    //- Fraction of power from neutronics that goes to this powerModel
    //  This take sinto account for instance cases where some of the 
    //  heating goes to the moderator
    scalarList fractionOfPowerFromNeutronics_;

    //- Number of cells in the 1-D radial fuel mesh
    labelList fuelMeshSize_;

    //- Number of cells in the 1-D radial clad mesh
    labelList cladMeshSize_;

    //- Sum of fuelMeshSize and cladMeshSize
    labelList meshSize_;

    //- Radial position of the mesh nodes
    scalarListList r_;

    //- Fuel inner radius
    scalarList rfi_;

    //- Fuel outer radius
    scalarList rfo_;

    //- Cladding inner radius
    scalarList rci_;

    //- Cladding outer radius
    scalarList rco_;

    //- Delta radius in fuel (i.e. mesh cell length)
    scalarList drf_;

    //- Delta radius in cladding (i.e. mesh cell length)
    scalarList drc_;

    //- Gap width
    scalarList drg_;

    //- Fuel volumetric heat capacity
    scalarList rhoCpf_;

    //- Cladding volumetric heat capacity
    scalarList rhoCpc_;

    //- Fuel thermal conductivity
    scalarList kf_;

    //- Cladding thermal conductivity
    scalarList kc_;

    //- Gap heat transfer coefficient
    scalarList gapH_;

    //- True if fuel has hole
    boolList hollowFuel_;

    //- List that maps from cell index to region index
    labelList cellToRegion_;

    //- List that maps the region index to the region name
    wordList regionIndexToRegionName_;

    //-
    wordList regionNames_;

    //- Just to avoid typing constant::mathematical::pi every time...
    scalar pi_;

    //- Mesh ring areas
    scalarListList dA_;

    //- Optional table to describe the dependence of the gap condunctance on
    //  the volumetric power density
    HashTable<autoPtr<Function1<scalar>>, word, word::hash> 
        gapHPowerDensityTable_;

    //-
    boolList useGapHPowerDensityTable_;

public:

    //- Runtime type information
    TypeName("nuclearFuelPinTest");

    // Constructors

        nuclearFuelPinTest
        (
            structure& structure,
            const dictionary& dicts
        );

        //- Destructor
        virtual ~nuclearFuelPinTest();

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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace thermoSubstructures
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //