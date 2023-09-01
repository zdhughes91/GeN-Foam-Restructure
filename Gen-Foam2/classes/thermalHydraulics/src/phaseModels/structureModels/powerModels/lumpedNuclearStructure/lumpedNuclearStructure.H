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
    Foam::powerModels::lumpedNuclearStructure

Description
    Model for representing a lumped-parameter structure characterized
    by a user selectable number of nodes.
    \verbatim
    Tmax            T[0]           T[1]          T[n-1]         Tsurface      
      | --- H[0] --- | --- H[1] --- | --- H[2] --- |  --- H[n] --- |  
    \endverbatim
    Note that only T[] are state variables while Tmax and Tsurface can
    be found after matrix solution.
    Note also that one needs n+1 conductances when solvinf for n nodes
    The equation in the class can be obtained starting from a simple
    energy balance for the ith node:
    V_i*rho_i*cp_i * dT/dt = 
    V_i*Q_i + H[i](T_(i-1) - T_i) - H[i+1](T_(i) - T_(i+1)
    and dividing all terms by the volume of the structure.
    This means that the heat conductances  should be calculated as 
    the power flowing between two nodes divided by the diffence in 
    temperatures at steadystate, abd divided by the volume of the 
    structure. 
    Zero-gradient BCs are used at the inner boundary and convective
    BCs are used at the external boundary (in contact with the coolant)

Usage 
    The following sub-dictionary should be included as a powerModel 
    in phaseProperties.

    Example for one node:
    \verbatim
    powerModel 
    {
        type                lumpedNuclearStructure;
        volumetricArea      40;
        powerDensity        5e6; // Power density smeared over the entire
                                 // structure
        nodesNumber         1; // Number of nodes
        nodeFuel            0; // Define which node temperature is used
                               // to paramtrize XS according to 
                               // nuclearDataFuelTemp
        nodeClad            0; // Define which node temperature is used
                               // to paramtrize XS according to 
                               // nuclearDataCladExp. Of course, it does 
                               // not have to represent a cladding. It
                               // could be used to paramtrize over the
                               // graphite temperature in PBRs
        heatConductances    (450000 225000); // Heat conductances from
                                             // average to max, 
                                             // and average to surface.
                                             // Note that heat conductances (W/K) 
                                             // must be divided by the volume
                                             // of the entire structure.
                                             // There must be one more conductance 
                                             // than the number of nodes 
                                             // (see explanation in class)
        rhoCp               (53576.92); // for each node
        volumeFractions     (1); // Frction of the volume of the structure 
                                 // occupied by each node 
        powerFractions      (1); // Fraction of total power in the strucure
                                 // that goes to each node
        T0                  900; // Initial temperature, if no found in
                                 // the time folder
    }
    \endverbatim

    Example for 3 nodes:
    \verbatim
    powerModel 
    {
        type                lumpedNuclearStructure;
        volumetricArea      40;
        powerDensity        5e6;
        nodesNumber         3;
        nodeFuel            0;
        nodeClad            1;
        heatConductances    (600000 240000 600000 1200000);
        rhoCp               (53576.92 53576.92 53576.92);
        volumeFractions     (0.75 0.125 0.125);
        powerFractions      (1 0 0);
        T0                  900;
    }
    \endverbatim

    Authors: 
    Carlo Fiorina (carlo.fiorina@outlook.com / carlo.fiorina@epfl.ch)
    
SourceFiles
    lumpedNuclearStructure.C

\*---------------------------------------------------------------------------*/

#ifndef lumpedNuclearStructure_H
#define lumpedNuclearStructure_H

#include "powerModel.H"
#include "IOFieldFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace powerModels
{

/*---------------------------------------------------------------------------*\
                           Class lumpedNuclearStructure Declaration
\*---------------------------------------------------------------------------*/

class lumpedNuclearStructure
:
    public powerModel
{
protected:

    //- Field (over the global mesh) of scalarFields (over a 1-D mesh of size
    //  subMeshSize_) representing the temperatures in each node
    IOFieldField<Field, scalar> T_;

    //- Fields representing max and surface temp
    volScalarField Tmax_;
    volScalarField Tsurface_;

    //- Fraction of power from neutronics that goes to this powerModel
    //  This take sinto account for instance cases where some of the 
    //  heating goes to the moderator
    scalarList fractionOfPowerFromNeutronics_;

    //- Number of nodes
    labelList nodesNumber_;

    //- Node representing fuel for couling purposes
    labelList nodeFuel_;

    //- Node representing clad (or other structure) for coupling purposes
    labelList nodeClad_;

    //- Heat conductances between nodes (including first and last 
    //  "half" resitances)
    scalarListList Hs_;

    //- Volumetric heat capacity for each node
    scalarListList rhoCp_;

    //- Fraction of volume of the structure of each node
    scalarListList volFraction_;

    //- Fraction of total power from neutronics that goes to each node
    scalarListList qFraction_;

    //- List that maps from cell index to region index
    labelList cellToRegion_;

    //- List that maps the region index to the region name
    wordList regionIndexToRegionName_;


public:

    //- Runtime type information
    TypeName("lumpedNuclearStructure");

    // Constructors

        lumpedNuclearStructure
        (
            structure& structure,
            const dictionary& dicts
        );

        //- Destructor
        virtual ~lumpedNuclearStructure();

    // Member Functions

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