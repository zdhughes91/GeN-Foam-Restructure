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
    Foam::linearElasticThermoMechanics

Description
    Derived class for linear elastic thermoMechanics

    Reference publications:

    Carlo Fiorina, Ivor Clifford, Manuele Aufiero, Konstantin Mikityuk,
    "GeN-Foam: a novel OpenFOAM® based multi-physics solver for 2D/3D transient
    analysis of nuclear reactors", Nuclear Engineering and Design, submitted

    Carlo Fiorina, Konstantin Mikityuk, " Application of the new GeN-Foam multi-physics
    solver to the European Sodium Fast Reactor and verification against available codes",
    Proceedings of ICAPP 2015, May 03-06, 2015 - Nice (France), Paper 15226

    Author of this file (and associated .C or included .H files): 
    Carlo Fiorina <carlo.fiorina@outlook.com; carlo.fiorina@epfl.ch;>, EPFL (Switzerland)

SourceFiles
    linearElasticThermoMechanics.C

\*---------------------------------------------------------------------------*/

#ifndef linearElasticThermoMechanics_H
#define linearElasticThermoMechanics_H

#include "thermoMechanics.H"
#include "XS.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class linearElasticThermoMechanics Declaration
\*---------------------------------------------------------------------------*/

class linearElasticThermoMechanics
:
    public thermoMechanics//class inheritance
{
protected:

    // Protected data

        //- Properties
        Switch planeStress_;
        Switch linkedFuel_;
        vector fuelOrientation_;
        dimensionedScalar TStructRef_;
        PtrList<entry> TMEntries_;
        label TMZoneNumber_;
        volScalarField rho_;
        volScalarField rhoE_;
        volScalarField nu_;
        volScalarField threeKalpha_;
        volScalarField DT_;
        volScalarField C_;
        volScalarField rhoK_;
        volScalarField alpha_;
        volScalarField E_;
        volScalarField mu_;
        volScalarField lambda_;
        volScalarField threeK_;
        volScalarField TFuelRef_;
        volScalarField alphaFuel_;
        volScalarField TCRRef_;
        volScalarField alphaCR_;

        //- Fields
        volVectorField disp_;
        volScalarField fuelDisp_;
        volScalarField CRDisp_;
        volVectorField fuelDispVector_;
        volScalarField TStruct_;
        volSymmTensorField sigmaD_;
        volVectorField divSigmaExp_;
        volScalarField TFuel_;
        volScalarField gapWidth;

        //- Ptrs to coupling fields
        const volScalarField* TFuelOrig_;
        const volScalarField* TStructOrig_;

        //- Control
        int nCorr_;
        scalar convergenceTolerance_;
        Switch compactNormalStress_;

private:

    // Private data

    // Private Member Functions

        //- Disallow default bitwise copy construct
        linearElasticThermoMechanics(const linearElasticThermoMechanics&);

        //- Disallow default bitwise assignment
        void operator=(const linearElasticThermoMechanics&);


public:

    //- Runtime type information
        TypeName("linearElasticThermoMechanics");


    // Constructors

        //- Construct from mesh, read data from IOdictionaries
        linearElasticThermoMechanics(fvMesh& mesh);

    //- Destructor
        virtual ~linearElasticThermoMechanics();


    // Member Functions

        //- Substitutes the old getFields. Set local class references to global
        //  coupling fields, to be called only once after construction. The 
        //  fields are set via lookup in the src objectRegistry (i.e. 
        //  fluidMeshNB for the current state of GeN-Foam)
        virtual void getCouplingFieldRefs
        (
            const objectRegistry& src,
            const meshToMesh& mechToFluid
        );

        //- The old getFields, yet operates on the local refs without having
        //  the need for all the fields to be passed. This effectively bypasses
        //  the old limitation of getFields, i.e. if someone wanted to
        //  implement a new thermoMechanical model that requires more fields 
        //  for the coupling, then they'd need to modify the getFields function
        //  in ALL of the other thermoMechanical models too
        virtual void interpolateCouplingFields
        (
            const meshToMesh& mechToFluid
        );

        //- Correct/update the properties
        virtual void correct(scalar& residual);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
