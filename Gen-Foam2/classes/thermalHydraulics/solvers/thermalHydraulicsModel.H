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
    Foam::thermalHydraulicsModel

Description
    Run-time selectable base abstract class to represent a thermal-hydraulic 
    solver for an indefinite number of phases.

    Author: 
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

SourceFiles
    thermalHydraulicsModel.C

\*---------------------------------------------------------------------------*/

#ifndef thermalHydraulicsModel_H
#define thermalHydraulicsModel_H

#include "fvCFD.H"
#include "fvOptions.H"
#include "IOdictionary.H"
#include "customPimpleControl.H"
#include "runTimeSelectionTables.H"

//- It is pointless to include these in the .C files as any derived solver 
//  class will need all of these anyway, so I include them in the header too
#include "volFields.H"
#include "UniformField.H"
#include "surfaceFields.H"
#include "surfaceInterpolate.H"
#include "fixedValueFvsPatchFields.H"
#include "zeroGradientFvPatchFields.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvmDdt.H"
#include "fvcCurl.H"
#include "fvcFlux.H"
#include "fvMatrix.H"
#include "fvcSnGrad.H"
#include "fvmLaplacian.H"
#include "HashPtrTable.H"

#include "fluid.H"
#include "structure.H"
#include "FFPair.H"
#include "FSPair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class regimeMapModel;

/*---------------------------------------------------------------------------*\
                      Class thermalHydraulicsModel Declaration
\*---------------------------------------------------------------------------*/

class thermalHydraulicsModel
:
    public IOdictionary
{
public :

    enum momentumMode
    {
        cellCentered, 
        cellCenteredFaceReconstruction,
        faceCentered
    };

    static const Enum<momentumMode> momentumModeNames_;

protected:

    Time& runTime_;

    fvMesh& mesh_;

    customPimpleControl& pimple_;

    fv::options& fvOptions_;

    //- Optional regime maps that can be used by individual physics models
    typedef HashPtrTable
    <
        regimeMapModel,
        word,
        word::hash
    >   regimeMapTable;
    regimeMapTable regimeMapModels_;

    //- Total volumetric flux
    surfaceScalarField phi_;

    //- Gravity related fields
    uniformDimensionedVectorField g_;
    uniformDimensionedScalarField hRef_;
    volScalarField gh_;
    surfaceScalarField ghf_;

    //- Pressure
    volScalarField p_rgh_;
    volScalarField p_;

    //- Minimum allowable pressure
    dimensionedScalar pMin_;

    //- Reference cell and pressure if no fixedValue BCs provided
    label pRefCell_;
    scalar pRefValue_;
    bool forcePRef_;

    //- Initial fluid(s) mass
    dimensionedScalar initialFluidMass_;

    //- Maximum Courant number
    scalar CoNum_;

    //- Mean Courant number
    scalar meanCoNum_;

    //- Pressure solution options
    const momentumMode momentumMode_;

    //- Feature that removes velocity oscillations and/or smearing in the 
    //  main flow direction. However, if dealing with structures that have
    //  large porosity gradiens (i.e. > 0.5 between consecutive cells) or
    //  large gradients in the structure drag tensor, this feature may
    //  actually induce velocity oscillations in the directions 
    //  perpendicular to the main flow direction.
    //  This factor weights the "traditional" faceMomentum velocity 
    //  reconstruction (the one which always causes smearing/oscillations)
    //  in the main flow direction) with a moving phase fraction-adjusted 
    //  velocity reconstruction (which may induced the aforementioned
    //  transversal oscillations). This feature only has an affect on the 
    //  cellCentered and the cellCenteredFaceReconstruction momentum modes.
    //  Deafults to 0, which amounts to a "traditional" velocity reconstruction
    //  procedure from phi
    const scalar porousInterfaceSharpness_;

    //- Ptrs to coupling fields (i.e. defined on their original mesh)
    const volScalarField* powerDensityOrig_;
    const volScalarField* secondaryPowerDenistyOrig_; // Used for graphite in MSRs

    //- Mapped fields
    volScalarField powerDensityNeutronics_;
    volScalarField powerDensityNeutronicsToLiquid_;

public:

    //- Runtime type information
    TypeName("thermalHydraulicsModel");

    //- Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            thermalHydraulicsModel,
            thermalHydraulicsModels,
            (
                Time& time,
                fvMesh& mesh,
                customPimpleControl& pimple,
                fv::options& fvOpt
            ),
            (time, mesh, pimple, fvOpt)
        );

    //- Constructors

        thermalHydraulicsModel
        (
            Time& time,
            fvMesh& mesh,
            customPimpleControl& pimple,
            fv::options& fvOpt
        );

        //- Destructor
        virtual ~thermalHydraulicsModel();

    //- Selectors

        static autoPtr<thermalHydraulicsModel> New
        (
            Time& time,
            fvMesh& mesh,
            customPimpleControl& pimple,
            fv::options& fvOpt
        );

    //- Member Functions

        //- Substitutes the old getFields. Set local class references to global
        //  coupling fields, to be called only once after construction. The 
        //  fields are set via lookup in the src objectRegistry (i.e. 
        //  neutroMesh for the current state of GeN-Foam)
        virtual void getCouplingFieldRefs
        (
            const objectRegistry& src,
            const meshToMesh& fluidToNeutro
        );

        //- The old getFields, yet operates on the local refs without having
        //  the need for all the fields to be passed. This effectively bypasses
        //  the old limitation of getFields, i.e. if someone wanted to
        //  implement a new  model that requires more fields for the
        //  coupling. Splitting getFields into a getCouplingFieldRefs and 
        //  interpolateCouplingFields allows one to be as flexible as possible
        //  with new classes implementations
        virtual void interpolateCouplingFields
        (
            const meshToMesh& fluidToNeutro,
            label liquidFuel
        );

        //- Solve thermal-hydraulics, combination of the other corrects
        virtual void correct
        (
            scalar& residual,
            bool solveFluidDynamics, 
            bool solveEnergy
        ) = 0;

        virtual void correctRegimeMaps();

        //- Correct fluid and structure temperatures without
        //  updating the flow regime map
        virtual void correctEnergy(scalar& residual) = 0;

        //- Correct alphas, pressure, velocities without updating 
        //  the flow regime map
        virtual void correctFluidMechanics(scalar& residual) = 0;

        //- Correct flow regime map and all the coupling fields
        virtual void correctModels
        (
            bool solveFluidDynamics, 
            bool solveEnergy
        ) = 0;

        //- Update and print Courant number
        virtual void correctCourant() = 0;

        //- Stop execution (for development purposes)
        void stop();

        //- A function to allow for the implementation of custom adjustTimeStep
        //  criteria
        virtual void adjustTimeStep();

        // Access

            //- Return maximum Courant number
            scalar CoNum() const
            {
                return CoNum_;
            }

            //- Return mean Courant number
            scalar meanCoNum() const
            {
                return meanCoNum_;
            }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
