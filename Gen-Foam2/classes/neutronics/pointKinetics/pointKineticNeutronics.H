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
    Foam::pointKineticNeutronics

Group
    none

Description
    Run-time selectable point-kinetic based neutronics solver

    Author of this file (and associated .C or included .H files): 
    Stefan Radman <stefanradman92@gmail.com; stefan.radman@epfl.ch;>
    Carlo Fiorina <carlo.fiorina@outlook.com; carlo.fiorina@epfl.ch;>, EPFL (Switzerland)
    Thomas Guilbaud <thomas.guilbaud@epfl.ch>, EPFL/Transmutex SA (Switzerland)

SourceFiles
    pointKineticNeutronics.C

\*---------------------------------------------------------------------------*/

#ifndef pointKineticNeutronics_H
#define pointKineticNeutronics_H

#include "neutronics.H"
#include "Function1.H"
#include "interpolationTable.H"
#include "timeProfile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class pointKineticNeutronics Declaration
\*---------------------------------------------------------------------------*/

class pointKineticNeutronics
:
    public neutronics
{
protected:

    //-
    IOdictionary nuclearData_;

    //- Total power, sum of fission power and decay power
    scalar power_;

    //- Point kinetics operates on fissionPower
    scalar fissionPower_;

    //- Decay power value set by decayPowerPtr_
    scalar decayPower_;

    //- Optional decay power
    autoPtr<Function1<scalar>> decayPowerPtr_;

    //- Time offset for the calcuation of the decayPower based on the
    //  provided Function1, if any
    scalar decayPowerStartTime_;

    //- Due to whatever the user might please
    timeProfile externalReactivityTimeProfile_;

    //- Used to set an initial steady state for liquid fuel
    scalar precEquilibriumReactivity_;

    //- Estimated value of Beta for liquid fuel
    scalar liquidFuelBeta_;

    //-
    scalar totalReactivity_;

    //-
    scalar promptGenerationTime_;

    //- Total effective delayed neutron fraction
    scalar beta_;

    //- Effective delayed fractions for each delayed group
    scalarList betas_;

    //- Decay constants for each delayed group
    scalarList lambdas_;

    //- Precursor "power" for each delayed group in W to be compatible with
    //  point kinetics being solved for power rather than neutron density
    scalarList precursorPowers_;

    //- Number of delayed groups
    label delayedGroups_;

    //- Old time power, precursorPowers to be used to calculate derivatives,
    //  as well as timeIndex to check if time step has changed (can't think
    //  of a better way to get something equivalent to the finalInnerIter flag)
    label timeIndex_;
    scalar powerOld_;
    scalar fissionPowerOld_;
    scalarList precursorPowersOld_;

    // default precursor field. Only used for liquid fuel
    volScalarField defaultPrec_;

    //- Feedback coefficients in 1/K (i.e. 1e5 pcm/K)
    scalar coeffTFuel_;
    scalar coeffTClad_;
    scalar coeffTCool_;
    scalar coeffRhoCool_;
    scalar coeffTStruct_;

    //- Absolute driveline expansion coefficient in m/K
    scalar coeffDrivelineExp_;

    //- Optional Boron reactivity curve
    timeProfile boronReactivityTimeProfile_;

    //- Reference integral values for feedbacks
    scalar TFuelRef_;
    scalar TCladRef_;
    scalar TCoolRef_;
    scalar rhoCoolRef_;
    scalar TStructRef_;
    scalar TDrivelineRef_;

    //- Ptrs to coupling fields (i.e. defined on their original mesh)
    const volScalarField* TFuelOrig_;
    const volScalarField* TCladOrig_;
    const volScalarField* TCoolOrig_;
    const volScalarField* rhoCoolOrig_;
    const volScalarField* TStructOrig_;
    const volScalarField* powerDensityOrig_;
    const volScalarField* powerDensityToLiquidOrig_;

    //- These have no use for now, will be needed by the liquidFuel version
    const volVectorField* UOrig_;
    const volScalarField* alphaOrig_;
    const volScalarField* alphatOrig_;
    const volScalarField* muOrig_;

    //- Mapped fields
    volScalarField TFuel_;
    volScalarField TClad_;
    volScalarField TCool_;
    volScalarField rhoCool_;
    volScalarField TStruct_; //- This is also used to compute TDriveline

    //- Also mapped, but these are autoPtr as these field are not needed at all
    //  unless doing liquidFuel. Thus, avoid wasting memory by only setting
    //  these fields if they are actually needed
    autoPtr<volVectorField> UPtr_;
    autoPtr<volScalarField> alphaPtr_;
    autoPtr<volScalarField> alphatPtr_;
    autoPtr<volScalarField> muPtr_;

    //- Fields derived from mapped fields
    autoPtr<surfaceScalarField> phiPtr_;
    autoPtr<volScalarField> diffCoeffPrecPtr_;

    //- Fields necessary when an albedo boundary condition has been used
    //  to calculate the initial fluxes
    volScalarField Dalbedo_;
    volScalarField fluxStarAlbedo_;

    //-
    volScalarField oneGroupFlux_;

    //- Initial flux
    volScalarField initOneGroupFlux_;

    //- Normalized initial flux
    volScalarField initOneGroupFluxN_;

    //- Integral of the square of the normalized flux
    scalar domainIntegratedInitOneGroupFluxN_;

    //- Number of energy groups (for scaling fluxes of a possibly existing
    //  previous simulation)
    label energyGroups_;

    //- Group-wise fluxes
    PtrList<volScalarField> fluxes_;

    //- Delayed-group-wise precursor concentrations
    PtrList<volScalarField> precursors_;

    //- Delayed-group-wise precursor concentrations resulting from point
    //  kinetics calculations with liquid fuel. This second set of fields is
    //  necessary because the unit for precursors in point kinetics is power,
    //   instead of 1/m3
    PtrList<volScalarField> precPK_;//reconstructed based on porosity
    PtrList<volScalarField> precPKStar_;//solved for

    //- If true, fuel temperature feedback (Doppler) computed as
    //  deltaReactivity = -coeffFastDoppler*ln(TFuel/TFuelRef) with
    //  coeffFastDoppler being the Doppler constant. If false, Doppler feedback
    //  computed as deltaReactivity = coeffTFuel*(TFuel-TFuelRef)
    bool fastNeutrons_;

    //- Doppler constant if doing fastNeutrons
    scalar coeffFastDoppler_;

    //- Fields to mark which cells are part of X and which aren't, used for
    //  computing the 0-D feedback values (TFuel, TStructures, etc...).
    //  Indeed, the spatial extent of these marker fields should reflect
    //  the spatial extent of the domain that was perturbed in neutronics
    //  calculations to compute feedback coefficients

    //- For TFuel and TCladding
    volScalarField fuelFeedbackCellField_;

    //- For rhoCoolant and TCoolant
    volScalarField coolFeedbackCellField_;

    //- For TStructures
    volScalarField structFeedbackCellField_;

    //- For the control rod driveline
    volScalarField drivelineFeedbackCellField_;

    //- Pairs of control rod insertion value (in m w.r.t. a reference) vs
    // reactivity in absolute units (i.e. 1e5 pcm)
    List<Pair<scalar>> controlRodReactivityMap_;

    //- Schmidt number. Only used for liquid fuel
    scalar ScNo_;

    //- Flag to initialize the concentration of precursors to their
    //  equilibtium value at the beginning of the simulation
    bool initPrecursorsLiquidFuel_;

    //- GEM related

    //- Reactivity map
    List<Pair<scalar>> GEMReactivityMap_;

    //- Reference height of sodium in GEM
    scalar GEMSodiumLevelRef_;

    //- Ref to volumetric flux on fluid mesh
    const surfaceScalarField* phiOrig_;

    //- Reference surface-integrated value of phi
    scalar intPhiRef_;

    //- List of face labels over which to compute integral flow values
    labelList phiFaces_;

    //- Surface areas of the phiFaces_ (I cannot access fvMesh from polyMesh,
    //  which is the only mesh-related object I can access from the mapping,
    //  thus magSf needs to be computed by calling the faceArea() method
    //  of polyMesh. To avoid that, I calculated the areas only once and
    //  store them in pgiMagSf_)
    scalarList phiMagSf_;

private:

    //- Disallow default bitwise copy construct
    pointKineticNeutronics(const pointKineticNeutronics&);

    //- Disallow default bitwise assignment
    void operator=(const pointKineticNeutronics&);

public:

    //- Runtime type information
        TypeName("pointKinetics");

    //- Constructors

        //- Construct from mesh, read data from IOdictionaries
        pointKineticNeutronics(fvMesh& mesh);

    //- Destructor
        virtual ~pointKineticNeutronics(){}

    //- Member Functions

        //- Calc keff from reactivity
        virtual scalar keff() const override
        {
            return 1.0/(1.0-totalReactivity_);
        }

        //- Update all the oneGroupFlux related quantities following changes in
        //  the oneGroupFlux after proper initialization
        void setInitOneGroupFlux();

        //- Used to set the *FeedbackCellField fields
        void setFeedbackCellField
        (
            volScalarField& feedbackCellField,
            const word& keyword
        );

        //- Compute driveline expansion reactivity from controlRod reactivity
        //  map given the absolute driveline expansion in meters
        scalar calcDrivelineReactivity(const scalar& drivelineExpansion);

        //- Compute GEM sodium height and their reactivity contributions
        Pair<scalar> calcGEMLevelAndReactivity();

        //- Used to set all the pointer member variables (the members that end
        //  with *Orig_) after class initialization (needs to be called by
        //  solver main)
        virtual void getCouplingFieldRefs
        (
            const objectRegistry& src,
            const meshToMesh& neutroToFluid
        );

        //- Update local (i.e. on local neutronicMesh) instances of the Orig_
        //  fields defined on the (baffleless) thermalHydraulic mesh
        virtual void interpolateCouplingFields
        (
            const meshToMesh& neutroToFluid
        );

        //- Solve pointKinetics
        virtual void correct
        (
            scalar& residual,
            label couplingIter
        );

        //- To split strings. The file where I got this function originally
        //  defined is in classes/thermalHydraulics/models/include. However,
        //  including things from thermalHydraulics in neutronics, even if it's
        //  just a helper function such as this one, is tremendously
        //  un-elegant. I don't want the code structure to be spaghetti, so
        //  I am re-implementing the function here

        //- The deliminter is a char, so it must be provided in single-
        //  quotation marks. Type can be either word or string
        template<class Type>
        List<Type> split(const Type& compoundItem, char d)
        {
            List<Type> items(1, "");
            int i = 0;
            for (char const &c : compoundItem)
            {
                if (c == d)
                {
                    i += 1;
                    items.append("");
                }
                else items[i] += c;
            }
            return items;
        }

        //- Access

            //- Return fluxes
            virtual const PtrList<volScalarField> fluxes() const
            {
                return fluxes_;
            };

            //- Return precursors
            virtual const PtrList<volScalarField> precursors() const
            {
                return precursors_;
            };

            //- Return one group flux
            virtual const volScalarField& oneGroupFlux() const
            {
                return oneGroupFlux_;
            }

            //- Return power
            virtual scalar power() const
            {
                return power_;
            };
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
