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
    Foam::fluid

Description
    Class to model a generic fluid phase. It mainly acts as a variable
    placeholder, as little functionality is deferred to the class itself.

    Author: 
    Stefan Radman (stefanradman92@gmail.com, stefan.radman@epfl.ch)

SourceFiles
    fluid.C

\*---------------------------------------------------------------------------*/

#ifndef fluid_H
#define fluid_H

#include "phaseBase.H"
#include "rhoThermo.H"
#include "phaseCompressibleTurbulenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fluidDiameterModel;

/*---------------------------------------------------------------------------*\
                           Class fluid Declaration
\*---------------------------------------------------------------------------*/

class fluid
:
    public phaseBase
{
public:

    enum stateOfMatter
    {
        undetermined,
        liquid,
        gas
    };

    static const Enum<stateOfMatter> stateOfMatterNames_;

private :

    //- Thermophysical properties
    autoPtr<rhoThermo> thermo_;

    //- State of matter
    const stateOfMatter stateOfMatter_;

    //- Velocity
    volVectorField U_;

    //- Velocity magnitude;
    volScalarField magU_;

    //- Superficial flux of the phase
    autoPtr<surfaceScalarField> phiPtr_;

    //- Volumetric flux of the phase
    surfaceScalarField alphaPhi_;

    //- Mass flux of the phase
    surfaceScalarField alphaRhoPhi_;

    //- Mass flux (but non-directional and stored in cell centers, needed only
    //  for regimeMap parametrization (if any is wanted/defined)
    volScalarField alphaRhoMagU_;

    //- Dilation term of the phase.
    //  It encompasses the "compressible" part of the continuity
    //  equation, including mass sources/sinks, if any.
    //  Example given in meta-code:
    /*
        ddt(rho*alpha) + div(alphaRhoPhi) = dmdt

        -> (expand the derivatives)

        rho1*ddt(alpha) + rho1*div(alphaPhi) +
        alpha1*ddt(rho) + alphaPhi1&grad(rho) = dmdt

        -> (re-arrange)

        ddt(alpha) + div(alphaPhi) 
        =
        (1/rho)*
        (
            alpha*ddt(rho) + alphaPhil&grad(rho) + dmdt
        ) 
        =
        dgdt
    */
    volScalarField dgdt_;

    //- A field that represents the sign of phase change in each cell, i.e.
    //  1 if the phase is appearning, -1 if disappearing, 0 else. This field
    //  is currently only used by the nucleateBoiling heatTransfer model
    autoPtr<scalarField> phaseChangeSignPtr_; 

    //- Flow quality of the phase. While these are not intrinsic
    //  properties to the fluid itself (they depend on the properties of the
    //  other fluid too), their are placed here for programming convenience.
    //  Handled by an autoPtr as it is non-trivial only in two-phase scenarios
    autoPtr<volScalarField> flowQualityPtr_;

    //- Place-holder for the Lockhart-Martinelli parameter. Logically, just as
    //  flow quality, it should belong to FFPair. However, if I place it there
    //  I loose the ability to compute it differently cell-by-cell based on
    //  the fluid multiplier that was selected for the regime that exists in
    //  said cell. Conversely, the twoPhaseDragMultiplierModels access a phase
    //  directly, so by placing it here, I avoid a lot of issues. I generally
    //  favour coding simplicity over philosophical consistency.
    //  Handled by an autoPtr as it is useless for the onePhase solver
    autoPtr<volScalarField> XLMPtr_;

    //- Minimum allowable value of the Lockhart-Martinelli parameter
    scalar minXLM_;

    //- Maximum allowable value of the Lockhart-Martinelli parameter
    scalar maxXLM_;

    //- Normalized phase fraction
    volScalarField normalized_;

    //- Molecular thermal conductivity and heat capacitiy. These are computed
    //  on the go by the thermo() model, and I need a ref to these that
    //  simply returns the current values, instead of having to re-compute them
    //  every time. These fields are thus place holders that need to be updated
    //  by the user
    volScalarField kappa_;
    volScalarField Cp_;
    volScalarField mu_;

    //- Prandtl of the fluid
    volScalarField Pr_;

    //- Additional powerDensity deposited directly to the fluid
    //autoPtr<volScalarField> powerDensityPtr_;

    //- Ref to powerDensity from neutronics (projected on fluid mesh)
    volScalarField& powerDensityNeutronicsToLiquid_;

    //- Continuity error of the phase 
    //  (== fvc::ddt(alpha, rho) + fvc::div(alphaRhoPhi)) (kg/m3/s)
    volScalarField contErr_;

    //- Cumulat-ive/ed continuity error integrated in time and volume (kg)
    scalar cumulContErr_;

    //- 
    autoPtr<fluidDiameterModel> diameterPtr_;

    //- Characteristic dimension of the fluid. This is equal to the structure
    //  hydraulic diameter if this fluid is continuous, else it is specified
    //  by the fluidDiameterModel within fluidGeometry. Information on whether
    //  a fluid is continuous or dispersed in a certain regime is specified in
    //  the fluidGeometry dict as well. Since this is regime-dependent, this 
    //  field is manipulated by reference by diameterPtr_ and can be different
    //  from the structure Dh only in two-phase simulations
    volScalarField Dh_;

    //- Marker field that is 1 where the phase is dispersed, 0 where it is 
    //  continuous, and between 0 and 1 in interpolation regimes.
    //  Handled by an autoPtr as it is non-trivial only in two-phase scenarios
    autoPtr<scalarField> dispersionPtr_;

    //- Turbulence model
    autoPtr<phaseCompressibleTurbulenceModel> turbulence_;

    //- Phase fraction below which the fluid energy equation is not solved
    //  and the phase temperature is set to the interfacial temperature
    dimensionedScalar thermoResidualAlpha_;

    //- Marker fields that are either 1.0 or 0.0 depending on whether the fluid
    //  is above or below thermoResidualAlpha_
    autoPtr<volScalarField> aboveThermoResidualAlphaPtr_;
    autoPtr<volScalarField> belowThermoResidualAlphaPtr_;

    //- True if equationOfState is Boussinesq
    bool Boussinesq_;

    //- Reference density field when using Boussinesq (never used otherwise)
    volScalarField* rho0Ptr_;

public:

    // Constructors

        fluid
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const word& phaseName,
            volScalarField& powerDensityNeutronicsToLiquid,
            bool readIfPresentAndWrite = true
        );

        //- Destructor
        virtual ~fluid();

    //- Member Functions

        //- Set all the autoPtr to the fields required for two-phase 
        //  simulations
        void initTwoPhaseFields();

        //- Initialize alphaPhi and alphaRhoPhi, to be called as a separate
        //  function in the solver main after adequate phase fraction
        //  normalization
        void initAlphaPhis();

        //- Construct diameter model. Only relevant in twoPhase simulations
        void constructDiameterModel();

        //- Construct turbulence model. Required as the turbulence model might
        //  require refs to fields (via objectRegistry lookup) that do not or
        //  cannot exist before the creation of the fluid model. An example is
        //  the fluid-structure Reynolds number, which exists in FSPair, that
        //  is created only after both the structure and fluid are
        //  created. Of course, details depend on the actual run time
        //  selected model, but in principle this is the reason why it is
        //  desireable to have a separate "constructur function" for the
        //  turbulenceModel
        void constructTurbulenceModel();

        //- I recall that alphaRhoMagU serves no purpose other than being a
        //  placeholder for the mass flux if requested by any regimeMap model
        void correctAlphaRhoMagU();

        //- Update diameter based on the diameterModel
        void correctDiameter();

        //- Update/set the belowThermoResidualAlpha/aboveThermoResidualAlpha 
        //  fields. Only useful in twoPhase simulations
        void correctThermoResidualMarkers();

        //- Boussinesq-sensitive density. If both the Boussinesq_ and 
        //  isConstant flags are true, the returned density is equal
        //  to rho0 provided in thermoPhysicalProperties in equationOfState
        volScalarField& rho(bool isVariableIfBoussinesq=false);
        const volScalarField& rho(bool isVariableIfBoussinesq=false) const;

        //- Related to stateOfMatter
        bool isLiquid() const
        {
            return stateOfMatter_ == stateOfMatter::liquid;
        }
        bool isGas() const
        {
            return stateOfMatter_ == stateOfMatter::gas;
        }

        //- Access (const and non-const)

            const phaseCompressibleTurbulenceModel& turbulence() const
            {
                return turbulence_();
            }

            phaseCompressibleTurbulenceModel& turbulence()
            {
                return turbulence_();
            }

            const rhoThermo& thermo() const
            {
                return thermo_();
            }

            rhoThermo& thermo()
            {
                return thermo_();
            }

            const volScalarField& T() const
            {
                return thermo_->T();
            }

            volScalarField& T()
            {
                return thermo_->T();
            }

            const volVectorField& U() const
            {
                return U_;
            }

            volVectorField& U()
            {
                return U_;
            }

            const volScalarField& magU() const
            {
                return magU_;
            }

            volScalarField& magU()
            {
                return magU_;
            }

            const surfaceScalarField& phi() const
            {
                return phiPtr_();
            }

            surfaceScalarField& phi()
            {
                return phiPtr_();
            }

            const surfaceScalarField& alphaPhi() const
            {
                return alphaPhi_;
            }

            surfaceScalarField& alphaPhi()
            {
                return alphaPhi_;
            }

            const surfaceScalarField& alphaRhoPhi() const
            {
                return alphaRhoPhi_;
            }

            surfaceScalarField& alphaRhoPhi()
            {
                return alphaRhoPhi_;
            }

            const volScalarField& dgdt() const
            {
                return dgdt_;
            }

            volScalarField& dgdt()
            {
                return dgdt_;
            }

            const scalarField& phaseChangeSign() const
            {
                return phaseChangeSignPtr_();
            }

            scalarField& phaseChangeSign()
            {
                return phaseChangeSignPtr_();
            }

            const volScalarField& flowQuality() const
            {
                return flowQualityPtr_();
            }

            volScalarField& flowQuality()
            {
                return flowQualityPtr_();
            }

            const volScalarField& XLM() const
            {
                return XLMPtr_();
            }

            volScalarField& XLM()
            {
                return XLMPtr_();
            }

            const scalar& minXLM() const
            {
                return minXLM_;
            }

            const scalar& maxXLM() const
            {
                return maxXLM_;
            }

            const volScalarField& normalized() const
            {
                return normalized_;
            }

            volScalarField& normalized()
            {
                return normalized_;
            }

            const volScalarField& kappa() const
            {
                return kappa_;
            }

            volScalarField& kappa()
            {
                return kappa_;
            }

            const volScalarField& Cp() const
            {
                return Cp_;
            }

            volScalarField& Cp()
            {
                return Cp_;
            }

            const volScalarField& mu() const
            {
                return mu_;
            }

            volScalarField& mu()
            {
                return mu_;
            }

            const volScalarField& Pr() const
            {
                return Pr_;
            }

            volScalarField& Pr()
            {
                return Pr_;
            }

            const volScalarField& powerDensity() const
            {
                return powerDensityNeutronicsToLiquid_;
            }

            volScalarField& powerDensity()
            {
                return powerDensityNeutronicsToLiquid_;
            }

            const volScalarField& contErr() const
            {
                return contErr_;
            }

            volScalarField& contErr()
            {
                return contErr_;
            }

            const scalar& cumulContErr() const
            {
                return cumulContErr_;
            }

            scalar& cumulContErr()
            {
                return cumulContErr_;
            }

            const volScalarField& Dh() const
            {
                return Dh_;
            }

            volScalarField& Dh()
            {
                return Dh_;
            }

            const scalarField& dispersion() const
            {
                return dispersionPtr_();
            }

            scalarField& dispersion()
            {
                return dispersionPtr_();
            }

            const dimensionedScalar thermoResidualAlpha() const
            {
                return thermoResidualAlpha_;
            }

            autoPtr<volScalarField>& aboveThermoResidualAlphaPtr()
            {
                return aboveThermoResidualAlphaPtr_;
            }

            autoPtr<volScalarField>& belowThermoResidualAlphaPtr()
            {
                return belowThermoResidualAlphaPtr_;
            }

            const volScalarField& aboveThermoResidualAlpha() const
            {
                return aboveThermoResidualAlphaPtr_();
            }

            const volScalarField& belowThermoResidualAlpha() const
            {
                return belowThermoResidualAlphaPtr_();
            }

            bool Boussinesq() const
            {
                return Boussinesq_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
