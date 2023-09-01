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
    Foam::RASModels::porousKEpsilon

Group
    grpRASTurbulence

Description
    Same as standard OpenFOAM. porousKEpsilon is provided as additional model.
    The only difference is that it forces k and epsilon to equilibrium values
    inside the porous zones. These equilibium values can be set in the
    porousKepsionProperties sub-dictionary here below. k and epsilon are 
    determined based on correlations for tubulent intensity (I) and lengh 
    scale (L). Turbulent intensity correlation in the form 
    turbulenceIntensityCoeff*Reynolds^turbulenceIntensityExp,
    with Reynolds number calculated by the thermal-hydraulic class, according 
    to the input data in phaseProperties. 

    Please notice that a porous medium simulation using the porousKEpsilon 
    model entails the risk of an unstable solution. This is due to the fact 
    that the turbulent viscosity will be that of the sub-scale structure, and 
    thus not enough to stabilize a solution on the length scale of the coarse 
    mesh. To address this problem, one can define the keyword DhStruct in 
    constant/fluidRegion/phaseProperties/dragModels.(nameOfPhase).structure.
    (nameOfCellZones). This keyword defines the hydraulic diameter of the whole
    porous structure. The code uses it to make sure the turbulent viscosity 
    results in a  laminar Reynolds number (defaulted to 500).

Usage
    The following sub-dictionary should be included in the turbulenceProperties 
    dictionary:
    \verbatim
    porousKEpsilonProperties
    {
        "zones of application"
        {
        
            convergenceLength           0.5; // k and epsilon will exponentially 
                         // converge to equilibrium according to this exponent
            turbulenceIntensityCoeff    0.16; 
            turbulenceIntensityExp      -0.125;
            turbulenceLengthScaleCoeff  0.07;
        }
    }
    \endverbatim

SourceFiles
    porousKEpsilon.C

\*---------------------------------------------------------------------------*/

#ifndef porousKEpsilon_H
#define porousKEpsilon_H

#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Forward declarations
class structure;
class FSPair;

namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class porousKEpsilon Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class porousKEpsilon
:
    public eddyViscosity<RASModel<BasicTurbulenceModel>>
{
    // Private Member Functions

        //- No copy construct
        porousKEpsilon(const porousKEpsilon&) = delete;

        //- No copy assignment
        void operator=(const porousKEpsilon&) = delete;


protected:

    // Protected data

        // Model coefficients

            dimensionedScalar Cmu_;
            dimensionedScalar C1_;
            dimensionedScalar C2_;
            dimensionedScalar C3_;
            dimensionedScalar sigmak_;
            dimensionedScalar sigmaEps_;

        // Fields

            volScalarField k_;
            volScalarField epsilon_;

        //- FFS related

            const word fluidName_;

            //- Reference to structure
            const structure& structure_;

            //- Reference to fluid-structure pair
            const FSPair& FSPair_;

            //- Sub-dict
            dictionary porousKEpsilonDict_;

            //- Characteristic length over which k and epsilon converge to 
            //  their equilibrium values
            volScalarField convergenceLength_;

            //- A in I_t = A*Re^B
            volScalarField turbulenceIntensityCoeff_;

            //- B in I_t = A*Re^B
            volScalarField turbulenceIntensityExp_;

            //- C in L_t = C*Dh
            volScalarField turbulenceLengthScaleCoeff_;

            //- ...
            volScalarField equilibriumEpsilon_;
            volScalarField equilibriumK_;

            //- Fields for nutStabilization as seen in the old GeN-Foam 
            //  porousKEpsilon
            bool nutStabilization_;
            autoPtr<volScalarField> DhStructPtr_;
            autoPtr<volScalarField> laminarReStructPtr_;

        //- Precomputed constants

            dimensionedScalar Cmu3by4_;


    // Protected Member Functions

        virtual void correctNut();
        virtual tmp<fvScalarMatrix> kSource() const;
        virtual tmp<fvScalarMatrix> epsilonSource() const;


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("porousKEpsilon");


    // Constructors

        //- Construct from components
        porousKEpsilon
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~porousKEpsilon()
    {}


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DkEff",
                    (this->nut_/sigmak_ + this->nu())
                )
            );
        }

        //- Return the effective diffusivity for epsilon
        tmp<volScalarField> DepsilonEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DepsilonEff",
                    (this->nut_/sigmaEps_ + this->nu())
                )
            );
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }

        //- Return the (estimated) specific dissipation rate. This function
        //  is required by the new v2006 OpenFOAM version, as it is pure
        //  virtual in EddyViscosity
        virtual tmp<volScalarField> omega() const
        {
            return tmp<volScalarField>::New
            (
                IOobject
                (
                    IOobject::groupName("omega", this->alphaRhoPhi_.group()),
                    this->runTime_.timeName(),
                    this->mesh_
                ),
                epsilon_/(Cmu_*k_)
            );
        }

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "porousKEpsilon.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
