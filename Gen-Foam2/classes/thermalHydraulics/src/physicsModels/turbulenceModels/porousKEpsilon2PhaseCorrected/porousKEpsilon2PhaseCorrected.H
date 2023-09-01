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
    Foam::RASModels::porousKEpsilon2PhaseCorrected

Group
    grpRASTurbulence

Description
   Same as porousKEpsiolon but turbulent intensity correlation in the form 
   turbulenceIntensityCoeff*Reynolds^turbulenceIntensityExp + 
   turbulenceIntensityAlphaCoeff * 
   * (1-alpha), where alpha is the phase fraction. 

SourceFiles
    porousKEpsilon2PhaseCorrected.C

\*---------------------------------------------------------------------------*/

#ifndef porousKEpsilon2PhaseCorrected_H
#define porousKEpsilon2PhaseCorrected_H

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
                           Class porousKEpsilon2PhaseCorrected Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class porousKEpsilon2PhaseCorrected
:
    public eddyViscosity<RASModel<BasicTurbulenceModel>>
{
    // Private Member Functions

        //- No copy construct
        porousKEpsilon2PhaseCorrected(const porousKEpsilon2PhaseCorrected&) = delete;

        //- No copy assignment
        void operator=(const porousKEpsilon2PhaseCorrected&) = delete;


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
            dictionary porousKEpsilon2PhaseCorrectedDict_;

            //- Characteristic length over which k and epsilon converge to 
            //  their equilibrium values
            volScalarField convergenceLength_;

            //- A in I_t = A*Re^B + D*(1-alpha)
            volScalarField turbulenceIntensityCoeff_;

            //- B in I_t = A*Re^B + D*(1-alpha)
            volScalarField turbulenceIntensityExp_;

            //- D in I_t = A*Re^B + D*(1-alpha)
            volScalarField turbulenceIntensityAlphaCoeff_;

            //- C in L_t = C*Dh
            volScalarField turbulenceLengthScaleCoeff_;

            //- ...
            volScalarField equilibriumEpsilon_;
            volScalarField equilibriumK_;

            //- Fields for nutStabilization as seen in the old GeN-Foam 
            //  porousKEpsilon2PhaseCorrected
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
    TypeName("porousKEpsilon2PhaseCorrected");


    // Constructors

        //- Construct from components
        porousKEpsilon2PhaseCorrected
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
    virtual ~porousKEpsilon2PhaseCorrected()
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
    #include "porousKEpsilon2PhaseCorrected.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
