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
    Foam::RASModels::mixtureKEpsilon

Group
    grpRASTurbulence

Description
    Mixture k-epsilon turbulence model for two-phase gas-liquid systems

    The basic structure of the model is based on:
    \verbatim
        Behzadi, A., Issa, R. I., & Rusche, H. (2004).
        Modelling of dispersed bubble and droplet flow at high phase fractions.
        Chemical Engineering Science, 59(4), 759-770.
    \endverbatim

    but an effective density for the gas is used in the averaging and an
    alternative model for bubble-generated turbulence from:
    \verbatim
        Lahey Jr, R. T. (2005).
        The simulation of multidimensional multiphase flows.
        Nuclear Engineering and Design, 235(10), 1043-1060.
    \endverbatim

    The default model coefficients are
    \verbatim
        mixtureKEpsilonCoeffs
        {
            Cmu         0.09;
            C1          1.44;
            C2          1.92;
            C3          C2;
            sigmak      1.0;
            sigmaEps    1.3;
        }
    \endverbatim

    Last modified by Stefan Radman for compatibility with FFSEulerFoam
    (stefan.radman@epfl.ch, stefanradman92@gmail.com)

SourceFiles
    mixtureKEpsilon.C


\*---------------------------------------------------------------------------*/

#ifndef mixtureKEpsilon_H
#define mixtureKEpsilon_H

#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fluid;
class FFPair;

namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class mixtureKEpsilon Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class mixtureKEpsilon
:
    public eddyViscosity<RASModel<BasicTurbulenceModel>>
{
    // Private data

        //- Phase names
        mutable word liquidName_;
        mutable word gasName_;

        //- True if this turbulence models is of the gas phase
        bool isGas_;

        //- Ptrs to access global solver field
        mutable const fluid* liquidPtr_;
        mutable const fluid* gasPtr_;
        mutable const FFPair* pairPtr_;

    //- Private functions

        //- No copy assignment
        void operator=(const mixtureKEpsilon&) = delete;

        //- No copy construct
        mixtureKEpsilon(const mixtureKEpsilon&) = delete;

        //- Set liquidName_, gasName_
        void setFluidNames() const;

        //- Access-like

            const fluid& gas() const;
            
            const fluid& liquid() const;
            
            const FFPair& pair() const;

            //- This looks up the liquid-vapour friction factor Kd and
            //  starting from it, it inverts it to obtain what was the
            //  drag coefficient Cd. In particular
            //  Kd = 
            //      0.5*(alpha_c*alpha_d)/(alpha_c+alpha_d)*
            //      (rho_c/Dh_d)*mag(U_c-U_d)*Cd
            //
            //  with c and d standing for the continuous and dispersed
            //  phase respectively
            const tmp<volScalarField> Cd() const;


protected:

    // Protected data

        // Model coefficients

            dimensionedScalar Cmu_;
            dimensionedScalar C1_;
            dimensionedScalar C2_;
            dimensionedScalar C3_;
            dimensionedScalar Cp_;
            dimensionedScalar sigmak_;
            dimensionedScalar sigmaEps_;

        // Fields

            volScalarField k_;
            volScalarField epsilon_;

        // Mixture fields

            autoPtr<volScalarField> Ct2_;
            autoPtr<volScalarField> rhom_;
            autoPtr<volScalarField> km_;
            autoPtr<volScalarField> epsilonm_;


    // Protected Member Functions

        wordList epsilonBoundaryTypes(const volScalarField& epsilon) const;

        void correctInletOutlet
        (
            volScalarField& vsf,
            const volScalarField& refVsf
        ) const;

        void initMixtureFields();

        virtual void correctNut();

        tmp<volScalarField> Ct2() const;

        tmp<volScalarField> rholEff() const;
        tmp<volScalarField> rhogEff() const;
        tmp<volScalarField> rhom() const;

        tmp<volScalarField> mix
        (
            const volScalarField& fc,
            const volScalarField& fd
        ) const;

        tmp<volScalarField> mixU
        (
            const volScalarField& fc,
            const volScalarField& fd
        ) const;

        tmp<surfaceScalarField> mixFlux
        (
            const surfaceScalarField& fc,
            const surfaceScalarField& fd
        ) const;

        tmp<volScalarField> bubbleG() const;
        virtual tmp<fvScalarMatrix> kSource() const;
        virtual tmp<fvScalarMatrix> epsilonSource() const;

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff(const volScalarField& nutm) const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DkEff",
                    nutm/sigmak_
                )
            );
        }

        //- Return the effective diffusivity for epsilon
        tmp<volScalarField> DepsilonEff(const volScalarField& nutm) const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DepsilonEff",
                    nutm/sigmaEps_
                )
            );
        }

public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("mixtureKEpsilon");


    // Constructors

        //- Construct from components
        mixtureKEpsilon
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
    virtual ~mixtureKEpsilon()
    {}


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

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

        //- Return the (estimated) specific dissipation rate
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
    #include "mixtureKEpsilon.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
