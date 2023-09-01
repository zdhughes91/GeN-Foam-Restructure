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
    Foam::RASModels::LaheyKEpsilon

Group
    grpRASTurbulence

Description
    Continuous-phase k-epsilon model including bubble-generated turbulence.

    Reference:
    \verbatim
        Lahey Jr, R. T. (2005).
        The simulation of multidimensional multiphase flows.
        Nuclear Engineering and Design, 235(10), 1043-1060.
    \endverbatim

    The default model coefficients are
    \verbatim
        LaheyKEpsilonCoeffs
        {
            Cmu             0.09;
            C1              1.44;
            C2              1.92;
            C3              -0.33;
            sigmak          1.0;
            sigmaEps        1.3;
            Cp              0.25;
            Cmub            0.6;
            alphaInversion  0.3;
        }
    \endverbatim

    Last modified by Stefan Radman for compatibility with FFSEulerFoam
    (stefan.radman@epfl.ch, stefanradman92@gmail.com)

SourceFiles
    LaheyKEpsilon.C

Note
    This model will work as expected only in scenarios where the gas is
    dispersed and the liquid is continuous! While it will still technically
    run also in opposite scenarios, the validity of the results is 
    questionable. This is due to the fact that the fluid-fluid drag coefficient
    Cd is computed by taking into consideration the dispersed Dh and continuous
    rho, while the source terms in this model are hard-coded to use the
    gas Dh and liquid rho in the different terms. Why? Read the paper from
    Lahey, this model was devised only for these dipsersed gas / continuous
    liquid scenarios

\*---------------------------------------------------------------------------*/

#ifndef LaheyKEpsilon_H
#define LaheyKEpsilon_H

#include "kEpsilon.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fluid;
class FFPair;

namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class LaheyKEpsilon Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class LaheyKEpsilon
:
    public kEpsilon<BasicTurbulenceModel>
{
    // Private data
    
        //- Phase names
        mutable word liquidName_;
        mutable word gasName_; 
        
        //- Field pointers. As these fields are expected to be constructed at 
        //  different times in the solver, they can't be looked-up at 
        //  construction. The private access functions serve that purpose
        mutable const fluid* gasPtr_;
        mutable const fluid* liquidPtr_;
        mutable const FFPair* pairPtr_;

    //- Private functions

        //- No copy construct
        LaheyKEpsilon(const LaheyKEpsilon&) = delete;

        //- No copy assignment
        void operator=(const LaheyKEpsilon&) = delete;

        //- Set liquidName_ and gasName_
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

            dimensionedScalar alphaInversion_;
            dimensionedScalar Cp_;
            dimensionedScalar C3_;
            dimensionedScalar Cmub_;


    // Protected Member Functions

        virtual void correctNut();
        tmp<volScalarField> bubbleG() const;
        tmp<volScalarField> phaseTransferCoeff() const;
        virtual tmp<fvScalarMatrix> kSource() const;
        virtual tmp<fvScalarMatrix> epsilonSource() const;


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("LaheyKEpsilon");


    // Constructors

        //- Construct from components
        LaheyKEpsilon
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
    virtual ~LaheyKEpsilon()
    {}


    // Member Functions

        //- Read model coefficients if they have changed
        virtual bool read();

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "LaheyKEpsilon.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
