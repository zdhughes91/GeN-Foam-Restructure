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
    Foam::FFPair

Description
    Class to handle dimensionless numbers that characterise the system, e.g.
    Reynolds, etc. While some can (rightfully) think that some of these numbers
    should be directly members of other classes (e.g., Reynolds of fluids
    should be in the fluid class), in general these numbers depend on a
    multitude of parameters. You cannot define a Reynolds for a fluid by
    itself, it needs to be geometrically constratined in a certain domain 
    (so that a hydraulic diameter can be defined), so it is not properly
    "just" of the fluid. While this logic can of course be debate, I wanted
    to keep it simple and not start dispersing dimensionless numbers randomly
    around the code. 

    Author: 
    Stefan Radman (sradman@protonmail.com; stefan.radman@epfl.ch)

SourceFiles
    FFPair.C

\*---------------------------------------------------------------------------*/

#ifndef FFPair_H
#define FFPair_H

#include "customPimpleControl.H"
#include "fluid.H"
#include "virtualMass.H"
#include "FFDragFactor.H"
#include "dispersionModel.H"
#include "interfacialAreaModel.H"
#include "FFHeatTransferCoefficientModel.H"
#include "phaseChangeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class FFPair Declaration
\*---------------------------------------------------------------------------*/

class FFPair
    : public regIOobject
{
protected :

    //- Ref to mesh
    const fvMesh& mesh_;

    //- Ref to phaseProperties dict
    const dictionary& dict_;

    //-
    const customPimpleControl& pimple_;

    //- Ref to fluid1
    fluid& fluid1_;

    //- Ref to fluid2
    fluid& fluid2_;

    //- Phase fraction of the dispersed phase
    volScalarField alphaDispersed_;

    //- Phase fraction of the continuous phase
    volScalarField alphaContinuous_;

    //- Hydraulic diameter of the dispersed phase
    volScalarField DhDispersed_;

    //- Hydraulic diameter of the continuous phase
    volScalarField DhContinuous_;

    //- Density of the continuous phase
    volScalarField rhoContinuous_;

    //- Effective momentum diffusivity of the continuous phase
    volScalarField nuContinuous_;

    //- Pr of the continuous phase
    volScalarField PrContinuous_;

    //- Reynolds of the pair
    volScalarField Re_;

    //- Minimum allowable Reynolds
    dimensionedScalar minRe_;

    //- Magnitude of relative velocity difference between phases
    volScalarField magUr_;

    //- Fluid-fluid interfacial drag factor model
    autoPtr<FFDragFactor> KdPtr_;

    //- Fluid-fluid interfacial drag factor
    volScalarField Kd_;

    //-
    autoPtr<virtualMass> virtualMassPtr_;

    //- Heat transfer coefficient models
    autoPtr<FFHeatTransferCoefficientModel> htc1Ptr_;
    autoPtr<FFHeatTransferCoefficientModel> htc2Ptr_;

    //- Table of volScalarFields (for htcs)
    //- Heat transfer coefficients between bulk of each fluid and interface
    typedef HashTable<volScalarField*, word, word::hash> 
        volScalarFieldPtrTable;
    volScalarFieldPtrTable htcs_;

    //- Linearized semi-implicit heat sources due to interfacial heat and mass
    //  transfer (inclusive of intrinsic mass transfer enthalpy contribution)
    typedef HashTable<fvScalarMatrix*, word, word::hash> 
        fvScalarMatrixPtrTable;
    fvScalarMatrixPtrTable heSources_;

    //- Fluid dispersion model (controls which phase is dispersed and which is
    //  continuous)
    autoPtr<dispersionModel> dispersion1Ptr_;

    //- Phase change model
    autoPtr<phaseChangeModel> phaseChangePtr_;

    //- Interfacial temperature
    volScalarField iT_;

    //- Interfacial area model
    autoPtr<interfacialAreaModel> iAPtr_;

    //- Interfacial area
    volScalarField iA_;

    //- This sucks... A flag just to avoid u-relaxing the interfacial T (if not
    //  set by a phasChangeModel) on the first time step. It is by far the
    //  simplest but also less elegant way to do it
    bool firstTimeStepAndIter_;

public:

    TypeName("FFPair");

    //- Constructors

        FFPair
        (
            fluid& fluid1,
            fluid& fluid2,
            const dictionary& dict
        );

        //- Destructor
        virtual ~FFPair(){}

    //- Member Functions

        //- To remove abstraction from regIOobject
        virtual bool writeData(Ostream& os) const
        {
            return os.good();
        }

        //- Add virtual mass contribution to UEqn
        void addVirtualMassForce(fvVectorMatrix& UEqn) const;

        //- 
        void correct
        (
            const bool& correctFluidDynamics, 
            const bool& correctEnergy
        );

        //- Access

            const fvMesh& mesh() const
            {
                return mesh_;
            }

            const customPimpleControl& pimple() const
            {
                return pimple_;
            }

            const fluid& fluid1() const
            {
                return fluid1_;
            }

            const fluid& fluid2() const
            {
                return fluid2_;
            }

            word reversedName() const
            {
                return fluid2_.name()+"."+fluid1_.name();
            }

            const fluid* getFluidPtr(const word& name) const
            {
                if (fluid1_.name() == name)
                    return &fluid1_;
                else if (fluid2_.name() == name)
                    return &fluid2_;
                else
                    return nullptr;
            }

            const fluid* getOtherFluidPtr(const word& name) const
            {
                if (fluid1_.name() == name)
                    return &fluid2_;
                else if (fluid2_.name() == name)
                    return &fluid1_;
                else
                    return nullptr;
            }

            bool foundUnorderedPairSubDict
            (
                const dictionary& dict
            ) const
            {
                return 
                    (
                        dict.found(this->name())
                    or  dict.found(this->reversedName())
                    );
            }

            const dictionary& getUnorderedPairSubDict
            (
                const dictionary& dict
            ) const
            {
                return 
                    (dict.found(this->name())) ?
                    dict.subDict(this->name()) :
                    dict.subDict(this->reversedName());
            }

            const volScalarField& alphaDispersed() const
            {
                return alphaDispersed_;
            }

            const volScalarField& alphaContinuous() const
            {
                return alphaContinuous_;
            }

            const volScalarField& DhDispersed() const
            {
                return DhDispersed_;
            }

            const volScalarField& DhContinuous() const
            {
                return DhContinuous_;
            }

            const volScalarField& rhoContinuous() const
            {
                return rhoContinuous_;
            }

            const volScalarField& nuContinuous() const
            {
                return nuContinuous_;
            }

            const volScalarField& PrContinuous() const
            {
                return PrContinuous_;
            }

            const volScalarField& Re() const
            {
                return Re_;
            }       

            const volScalarField& magUr() const
            {
                return magUr_;
            }

            const volScalarField& Kd() const
            {
                return Kd_;
            }

            //- Virtual mass coefficient
            const volScalarField& Vm() const
            {
                return virtualMassPtr_->Vm();
            }

            volScalarField& htc(const word& fluidName)
            {
                return *htcs_[fluidName];
            }

            const volScalarField& htc(const word& fluidName) const
            {
                return *htcs_[fluidName];
            }

            const fvScalarMatrix& heSource(const word& heName) const
            {
                return *heSources_[heName];
            }

            bool phaseChange() const
            {
                return phaseChangePtr_.valid();
            }

            //- Latent heat
            const volScalarField& L() const
            {
                return phaseChangePtr_->L();
            }

            const volScalarField& dmdt() const
            {
                return phaseChangePtr_->dmdt();
            }

            volScalarField& iT()
            {
                return iT_;
            }

            const volScalarField& iT() const
            {
                return iT_;
            }

            volScalarField& iA()
            {
                return iA_;
            }

            const volScalarField& iA() const
            {
                return iA_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
