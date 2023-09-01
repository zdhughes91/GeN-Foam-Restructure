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
    Foam::compressible::NusseltThermalBaffle1DFvPatchScalarField

Description
    This BC models heat conduction through a thin 1-D baffle thermally coupling
    two parts of a mesh. By 'thin' we mean that the thermal inertia of the
    baffle is neglected and that baffle wall temperatures on both sides are
    predicted via an instantaneous heat flux balance. The fluid-baffle heat
    transfer coefficients for wall-fluid (and the other wall-fluid pair if 
    doing twoPhase simulations) are calculated starting from the Nusselt number
    calculated in a form Nu = const + coeff*(Re^expRe)*(Pr^expPr) with user
    provided coefficients. The update for the fields on both sides of the
    baffle is implicit, meaning that the owner patchField also sets the 
    valueFraction/refValue/refGradient of the nbrPatch. In two-phase scenarios,
    the owner of the patchField associated with the first instantiation of a 
    fluid type object also sets the aforementioned fields for the owner
    and nbr patchFields of the other fluid temperature. This is done to
    preserve the implicitness of the method and to avoid useless extra 
    caculations. Encapsulation flies out of the window but encapsulation is a 
    C++ idea that, if followed with pseudo-religious fervor, leads to shit
    anyways!

    Author: 
    Stefan Radman (sradman@pm.me / stefan.radman@epfl.ch)

Usage
    \verbatim
    <masterPatchName>
    {
        type            NusseltThermalBaffle1D;
        value           $internalField;
        samplePatch     <slavePatchName>;

        thickness       0.005;  //- Required, baffle thickness
        kappa           20;     //- Required, baffle thermal conductivity
        const           4;      //- Required
        coeff           0.018;  //- Required
        expRe           0.8;    //- Required
        expPr           0.8;    //- Required
    }

    <slavePatchName>
    {
        type            NusseltThermalBaffle1D;
        value           $internalField;
        samplePatch     <masterPatchName>;
        
        const           4;      //- If not provided, defaults to master value
        coeff           0.018;  //- If not provided, defaults to master value
        expRe           0.8;    //- If not provided, defaults to master value
        expPr           0.8;    //- If not provided, defaults to master value
    }
    \endverbatim

SourceFiles
    NusseltThermalBaffle1DFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef NusseltThermalBaffle1DFvPatchScalarField_H
#define NusseltThermalBaffle1DFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "autoPtr.H"
#include "mappedPatchBase.H"

#include "FSPair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
            Class NusseltThermalBaffle1DFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class NusseltThermalBaffle1DFvPatchScalarField
:
    public mappedPatchBase,
    public mixedFvPatchScalarField
{
    // Private data

        //- Name of the temperature field
        word TName_;

        //- True for twoPhase simulations, false for onePhase ones
        bool twoPhase_;

        //- True if this fvPatchScalarField is supposed to set the
        //  valueFraction, refValue and refGradient of ALL the involved
        //  fvPatchScalarFields, namely on both sides of the baffle for this
        //  fluid (i.e. this and nbrPatchFieldPtr_) and the other fluid (i.e.
        //  otherFluidPatchFieldPtr_ and otherFluidNbrPatchFieldPtr_
        bool twoPhaseOwner_;

        //- Ptrs
        const fluid* fluidPtr_;
        const fluid* otherFluidPtr_;
        const FSPair* pairPtr_;
        const FSPair* otherPairPtr_;
        NusseltThermalBaffle1DFvPatchScalarField* nbrPatchFieldPtr_;
        NusseltThermalBaffle1DFvPatchScalarField* otherFluidPatchFieldPtr_;
        NusseltThermalBaffle1DFvPatchScalarField* otherFluidNbrPatchFieldPtr_;

        //- Nusselt correlation parameters
        scalar const_;
        scalar coeff_;
        scalar expRe_;
        scalar expPr_;

        //- Wall properties
        scalar tw_;
        scalar kw_;

        //- Wall heat transfer coefficient
        scalar hw_;

    // Private members

        //-
        void setPtrs();

        //- Is Owner
        bool owner() const;


public:

    //- Runtime type information
    TypeName("NusseltThermalBaffle1D");


    // Constructors

        //- Construct from patch and internal field
        NusseltThermalBaffle1DFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        NusseltThermalBaffle1DFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  NusseltThermalBaffle1DFvPatchScalarField onto a new patch
        NusseltThermalBaffle1DFvPatchScalarField
        (
            const NusseltThermalBaffle1DFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        NusseltThermalBaffle1DFvPatchScalarField
        (
            const NusseltThermalBaffle1DFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new NusseltThermalBaffle1DFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        NusseltThermalBaffle1DFvPatchScalarField
        (
            const NusseltThermalBaffle1DFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new NusseltThermalBaffle1DFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        //-
        void readDict();

        // Access

        bool twoPhaseOwner() const
        {
            return twoPhaseOwner_;
        }
        bool& twoPhaseOwner()
        {
            return twoPhaseOwner_;
        }
        scalar cnst() const
        {
            return const_;
        }
        scalar coeff()  const
        {
            return coeff_;
        }
        scalar expRe()  const
        {
            return expRe_;
        }
        scalar expPr() const
        {
            return expPr_;
        }
        scalar tw() const
        {
            return tw_;
        }
        scalar kw() const
        {
            return kw_;
        }
        scalar hw() const
        {
            return hw_;
        }

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );

        //- Evaluate the patch field
        virtual void evaluate
        (
            const Pstream::commsTypes commsType =
                Pstream::commsTypes::blocking
        ) override;

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        tmp<scalarField> calcH
        (
            const FSPair& pair,
            const NusseltThermalBaffle1DFvPatchScalarField& patchField            
        );

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
