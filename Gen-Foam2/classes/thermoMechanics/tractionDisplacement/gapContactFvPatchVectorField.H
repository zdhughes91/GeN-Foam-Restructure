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
    Foam::gapContactFvPatchVectorField

Description
    Fixed traction boundary condition for the standard linear elastic, fixed
    coefficient displacement equation.

    Modification of this file (and associated .C or included .H files)
    Carlo Fiorina <carlo.fiorina@outlook.com; carlo.fiorina@epfl.ch;>, EPFL (Switzerland)

SourceFiles
    gapContactFvPatchVectorField.C


\*---------------------------------------------------------------------------*/

#ifndef gapContactFvPatchVectorField_H
#define gapContactFvPatchVectorField_H

#include "tractionDisplacementFvPatchVectorField.H"
#include "regionCoupledBaseFvPatch.H"
//#include "gapGasModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class tractionDisplacementFvPatch Declaration
\*---------------------------------------------------------------------------*/

class gapContactFvPatchVectorField
:
    public tractionDisplacementFvPatchVectorField
{

    // Private Data

        //- Local reference to region couple patch
        const regionCoupledBaseFvPatch& regionCoupledPatch_;
        
        //- Penalty factor
        scalar penaltyFact_;

        //- POffse
        scalar offset_;
        
        //- Gap width [m]
        mutable scalarField gapWidth_;

        //- Interface pressure
        mutable scalarField interfaceP_;

    // Private functions
        
protected:

    // Protected functions
    
            //- Return the neighbour patchField
            inline const gapContactFvPatchVectorField& neighbour() const
            {
                return refCast<const gapContactFvPatchVectorField>
                (
                    regionCoupledPatch_.neighbFvPatch().lookupPatchField
                    <volVectorField, vector>
                    (
                        internalField().name()
                    )
                );
            }

            //- Return the gap width
            tmp<scalarField> gapWidth() const;

            //- Minimum patch stiffness for interface pressure evaluation
            scalar boundaryStiffness() const;
            
public:

    //- Runtime type information
    TypeName("gapContact");


    // Constructors

        //- Construct from patch and internal field
        gapContactFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        gapContactFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  gapContactFvPatchVectorField onto a new patch
        gapContactFvPatchVectorField
        (
            const gapContactFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        gapContactFvPatchVectorField
        (
            const gapContactFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new gapContactFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        gapContactFvPatchVectorField
        (
            const gapContactFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new gapContactFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access
        
            //- Interface pressure
            inline const scalarField& interfaceP() const
            {
                return interfaceP_;
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
                const fvPatchVectorField&,
                const labelList&
            );
            
        // Modify

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
