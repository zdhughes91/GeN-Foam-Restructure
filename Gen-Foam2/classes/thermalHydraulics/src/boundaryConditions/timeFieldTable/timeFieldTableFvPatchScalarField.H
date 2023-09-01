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
    Foam::timeFieldTableFvPatchScalarField

Group
    grpInletBoundaryConditions

Description
    Boundary condition that allows to specify a scalarField over the boundary
    that changes in time via a table entry in the form
    
    \verbatim
    table
    (
        (t0 (v0_0 v1_0 v2_0 ...))
        ...
        (tN (v0_N v1_N v2_N ...))
    );
    \endverbatim
    
    with ti being the time value at whiche the field vj_i is to be applied, 
    with j looping over face cells according to the patch indices. Field
    vj_0 will be applied for times t < t0 and field vj_N will be applied
    for times t > tN.
    This functionality is currently unavailable in OpenFOAM and I could not
    even resort to the Table class while coding this, as it is templated and 
    does not support the scalarField type as template parameter.

    Author: 
    Stefan Radman (stefan.radman@epfl.ch; stefanradman92@gmail.com)

Usage
    \verbatim
    boundaryName
    {
        type timeFieldTable;
        table
        (
            (t0 (v0_0 v1_0 v2_0 ...))
            (t1 (v0_1 v1_1 v2_1 ...))
            ...
        );
    }
    \endverbatim

SourceFiles
    timeFieldTableFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef timeFieldTableFvPatchScalarField_H
#define timeFieldTableFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class timeFieldTableFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class timeFieldTableFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

    int currentIndex_;

    List<Tuple2<scalar, scalarField>> table_;

    scalar tStart_;

public:

   //- Runtime type information
   TypeName("timeFieldTable");


   // Constructors

        //- Construct from patch and internal field
        timeFieldTableFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        timeFieldTableFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  timeFieldTableFvPatchScalarField
        //  onto a new patch
        timeFieldTableFvPatchScalarField
        (
            const timeFieldTableFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        timeFieldTableFvPatchScalarField
        (
            const timeFieldTableFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new timeFieldTableFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        timeFieldTableFvPatchScalarField
        (
            const timeFieldTableFvPatchScalarField&,
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
                new timeFieldTableFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        //- Update patchField with table values and linearly interpolate fields
        //  if necessary
        void updateField();

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
