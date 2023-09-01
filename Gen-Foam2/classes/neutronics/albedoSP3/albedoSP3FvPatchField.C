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

\*---------------------------------------------------------------------------*/

#include "albedoSP3FvPatchField.H"
#include "dictionary.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "uniformDimensionedFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
albedoSP3FvPatchField<Type>::albedoSP3FvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(p, iF),
    gamma_(scalar(0.5)),
    diffCoeffName_("diffCoeffName"),
    fluxStarAlbedo_("fluxStarAlbedo"),
    forSecondMoment_(false)
{}


template<class Type>
albedoSP3FvPatchField<Type>::albedoSP3FvPatchField
(
    const albedoSP3FvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fvPatchField<Type>(ptf, p, iF, mapper),
    gamma_(ptf.gamma_),
    diffCoeffName_(ptf.diffCoeffName_),
    fluxStarAlbedo_(ptf.fluxStarAlbedo_),
    forSecondMoment_(ptf.forSecondMoment_)
{}

 
template<class Type>
albedoSP3FvPatchField<Type>::albedoSP3FvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fvPatchField<Type>(p, iF, dict),
    gamma_(dict.lookupOrDefault<scalar>("gamma", scalar(0.5))),
    diffCoeffName_(dict.lookupOrDefault<word>("diffCoeffName", "diffCoeffName")),
    fluxStarAlbedo_(dict.lookupOrDefault<word>("fluxStarAlbedo", "fluxStarAlbedo")),
    forSecondMoment_(dict.lookupOrDefault<bool>("forSecondMoment", false))
{
    evaluate();
}


template<class Type>
albedoSP3FvPatchField<Type>::albedoSP3FvPatchField
(
    const albedoSP3FvPatchField<Type>& ptf
)
:
    fvPatchField<Type>(ptf),
    gamma_(ptf.gamma_),
    diffCoeffName_(ptf.diffCoeffName_),
    fluxStarAlbedo_(ptf.fluxStarAlbedo_),
    forSecondMoment_(ptf.forSecondMoment_)
{}


template<class Type>
albedoSP3FvPatchField<Type>::albedoSP3FvPatchField
(
    const albedoSP3FvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(ptf, iF),
    gamma_(ptf.gamma_),
    diffCoeffName_(ptf.diffCoeffName_),
    fluxStarAlbedo_(ptf.fluxStarAlbedo_),
    forSecondMoment_(ptf.forSecondMoment_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void albedoSP3FvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fvPatchField<Type>::autoMap(m);
}


template<class Type>
void albedoSP3FvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fvPatchField<Type>::rmap(ptf, addr);

    const albedoSP3FvPatchField<Type>& fgptf =
        refCast<const albedoSP3FvPatchField<Type> >(ptf);

}


template<class Type>
void albedoSP3FvPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    const Field<scalar>& diffCoeff =
        this->patch().template lookupPatchField<volScalarField, scalar>
        (
            diffCoeffName_
        );
    const Field<scalar> fluxStarAlbedo =
        this->patch().template lookupPatchField<volScalarField, scalar>
        (
            fluxStarAlbedo_
        );

    if(forSecondMoment_)
    {
        Field<Type>::operator=
        (
            this->patchInternalField() 
        + (this->patchInternalField()*gamma_/diffCoeff*21.0/20.0)/this->patch().deltaCoeffs()
        -  (Type(pTraits<Type>::one)*fluxStarAlbedo*gamma_/(diffCoeff*27.0/35.0)*3.0/20.0)/this->patch().deltaCoeffs()
        );
    }
    else
    {
        Field<Type>::operator=
        (
            this->patchInternalField() 
        + (this->patchInternalField()*gamma_/diffCoeff)/this->patch().deltaCoeffs()
        -  (Type(pTraits<Type>::one)*fluxStarAlbedo*gamma_/diffCoeff*3.0/4.0)/this->patch().deltaCoeffs()
        );
    }

    fvPatchField<Type>::evaluate();
}


template<class Type>
tmp<Field<Type> > albedoSP3FvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const

{

    const Field<scalar>& diffCoeff =
        this->patch().template lookupPatchField<volScalarField, scalar>
        (
            diffCoeffName_
        );
    const Field<scalar>& fluxStarAlbedo =
        this->patch().template lookupPatchField<volScalarField, scalar>
        (
            fluxStarAlbedo_
        );

    if(forSecondMoment_)
    {
        return Type(pTraits<Type>::one) * (1 + gamma_/(diffCoeff*27.0/35.0)*21.0/20.0/this->patch().deltaCoeffs());

    }
    else
    {
        return Type(pTraits<Type>::one) * (1 + gamma_/diffCoeff/this->patch().deltaCoeffs());
    }

    
}

template<class Type>
tmp<Field<Type> > albedoSP3FvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const

{
    const Field<scalar>& diffCoeff =
        this->patch().template lookupPatchField<volScalarField, scalar>
        (
            diffCoeffName_
        );
    const Field<scalar>& fluxStarAlbedo =
        this->patch().template lookupPatchField<volScalarField, scalar>
        (
            fluxStarAlbedo_
        );
 
    if(forSecondMoment_)
    {
        return -Type(pTraits<Type>::one) * (fluxStarAlbedo*gamma_/(diffCoeff*27.0/35.0)*3.0/20.0/this->patch().deltaCoeffs());

    }
    else
    {
        return -Type(pTraits<Type>::one) * (fluxStarAlbedo*gamma_/diffCoeff*3.0/4.0/this->patch().deltaCoeffs());
    }


}


template<class Type>
tmp<Field<Type> > albedoSP3FvPatchField<Type>::
gradientInternalCoeffs() const
{
    const Field<scalar>& diffCoeff =
        this->patch().template lookupPatchField<volScalarField, scalar>
        (
            diffCoeffName_
        );
    const Field<scalar>& fluxStarAlbedo =
        this->patch().template lookupPatchField<volScalarField, scalar>
        (
            fluxStarAlbedo_
        );

    if(forSecondMoment_)
    {
        return Type(pTraits<Type>::one)*(- gamma_/(diffCoeff*27.0/35.0)*21.0/20.0);

    }
    else
    {
        return Type(pTraits<Type>::one)*(- gamma_/diffCoeff);
    }

}


template<class Type>
tmp<Field<Type> > albedoSP3FvPatchField<Type>::
gradientBoundaryCoeffs() const
{
    const Field<scalar>& diffCoeff =
        this->patch().template lookupPatchField<volScalarField, scalar>
        (
            diffCoeffName_
        );
    const Field<scalar>& fluxStarAlbedo =
        this->patch().template lookupPatchField<volScalarField, scalar>
        (
            fluxStarAlbedo_
        );

    if(forSecondMoment_)
    {
        return Type(pTraits<Type>::one)*fluxStarAlbedo*gamma_/(diffCoeff*27.0/35.0)*3.0/20.0;

    }
    else
    {
        return Type(pTraits<Type>::one)*fluxStarAlbedo*gamma_/diffCoeff*3.0/4.0;
    }
}

template<class Type>
tmp<Field<Type> > albedoSP3FvPatchField<Type>::snGrad() const
{
    const Field<scalar>& diffCoeff =
        this->patch().template lookupPatchField<volScalarField, scalar>
        (
            diffCoeffName_
        ); 
    const Field<scalar>& fluxStarAlbedo =
        this->patch().template lookupPatchField<volScalarField, scalar>
        (
            fluxStarAlbedo_
        );

    return (- (this->patchInternalField()*gamma_/diffCoeff));
}

template<class Type>
void albedoSP3FvPatchField<Type>::write(Ostream& os) const
{

    fvPatchField<Type>::write(os);

    os.writeKeyword("gamma")
        << gamma_ << token::END_STATEMENT << nl;   
    os.writeKeyword("diffCoeffName")
        << diffCoeffName_ << token::END_STATEMENT << nl;
    os.writeKeyword("fluxStarAlbedo")
        << fluxStarAlbedo_ << token::END_STATEMENT << nl;
    os.writeKeyword("forSecondMoment")
        << forSecondMoment_ << token::END_STATEMENT << nl;    

    this->writeEntry("value", os);

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeFieldTypedefs(albedoSP3);

makePatchFields(albedoSP3);
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam




// ************************************************************************* //
