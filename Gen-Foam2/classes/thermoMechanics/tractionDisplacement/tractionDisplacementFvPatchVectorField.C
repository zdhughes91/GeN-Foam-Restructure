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

#include "tractionDisplacementFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

tractionDisplacementFvPatchVectorField::
tractionDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(p, iF),
    traction_(p.size(), vector::zero),
    pressure_(p.size(), 0.0)
{
    fvPatchVectorField::operator=(patchInternalField());
    gradient() = vector::zero;
}


tractionDisplacementFvPatchVectorField::
tractionDisplacementFvPatchVectorField
(
    const tractionDisplacementFvPatchVectorField& tdpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchVectorField(tdpvf, p, iF, mapper),
    traction_(tdpvf.traction_, mapper),
    pressure_(tdpvf.pressure_, mapper)
{}


tractionDisplacementFvPatchVectorField::
tractionDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchVectorField(p, iF),
    traction_("traction", dict, p.size()),
    pressure_("pressure", dict, p.size())
{
    fvPatchVectorField::operator=(patchInternalField());
    gradient() = vector::zero;
}


tractionDisplacementFvPatchVectorField::
tractionDisplacementFvPatchVectorField
(
    const tractionDisplacementFvPatchVectorField& tdpvf
)
:
    fixedGradientFvPatchVectorField(tdpvf),
    traction_(tdpvf.traction_),
    pressure_(tdpvf.pressure_)
{}


tractionDisplacementFvPatchVectorField::
tractionDisplacementFvPatchVectorField
(
    const tractionDisplacementFvPatchVectorField& tdpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(tdpvf, iF),
    traction_(tdpvf.traction_),
    pressure_(tdpvf.pressure_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void tractionDisplacementFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchVectorField::autoMap(m);
    traction_.autoMap(m);
    pressure_.autoMap(m);
}


void tractionDisplacementFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchVectorField::rmap(ptf, addr);

    const tractionDisplacementFvPatchVectorField& dmptf =
        refCast<const tractionDisplacementFvPatchVectorField>(ptf);

    traction_.rmap(dmptf.traction_, addr);
    pressure_.rmap(dmptf.pressure_, addr);
}


void tractionDisplacementFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const dictionary& thermoMechanicalProperties =
        db().lookupObject<IOdictionary>("thermoMechanicalProperties");

    //const dictionary& thermalProperties =
    //    db().lookupObject<IOdictionary>("thermalProperties");


    const fvPatchField<scalar>& rho =
        patch().lookupPatchField<volScalarField, scalar>("rho");

    const fvPatchField<scalar>& rhoE =
        patch().lookupPatchField<volScalarField, scalar>("E");

    const fvPatchField<scalar>& nu =
        patch().lookupPatchField<volScalarField, scalar>("nu");

    scalarField E(rhoE/rho);
    scalarField mu(E/(2.0*(1.0 + nu)));
    scalarField lambda(nu*E/((1.0 + nu)*(1.0 - 2.0*nu)));
    scalarField threeK(E/(1.0 - 2.0*nu));

    Switch planeStress(thermoMechanicalProperties.lookup("planeStress"));

    if (planeStress)
    {
        lambda = nu*E/((1.0 + nu)*(1.0 - nu));
        threeK = E/(1.0 - nu);
    }

    scalarField twoMuLambda(2*mu + lambda);

    vectorField n(patch().nf());

    const fvPatchField<symmTensor>& sigmaD =
        patch().lookupPatchField<volSymmTensorField, symmTensor>("sigmaD");

    gradient() =
    (
        (traction_ - pressure_*n)/rho
      + twoMuLambda*fvPatchField<vector>::snGrad() - (n & sigmaD)
    )/twoMuLambda;

    //Switch thermalStress(thermalProperties.lookup("thermalStress"));


    //if (thermalStress)
    {
        const fvPatchField<scalar>&  threeKalpha=
            patch().lookupPatchField<volScalarField, scalar>("threeKalpha");

        const fvPatchField<scalar>& Tmech =
            patch().lookupPatchField<volScalarField, scalar>("TStruct");

		dimensionedScalar TrefStructure
        (
            "",
            dimTemperature,
            thermoMechanicalProperties.get<scalar>("TStructRef")
        ); //modified
		//fvPatchField<scalar> TrefStructurePatch =
		//patch().lookupPatchField<volScalarField, scalar>("Tmech");
		//fvPatchField<scalar> TrefStructurePatch = 
            //(Tmech/Tmech)*TrefStructure.value();

        gradient() += n*threeKalpha*(Tmech-TrefStructure.value())/twoMuLambda; //modified
        //gradient() += n*threeKalpha*(Tmech)/twoMuLambda;
    }

    fixedGradientFvPatchVectorField::updateCoeffs();
}


void tractionDisplacementFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    traction_.writeEntry("traction", os);
    pressure_.writeEntry("pressure", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    tractionDisplacementFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
