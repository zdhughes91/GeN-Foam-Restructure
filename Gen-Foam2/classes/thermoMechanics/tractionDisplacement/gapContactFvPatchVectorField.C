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

#include "gapContactFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "regionCoupledFvPatch.H"
//#include "thermoMechanicsSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Private members  * * * * * * * * * * * * *//
    
Foam::tmp<Foam::scalarField> 
Foam::gapContactFvPatchVectorField::gapWidth() const
{
    const regionCoupledFvPatch& patch
        = refCast<const regionCoupledFvPatch>(regionCoupledPatch_);
    
    const regionCoupledFvPatch& nbrPatch
        = refCast<const regionCoupledFvPatch>(patch.neighbPatch());
        
    //- Include the total displacement
/*
    const thermoMechanicsSolver& thermoMechanics
        = patch.boundaryMesh().mesh().lookupObject<thermoMechanicsSolver>("thermoMechanics");        

    word patchName = patch.name();
    label patchID = patch.boundaryMesh().findPatchID(patchName); 
    const fvPatchVectorField& totalDispPatch = thermoMechanics.totalDisplacement().boundaryField()[patchID]; 
    
    word nbrPatchName = nbrPatch.name();
    label nbrPatchID = nbrPatch.boundaryMesh().findPatchID(nbrPatchName);     
    const fvPatchVectorField& totalDispNbrPatch = thermoMechanics.totalDisplacement().boundaryField()[nbrPatchID];     
*/
    volVectorField totalDisp(patch.boundaryMesh().mesh().lookupObject<volVectorField>("meshDisp"));

    word patchName = patch.name();
    label patchID = patch.boundaryMesh().findPatchID(patchName); 
    const fvPatchVectorField& totalDispPatch = totalDisp.boundaryField()[patchID]; 

    word nbrPatchName = nbrPatch.name();
    label nbrPatchID = nbrPatch.boundaryMesh().findPatchID(nbrPatchName);   
    const fvPatchVectorField& totalDispNbrPatch = totalDisp.boundaryField()[nbrPatchID]; 

       
    vectorField nf = patch.Sf() / patch.magSf();

    vectorField Cf = patch.Cf()
                   + totalDispPatch;                   
                   
    vectorField nbrCf = regionCoupledPatch_.regionCoupledPatch().interpolate
                      (
                             nbrPatch.Cf()
                           + totalDispNbrPatch
                      );

    return ((nbrCf - Cf) & nf) - offset_;
}


Foam::scalar gapContactFvPatchVectorField::boundaryStiffness() const
{
    const regionCoupledFvPatch& patch
        = refCast<const regionCoupledFvPatch>(regionCoupledPatch_);
    
    const fvPatchField<scalar>& mu =
        patch.lookupPatchField<volScalarField, scalar>("mu");


    const fvPatchField<scalar>& lambda =
        patch.lookupPatchField<volScalarField, scalar>("lambda");

    return gMin(2*(2*mu+lambda)*patch.deltaCoeffs());
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

gapContactFvPatchVectorField::gapContactFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    tractionDisplacementFvPatchVectorField(p, iF),
    regionCoupledPatch_(refCast<const regionCoupledBaseFvPatch>(p)),
    penaltyFact_(0.1),
    offset_(0.0),
    gapWidth_(p.size(), 0),
    interfaceP_(p.size(), 0)
{}


gapContactFvPatchVectorField::gapContactFvPatchVectorField
(
    const gapContactFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    tractionDisplacementFvPatchVectorField(ptf, p, iF, mapper),
    regionCoupledPatch_(refCast<const regionCoupledBaseFvPatch>(p)),
    penaltyFact_(ptf.penaltyFact_),
    offset_(ptf.offset_),
    gapWidth_(ptf.gapWidth_, mapper),
    interfaceP_(ptf.interfaceP_, mapper)
{}


gapContactFvPatchVectorField::gapContactFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    tractionDisplacementFvPatchVectorField(p, iF,dict),
    regionCoupledPatch_(refCast<const regionCoupledBaseFvPatch>(p)),
    penaltyFact_(0.1),
    offset_(0.0),
    gapWidth_(p.size(), 0),
    interfaceP_(p.size(), 0)
{
    if (regionCoupledPatch_.owner())
    {
        penaltyFact_ = dict.lookupOrDefault<scalar>("penaltyFactor", 1e-2);
        offset_ = dict.lookupOrDefault<scalar>("offset", 0.0);
    }
    
    if (dict.found("interfaceP"))
    {
        interfaceP_ = scalarField("interfaceP", dict, p.size());
    }
}


gapContactFvPatchVectorField::gapContactFvPatchVectorField
(
    const gapContactFvPatchVectorField& ptf
)
:
    tractionDisplacementFvPatchVectorField(ptf),
    regionCoupledPatch_(ptf.regionCoupledPatch_),
    penaltyFact_(ptf.penaltyFact_),
    offset_(ptf.offset_),
    gapWidth_(ptf.gapWidth_),
    interfaceP_(ptf.interfaceP_)
{}


gapContactFvPatchVectorField::gapContactFvPatchVectorField
(
    const gapContactFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    tractionDisplacementFvPatchVectorField(ptf, iF),
    regionCoupledPatch_(ptf.regionCoupledPatch_),
    penaltyFact_(ptf.penaltyFact_),
    offset_(ptf.offset_),
    gapWidth_(ptf.gapWidth_),
    interfaceP_(ptf.interfaceP_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void gapContactFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    tractionDisplacementFvPatchVectorField::autoMap(m);
    gapWidth_.autoMap(m);
    interfaceP_.autoMap(m);
}


void gapContactFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    tractionDisplacementFvPatchVectorField::rmap(ptf, addr);

    const gapContactFvPatchVectorField& dmptf =
        refCast<const gapContactFvPatchVectorField>(ptf);

    gapWidth_.rmap(dmptf.gapWidth_, addr);
    interfaceP_.rmap(dmptf.interfaceP_, addr);
}


void gapContactFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const regionCoupledFvPatch& patch
        = refCast<const regionCoupledFvPatch>(regionCoupledPatch_);
    
    const regionCoupledFvPatch& nbrPatch
        = refCast<const regionCoupledFvPatch>(patch.neighbPatch());
        
    const gapContactFvPatchVectorField& nbr = this->neighbour();
    
    if (regionCoupledPatch_.owner())
    {
        // Update the gap width
        gapWidth_ = 0.5*
        (
            gapWidth()
          + patch.regionCoupledPatch().interpolate(nbr.gapWidth())
        );
        
        // Update the neighbour gap width
        nbr.gapWidth_ = nbrPatch.regionCoupledPatch().interpolate(gapWidth_);
        
        //- Calculate the interface pressure
        interfaceP_ = max
        (
           -penaltyFact_*min(boundaryStiffness(), nbr.boundaryStiffness())*gapWidth_,
            0.0
        );
        
        nbr.interfaceP_ = nbrPatch.regionCoupledPatch().interpolate(interfaceP_);
    }
        
    // Debugging

    const_cast<fvPatchScalarField&>
    (
        patch.lookupPatchField<volScalarField, scalar>("gapWidth")
    ) = gapWidth_;
  

    //- Include the gap gas pressure
    //const gapGasModel& gapGas
    //    = patch.boundaryMesh().mesh().lookupObject<gapGasModel>("gapGas");
  /*  
    const_cast<fvPatchScalarField&>
    (
        patch.lookupPatchField<volScalarField, scalar>("interfaceP")
    ) = interfaceP_;
 */           
    pressure() = interfaceP_ ;//+ gapGas.p();

    /*
    bool isIncremental = patch.boundaryMesh().mesh().lookupObject<thermoMechanicsSolver>("thermoMechanics").isIncremental();
    if(isIncremental)
    {

            word patchName = patch.name();
            label patchID = patch.boundaryMesh().findPatchID(patchName);
            const scalarField& interfacePOld = patch.boundaryMesh().mesh().lookupObject<volScalarField>("interfaceP").oldTime().boundaryField()[patchID];         
            pressure() = pressure() - (interfacePOld + gapGas.oldTime());
    }
    */

    tractionDisplacementFvPatchVectorField::updateCoeffs();
}


void gapContactFvPatchVectorField::write(Ostream& os) const
{
    tractionDisplacementFvPatchVectorField::write(os);
    
    if (regionCoupledPatch_.owner())
    {
        //os.writeEntryIfDifferent<scalar>("penaltyFactor", 0.1, penaltyFact_);
    }
    
    interfaceP_.writeEntry("interfaceP", os);
    os.writeKeyword("penaltyFactor")
        << penaltyFact_ << token::END_STATEMENT << nl;   
    os.writeKeyword("offset")
        << offset_ << token::END_STATEMENT << nl;

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    gapContactFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
