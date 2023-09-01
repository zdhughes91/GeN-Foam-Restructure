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

#include "volFields.H"
#include "surfaceFields.H"
#include "turbulentFluidThermoModel.H"
#include "mapDistribute.H"
#include "NusseltThermalBaffle1DFvPatchScalarField.H"
#include "myOps.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

defineTypeNameAndDebug(NusseltThermalBaffle1DFvPatchScalarField, 0);
addToPatchFieldRunTimeSelection
(
    fvPatchScalarField,
    NusseltThermalBaffle1DFvPatchScalarField
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

NusseltThermalBaffle1DFvPatchScalarField::
NusseltThermalBaffle1DFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mappedPatchBase(p.patch()),
    mixedFvPatchScalarField(p, iF),
    TName_(this->internalField().name()),
    twoPhase_(myOps::split<word>(TName_, '.').size() == 2),
    twoPhaseOwner_(false),
    fluidPtr_(nullptr),
    otherFluidPtr_(nullptr),
    pairPtr_(nullptr),
    otherPairPtr_(nullptr),
    nbrPatchFieldPtr_(nullptr),
    otherFluidPatchFieldPtr_(nullptr),
    otherFluidNbrPatchFieldPtr_(nullptr),
    const_(0),
    coeff_(0),
    expRe_(0),
    expPr_(0),
    tw_(0),
    kw_(0),
    hw_(0)
{}

NusseltThermalBaffle1DFvPatchScalarField::
NusseltThermalBaffle1DFvPatchScalarField
(
    const NusseltThermalBaffle1DFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mappedPatchBase(p.patch(), ptf),
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    TName_(ptf.TName_),
    twoPhase_(ptf.twoPhase_),
    twoPhaseOwner_(ptf.twoPhaseOwner_),
    fluidPtr_(ptf.fluidPtr_),
    otherFluidPtr_(ptf.otherFluidPtr_),
    pairPtr_(ptf.pairPtr_),
    otherPairPtr_(ptf.otherPairPtr_),
    nbrPatchFieldPtr_(ptf.nbrPatchFieldPtr_),
    otherFluidPatchFieldPtr_(ptf.otherFluidPatchFieldPtr_),
    otherFluidNbrPatchFieldPtr_(ptf.otherFluidNbrPatchFieldPtr_),
    const_(ptf.const_),
    coeff_(ptf.coeff_),
    expRe_(ptf.expRe_),
    expPr_(ptf.expPr_),
    tw_(ptf.tw_),
    kw_(ptf.kw_),
    hw_(ptf.hw_)
{}

NusseltThermalBaffle1DFvPatchScalarField::
NusseltThermalBaffle1DFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mappedPatchBase(p.patch(), NEARESTPATCHFACE, dict),
    mixedFvPatchScalarField(p, iF),
    TName_(this->internalField().name()),
    twoPhase_(myOps::split<word>(TName_, '.').size() == 2),
    twoPhaseOwner_
    (
        twoPhase_
    and this->internalField().mesh().lookupClass<fluid>().size() == 1
    ),
    fluidPtr_(nullptr),
    otherFluidPtr_(nullptr),
    pairPtr_(nullptr),
    otherPairPtr_(nullptr),
    nbrPatchFieldPtr_(nullptr),
    otherFluidPatchFieldPtr_(nullptr),
    otherFluidNbrPatchFieldPtr_(nullptr),
    const_(0),
    coeff_(0),
    expRe_(0),
    expPr_(0),
    tw_(0),
    kw_(0),
    hw_(0)
{
    //- If I am either in onePhase OR (in twoPhase AND I am the controller fluid)
    if (!twoPhase_ or twoPhaseOwner_)
    {
        if (this->owner())
        {
            const_ = dict.get<scalar>("const");
            coeff_ = dict.get<scalar>("coeff");
            expRe_ = dict.get<scalar>("expRe");
            expPr_ = dict.get<scalar>("expPr");
            kw_ = dict.get<scalar>("kappa");
            tw_ = dict.get<scalar>("thickness");
            hw_ = max(kw_/max(tw_, 1e-69), 1e-69);
        }
        else
        {
            const fvPatch& nbrPatch =
                patch().boundaryMesh()[samplePolyPatch().index()];
            const NusseltThermalBaffle1DFvPatchScalarField& nbrField =
            refCast<const NusseltThermalBaffle1DFvPatchScalarField>
            (
                nbrPatch.lookupPatchField<volScalarField, scalar>
                (
                    TName_
                )
            );
            const_ = dict.getOrDefault<scalar>("const", nbrField.cnst());
            coeff_ = dict.getOrDefault<scalar>("coeff", nbrField.coeff());
            expRe_ = dict.getOrDefault<scalar>("expRe", nbrField.expRe());
            expPr_ = dict.getOrDefault<scalar>("expPr", nbrField.expPr());
            kw_ = nbrField.kw();
            tw_ = nbrField.tw();
            hw_ = nbrField.hw();
        }
    }
    else //- If I am in twoPhase AND I am not the controller fluid
    {
        if (this->owner())
        {
            //- Set ptr to otherFluidPatchFieldPtr_ but I cannot do this with
            //  setPtrs as it might cause the whole constructor to fail (
            //  setPtrs does other things as well, and does other things would
            //  end up looking for patches that might not exist yet) 
            word fluidName(myOps::split<word>(TName_, '.')[1]); 
            HashTable<const fluid*> fluids
            (
                this->internalField().mesh().lookupClass<fluid>()
            );
            word otherFluidTName
            (
                "T."
            +   (
                    (fluids[fluids.toc()[0]]->name() == fluidName) ?
                    fluids[fluids.toc()[1]]->name() : 
                    fluids[fluids.toc()[0]]->name()
                )
            );
            otherFluidPatchFieldPtr_ =
            &(  
                const_cast<NusseltThermalBaffle1DFvPatchScalarField&>
                (
                    refCast<const NusseltThermalBaffle1DFvPatchScalarField>
                    (
                        patch().lookupPatchField
                        <
                            volScalarField, 
                            scalar
                        >
                        (
                            otherFluidTName
                        )
                    )
                )
            );
            const_ = dict.getOrDefault<scalar>
                ("const", otherFluidPatchFieldPtr_->cnst());
            coeff_ = dict.getOrDefault<scalar>
                ("coeff", otherFluidPatchFieldPtr_->coeff());
            expRe_ = dict.getOrDefault<scalar>
                ("expRe", otherFluidPatchFieldPtr_->expRe());
            expPr_ = dict.getOrDefault<scalar>
                ("expPr", otherFluidPatchFieldPtr_->expPr());
            tw_ = otherFluidPatchFieldPtr_->tw();
            kw_ = otherFluidPatchFieldPtr_->kw();
            hw_ = otherFluidPatchFieldPtr_->hw();
        }
        else
        {
            const fvPatch& nbrPatch =
                patch().boundaryMesh()[samplePolyPatch().index()];
            const NusseltThermalBaffle1DFvPatchScalarField& nbrField =
            refCast<const NusseltThermalBaffle1DFvPatchScalarField>
            (
                nbrPatch.lookupPatchField<volScalarField, scalar>
                (
                    TName_
                )
            );
            const_ = dict.getOrDefault<scalar>("const", nbrField.cnst());
            coeff_ = dict.getOrDefault<scalar>("coeff", nbrField.coeff());
            expRe_ = dict.getOrDefault<scalar>("expRe", nbrField.expRe());
            expPr_ = dict.getOrDefault<scalar>("expPr", nbrField.expPr());
            kw_ = nbrField.kw();
            tw_ = nbrField.tw();
            hw_ = nbrField.hw();
        }
    }
    hw_ = max(hw_, 1e-69);
    
    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    if (dict.found("refValue"))
    {
        // Full restart
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        // Start from user entered data. Assume zeroGradient.
        refValue() = *this;
        refGrad() = 0.0;
        valueFraction() = 0.0;
    }
}

NusseltThermalBaffle1DFvPatchScalarField::
NusseltThermalBaffle1DFvPatchScalarField
(
    const NusseltThermalBaffle1DFvPatchScalarField& ptf
)
:
    mappedPatchBase(ptf.patch().patch(), ptf),
    mixedFvPatchScalarField(ptf),
    TName_(ptf.TName_),
    twoPhase_(ptf.twoPhase_),
    twoPhaseOwner_(ptf.twoPhaseOwner_),
    fluidPtr_(ptf.fluidPtr_),
    otherFluidPtr_(ptf.otherFluidPtr_),
    pairPtr_(ptf.pairPtr_),
    otherPairPtr_(ptf.otherPairPtr_),
    nbrPatchFieldPtr_(ptf.nbrPatchFieldPtr_),
    otherFluidPatchFieldPtr_(ptf.otherFluidPatchFieldPtr_),
    otherFluidNbrPatchFieldPtr_(ptf.otherFluidNbrPatchFieldPtr_),
    const_(ptf.const_),
    coeff_(ptf.coeff_),
    expRe_(ptf.expRe_),
    expPr_(ptf.expPr_),
    tw_(ptf.tw_),
    kw_(ptf.kw_),
    hw_(ptf.hw_)
{}

NusseltThermalBaffle1DFvPatchScalarField::
NusseltThermalBaffle1DFvPatchScalarField
(
    const NusseltThermalBaffle1DFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mappedPatchBase(ptf.patch().patch(), ptf),
    mixedFvPatchScalarField(ptf, iF),
    TName_(ptf.TName_),
    twoPhase_(ptf.twoPhase_),
    twoPhaseOwner_(ptf.twoPhaseOwner_),
    fluidPtr_(ptf.fluidPtr_),
    otherFluidPtr_(ptf.otherFluidPtr_),
    pairPtr_(ptf.pairPtr_),
    otherPairPtr_(ptf.otherPairPtr_),
    nbrPatchFieldPtr_(ptf.nbrPatchFieldPtr_),
    otherFluidPatchFieldPtr_(ptf.otherFluidPatchFieldPtr_),
    otherFluidNbrPatchFieldPtr_(ptf.otherFluidNbrPatchFieldPtr_),
    const_(ptf.const_),
    coeff_(ptf.coeff_),
    expRe_(ptf.expRe_),
    expPr_(ptf.expPr_),
    tw_(ptf.tw_),
    kw_(ptf.kw_),
    hw_(ptf.hw_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool NusseltThermalBaffle1DFvPatchScalarField::owner() const
{
    const label patchi = patch().index();

    const label nbrPatchi = samplePolyPatch().index();

    return (patchi < nbrPatchi);
}

void NusseltThermalBaffle1DFvPatchScalarField::setPtrs()
{    
    if (nbrPatchFieldPtr_ == nullptr)
    {
        nbrPatchFieldPtr_ =
            &(  
                const_cast<NusseltThermalBaffle1DFvPatchScalarField&>
                (
                    refCast<const NusseltThermalBaffle1DFvPatchScalarField>
                    (
                        (patch().boundaryMesh()
                        [
                            samplePolyPatch().index()
                        ]).lookupPatchField<volScalarField, scalar>
                        (
                            TName_
                        )
                    )
                )
            );
    }
    if (twoPhase_)
    {
        if (fluidPtr_ == nullptr or otherFluidPtr_ == nullptr)
        {
            word fluidName(myOps::split<word>(TName_, '.')[1]); 
            HashTable<const fluid*> fluids
            (
                this->internalField().mesh().lookupClass<fluid>()
            );

            if (fluids[fluids.toc()[0]]->name() == fluidName)
            {
                fluidPtr_ = fluids[fluids.toc()[0]];
                otherFluidPtr_ = fluids[fluids.toc()[1]];
            }
            else if (fluids[fluids.toc()[1]]->name() == fluidName)
            {
                fluidPtr_ = fluids[fluids.toc()[1]];
                otherFluidPtr_ = fluids[fluids.toc()[0]];
            }
        }
        if (pairPtr_ == nullptr or otherPairPtr_ == nullptr)
        {
            pairPtr_ = 
                &(
                    this->internalField().mesh().lookupObject<FSPair>
                    (
                        fluidPtr_->name()+".structure"
                    )
                );
            otherPairPtr_ = 
                &(
                    this->internalField().mesh().lookupObject<FSPair>
                    (
                        otherFluidPtr_->name()+".structure"
                    )
                );
        }
        if (otherFluidPatchFieldPtr_ == nullptr)
        {
            word otherFluidTName("T."+otherFluidPtr_->name());

            //- Check that the BC for the other fluid is also of type
            //  NusseltThermalBaffle1D
            const mixedFvPatchScalarField& otherFluidPatchField =
                refCast<const mixedFvPatchScalarField>
                (
                    patch().lookupPatchField<volScalarField, scalar>
                    (
                        otherFluidTName
                    )
                );
            if (otherFluidPatchField.type() != this->type())
            {
                FatalErrorInFunction()
                << "Patch field for field " << otherFluidTName
                << " should be of type " << this->type() << " instead of "
                << otherFluidPatchField.type()
                << exit(FatalError);
            }

            //- Set ptr to other fluid patch owner
            otherFluidPatchFieldPtr_ =
            &(  
                const_cast<NusseltThermalBaffle1DFvPatchScalarField&>
                (
                    refCast<const NusseltThermalBaffle1DFvPatchScalarField>
                    (
                        patch().lookupPatchField
                        <
                            volScalarField, 
                            scalar
                        >
                        (
                            otherFluidTName
                        )
                    )
                )
            );

            //- Set ptr to other fluid patch nbr
            const label otherFluidNbrPatchi = 
                otherFluidPatchFieldPtr_->samplePolyPatch().index();
            const fvPatch& nbrPatch = 
                otherFluidPatchFieldPtr_->patch().boundaryMesh()
                [
                    otherFluidNbrPatchi
                ];
            otherFluidNbrPatchFieldPtr_ =
            &(  
                const_cast<NusseltThermalBaffle1DFvPatchScalarField&>
                (
                    refCast<const NusseltThermalBaffle1DFvPatchScalarField>
                    (
                        nbrPatch.boundaryMesh()
                        [
                            samplePolyPatch().index()
                        ].lookupPatchField<volScalarField, scalar>
                        (
                            otherFluidTName
                        )
                    )
                )
            );
        }
    }
    else
    {
        if (fluidPtr_ == nullptr)
        {
            HashTable<const fluid*> fluids
            (
                this->internalField().mesh().lookupClass<fluid>()
            );
            fluidPtr_ = fluids[fluids.toc()[0]];
        }
        if (pairPtr_ == nullptr)
        {
            HashTable<const FSPair*> pairs
            (
                this->internalField().mesh().lookupClass<FSPair>()
            );
            pairPtr_ = pairs[pairs.toc()[0]];
        }
    }

    /*
    Info << "Fluid name: " << fluidPtr_->name() << endl;
    if (otherFluidPtr_ != nullptr)
    {
        Info << "Other fluid name: " << otherFluidPtr_->name() << endl;
    }
    */
}

void NusseltThermalBaffle1DFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mappedPatchBase::clearOut();

    mixedFvPatchScalarField::autoMap(m);

    /*if (this->owner())
    {
        anyScalarFieldMember_.autoMap(m);
    }*/
}

void NusseltThermalBaffle1DFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);

    /*if (this->owner())
    {
        const NusseltThermalBaffle1DFvPatchScalarField& tiptf =
            refCast<const NusseltThermalBaffle1DFvPatchScalarField>(ptf);
        anyScalarFieldMember_.rmap(tiptf.anyScalarFieldMember_, addr);
    }*/
}

//- Some super-weird shit going on here... I had to override the evaluate
//  method of mixedFvPatchScalarField in order to limit the deltaCoeffs 
//  as I was occasionally getting floating point exceptions here in MPI
//  for very, very funky decomposed domains (e.g., random decomposition
//  scheme). Even weirder though, If I store 
//  max(this->patch().deltaCoeffs(), 1e-69) in a scalarField and use
//  that one instead in the denominator, it crashes as if I am not
//  limiting it. No clue what is happening here, though results between
//  funky MPI calculations and regular single-core ones agree well, so...
void NusseltThermalBaffle1DFvPatchScalarField::evaluate
(
    const Pstream::commsTypes
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }
    scalarField::operator=
    (
        valueFraction()*refValue()
    +
        (1.0 - valueFraction())*
        (
            this->patchInternalField()
        +   refGrad()/max(this->patch().deltaCoeffs(), 1e-69)
        )
    );
  
    fvPatchScalarField::evaluate();
}

void NusseltThermalBaffle1DFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;

    setPtrs();

    if (this->owner())
    {
        const label patchi = patch().index();
        const label nbrPatchi = samplePolyPatch().index();

        if (twoPhaseOwner_)
        {
            const label otherFluidPatchi = 
                otherFluidPatchFieldPtr_->patch().index();
            const label otherFluidNbrPatchi = 
                otherFluidNbrPatchFieldPtr_->patch().index();

            //- Cache fluid thermal conducitivities (molecular ones, not
            //  effective ones as turbulence is supposed to be captured by the
            //  Nusselt correlation)
            
            
            //- Calculate heat transfer coefficients
            scalarField hf
            (
                calcH
                (
                    *pairPtr_, 
                    *this
                )()
            );
            scalarField nbrHf
            (
                calcH
                (
                    *pairPtr_, 
                    *nbrPatchFieldPtr_
                )()
            );
            scalarField otherFluidHf
            (
                calcH
                (
                    *otherPairPtr_, 
                    *otherFluidPatchFieldPtr_
                )()
            );
            scalarField otherFluidNbrHf
            (
                calcH
                (
                    *otherPairPtr_, 
                    *otherFluidNbrPatchFieldPtr_
                )()
            );

            //- Internal field temperatures
            scalarField Ti
            (
                this->patchInternalField()()
            );
            scalarField nbrTi
            (
                nbrPatchFieldPtr_->patchInternalField()()
            );
            scalarField otherFluidTi
            (
                otherFluidPatchFieldPtr_->patchInternalField()()
            );
            scalarField otherFluidNbrTi
            (
                otherFluidNbrPatchFieldPtr_->patchInternalField()()
            );

            //- Effective turbulent thermal conductivities used only to adjust
            //  the gradient to account for the fact that this BC "sees" a
            //  cell-center - face-center heat transfer coefficient of hf
            //  (and the calculated variants), while the rest of the code 
            //  "sees" a heat transfer coefficient of k*patch().deltaCoeffs(),
            //  with k being the effective fluid thermal conductivity on the
            //  patch
            scalarField kfEff(fluidPtr_->turbulence().kappaEff(patchi));
            scalarField nbrKfEff(fluidPtr_->turbulence().kappaEff(nbrPatchi));
            scalarField otherFluidKfEff
            (
                otherFluidPtr_->turbulence().kappaEff(otherFluidPatchi)()
            );
            scalarField otherFluidNbrKfEff
            (
                otherFluidPtr_->turbulence().kappaEff(otherFluidNbrPatchi)()
            );
            
            //- These calls are necessary to sync fields on nbr sides if the
            //  nbr patch does not belong to the same processor as the owner 
            //  (only relevant in MPI)
            mappedPatchBase::map().distribute(nbrKfEff);
            mappedPatchBase::map().distribute(otherFluidNbrKfEff);
            mappedPatchBase::map().distribute(nbrTi);
            mappedPatchBase::map().distribute(otherFluidNbrTi);
            mappedPatchBase::map().distribute(nbrHf);
            mappedPatchBase::map().distribute(otherFluidNbrHf);

            //- Set mixedFvPatchScalarField quantities required for BC update
            //  As this runs only on the owner, I update the nbr values as well
            //  so to have a truly implicit update
            scalarField hm(hf+otherFluidHf);
            scalarField nbrHm(nbrHf+otherFluidNbrHf);
            scalarField A(1.0+hm/hw_);
            scalarField B((hf*Ti+otherFluidHf*otherFluidTi)/hw_);
            scalarField nbrA(1.0+nbrHf/hw_);
            scalarField nbrB
            (
                (nbrHf*nbrTi+otherFluidNbrHf*otherFluidNbrTi)/hw_
            );
            scalarField oneByAnbrAminOne(max(A*nbrA-1.0, 1e-69));
            scalarField Tw((B*nbrA+nbrB)/oneByAnbrAminOne); //- New T wall
            scalarField nbrTw((nbrB*A+B)/oneByAnbrAminOne); //- New nbr T wall
            
            //this->valueFraction() = 0;
            //this->refValue() = 0;
            this->refGrad() = (hf/kfEff)*(Tw-Ti);
            //nbrPatchFieldPtr_->valueFraction() = 0;
            //nbrPatchFieldPtr_->refValue() = 0;
            nbrPatchFieldPtr_->refGrad() = (nbrHf/nbrKfEff)*(nbrTw-nbrTi);
            //otherFluidPatchFieldPtr_->valueFraction() = 0;
            //otherFluidPatchFieldPtr_->refValue() = 0;
            otherFluidPatchFieldPtr_->refGrad() = 
                (otherFluidHf/otherFluidKfEff)*(Tw-otherFluidTi);
            //otherFluidNbrPatchFieldPtr_->valueFraction() = 0;
            //otherFluidNbrPatchFieldPtr_->refValue() = 0;
            otherFluidNbrPatchFieldPtr_->refGrad() = 
                (otherFluidNbrHf/otherFluidNbrKfEff)*(nbrTw-otherFluidNbrTi);
        }
        else if (!twoPhase_)
        {
            //- Cache thermal conductivities
            scalarField kf(fluidPtr_->turbulence().kappaEff(patchi));
            scalarField nbrKf(fluidPtr_->turbulence().kappaEff(nbrPatchi));
            
            //- Calculate heat transfer coefficients
            scalarField hf
            (
                calcH
                (
                    *pairPtr_, 
                    *this
                )()
            );
            scalarField nbrHf
            (
                calcH
                (
                    *pairPtr_, 
                    *nbrPatchFieldPtr_
                )()
            );

            //- Effective turbulent thermal conductivities used only to adjust
            //  the gradient to account for the fact that this BC "sees" a
            //  cell-center - face-center heat transfer coefficient of hf
            //  (and the calculated variants), while the rest of the code 
            //  "sees" a heat transfer coefficient of k*patch().deltaCoeffs(),
            //  with k being the effective fluid thermal conductivity on the
            //  patch
            scalarField kfEff(fluidPtr_->turbulence().kappaEff(patchi));
            scalarField nbrKfEff(fluidPtr_->turbulence().kappaEff(nbrPatchi));

            //- Internal temperature fields
            scalarField Ti(patchInternalField()());
            scalarField nbrTi(nbrPatchFieldPtr_->patchInternalField()());
            
            //- Distribute nbr quantities to avoid BS in MPI
            mappedPatchBase::map().distribute(nbrKfEff);
            mappedPatchBase::map().distribute(nbrTi);
            mappedPatchBase::map().distribute(nbrHf);

            //- Set mixedFvPatchScalarField quantities required for BC update
            //  As this runs only on the owner, I update the nbr values as well
            //  so to have a truly implicit update
            scalarField A(1.0+hf/hw_);
            scalarField B((hf/hw_)*Ti);
            scalarField nbrA(1.0+nbrHf/hw_);
            scalarField nbrB((nbrHf/hw_)*nbrTi);
            scalarField oneByAnbrAminOne(max(A*nbrA-1.0, 1e-69));
            scalarField Tw((B*nbrA+nbrB)/oneByAnbrAminOne);
            scalarField nbrTw((nbrB*A+B)/oneByAnbrAminOne);

            //this->valueFraction() = 0;
            //this->refValue() = 0;
            this->refGrad() = (hf/kfEff)*(Tw-Ti);
            //nbrPatchFieldPtr_->valueFraction() = 0;
            //nbrPatchFieldPtr_->refValue() = 0;
            nbrPatchFieldPtr_->refGrad() = (nbrHf/nbrKfEff)*(nbrTw-nbrTi);
        }
    }

    if (debug)
    {
        const phaseCompressibleTurbulenceModel& turbModel = 
            fluidPtr_->turbulence();
        scalarField kf(turbModel.kappaEff(patch().index()));
        scalar Q = gAverage(kf*snGrad());
            Pout<< patch().boundaryMesh().mesh().name() << ':'
                << patch().name() << ':'
                << this->internalField().name() << " <- "
                << patch().boundaryMesh()[samplePolyPatch().index()].name() 
                << ':' << this->internalField().name() << " :"
                << " heat[W]:" << Q
                << " walltemperature "
                << " min:" << gMin(*this)
                << " max:" << gMax(*this)
                << " avg:" << gAverage(*this)
                << endl;
    }

    // Restore tag
    UPstream::msgType() = oldTag;

    mixedFvPatchScalarField::updateCoeffs();
}

tmp<scalarField> NusseltThermalBaffle1DFvPatchScalarField::calcH
(
    const FSPair& pair,
    const NusseltThermalBaffle1DFvPatchScalarField& patchField            
)
{
    const label patchi(patchField.patch().index());
    //- Molecular thermal conducitivity for computing the htc from the Nu
    scalarField k(pair.fluidRef().thermo().kappa(patchi));
    const scalarField& Re(pair.Re().boundaryField()[patchi]);
    const scalarField& Pr(pair.fluidRef().Pr().boundaryField()[patchi]);
    const scalarField& Dh(pair.structureRef().Dh().boundaryField()[patchi]);
    if (pair.fPtr().valid()) //- Not valid in onePhase
    {
        const scalarField& f(pair.f().boundaryField()[patchi]);
        return
            (
                (
                    patchField.cnst() 
                +   patchField.coeff()*
                    pow(Re, patchField.expRe())*
                    pow(Pr, patchField.expPr())
                )*f*k/Dh
            );
    }
    else
    {
        return
            (
                (
                    patchField.cnst() 
                +   patchField.coeff()*
                    pow(Re, patchField.expRe())*
                    pow(Pr, patchField.expPr())
                )*k/Dh
            );
    }
}

void NusseltThermalBaffle1DFvPatchScalarField::write(Ostream& os) const
{
    mixedFvPatchScalarField::write(os);
    mappedPatchBase::write(os);

    os.writeKeyword("const")<< const_
        << token::END_STATEMENT << nl;
    os.writeKeyword("coeff")<< coeff_
        << token::END_STATEMENT << nl;
    os.writeKeyword("expRe")<< expRe_
        << token::END_STATEMENT << nl;
    os.writeKeyword("expPr")<< expPr_
        << token::END_STATEMENT << nl;

    if (this->owner() and ((twoPhase_ and twoPhaseOwner_) or !twoPhase_))
    {
        os.writeKeyword("thickness")<< tw_
            << token::END_STATEMENT << nl;
        os.writeKeyword("kappa")<< kw_
            << token::END_STATEMENT << nl;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
