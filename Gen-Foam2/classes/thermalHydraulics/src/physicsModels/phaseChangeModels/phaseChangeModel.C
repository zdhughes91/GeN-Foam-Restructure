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

#include "phaseChangeModel.H"
#include "FFPair.H"
#include "fvCFD.H"
#include "myOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(phaseChangeModel, 0);
    defineRunTimeSelectionTable(phaseChangeModel, phaseChangeModels);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseChangeModel::phaseChangeModel
(
    FFPair& pair,
    const dictionary& dict
)
:
    IOdictionary
    (
        IOobject
        (
            IOobject::groupName("phaseChangeModel", typeName),
            pair.mesh().time().constant(),
            pair.mesh()
        ),
        dict
    ),
    mesh_(pair.mesh()),
    pimple_(pair.pimple()),
    pair_(pair),
    fluid1_(pair.fluid1()),
    fluid2_(pair.fluid2()),
    liquid_
    (
        (pair.fluid1().isLiquid()) ? 
        pair.fluid1() : pair.fluid2()
    ),
    vapour_
    (
        (pair.fluid1().isGas()) ? 
        pair.fluid1() : pair.fluid2()
    ),
    p_(mesh_.lookupObject<volScalarField>("p")),
    T1_(fluid1_.thermo().T()),
    T2_(fluid2_.thermo().T()),
    htc1_(pair.htc(fluid1_.name())),
    htc2_(pair.htc(fluid2_.name())),
    iT_(pair.iT()),
    iA_(pair.iA()),
    L_
    (
        IOobject
        (
            IOobject::groupName("L", pair.name()),
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimEnergy/dimMass, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    dmdtI_
    (
        IOobject
        (
            IOobject::groupName("dmdtI", pair.name()),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimDensity/dimTime, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    dmdtW_
    (
        IOobject
        (
            IOobject::groupName("dmdtW", pair.name()),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimDensity/dimTime, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    dmdt_
    (
        IOobject
        (
            IOobject::groupName("dmdt", pair.name()),
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimDensity/dimTime, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    latentHeatPtr_
    (
        latentHeatModel::New
        (
            *this,
            this->subDict("latentHeatModel"),
            mesh_
        )
    ),
    saturationPtr_
    (
        saturationModel::New
        (
            *this,
            this->subDict("saturationModel"),
            mesh_
        )
    ),
    residualIA_
    (
        dict.getOrDefault<scalar>("residualInterfacialArea", 1e-3)
    ),
    residualIACells_(0),
    dmLostToLimiter_(0)
{
    //- 
    if 
    (
        !(fluid1_.isLiquid() and fluid2_.isGas())
    and !(fluid2_.isLiquid() and fluid1_.isGas())
    )
    {
        FatalErrorInFunction
            << "Either phase " << fluid1_.name() << " or " << fluid2_.name()
            << " have an undetermined stateOfMatter (should be specified in "
            << "phaseProperties." << fluid1_.name() << "Properties and/or "
            << "phaseProperties." << fluid2_.name() << "Properties)"
            << exit(FatalError);
    }

    //- Set initial interfacial temperature (to avoid problems with the
    //  under-relaxation at the first time-step in EEqns.H
    IOobject iTHeader
    (
        iT_.name(),
        mesh_.time().timeName(),
        mesh_,
        IOobject::NO_READ
    );
    if (!iTHeader.typeHeaderOk<volScalarField>(true))
    {
        saturationPtr_->correctField(iT_, "TSat");
    }

    //-
    latentHeatPtr_->correctField(L_);

    //- Read residualIACells from cellZones
    wordList residualIARegions
    (
        this->lookupOrDefault<wordList>
        (
            "residualInterfacialAreaRegions", 
            wordList()
        )
    );
    forAll(residualIARegions, i)
    {
        const labelList& regionCells
        (
            mesh_.cellZones()[residualIARegions[i]]
        );
        forAll(regionCells, j)
        {
            residualIACells_.append(regionCells[j]);
        }
    }

    //- It is dimPower instead of dimPower/dimVolume as the fvScalarMatrix
    //  expects objects that are either already integrated over the volume
    //  (e.g. what the fvm::Su, SuSp, etc. return) or objects that are not
    //  already volume integrated, yet whose dimensions are equal to the
    //  dimension of the fvScalarMatrix (dimPower in this case) divided by
    //  dimVolume (e.g. explicit terms)
    heSources_.set
    (
        fluid1_.thermo().he().name(), 
        new fvScalarMatrix(fluid1_.thermo().he(), dimPower)
    );
    heSources_.set
    (
        fluid2_.thermo().he().name(), 
        new fvScalarMatrix(fluid2_.thermo().he(), dimPower)
    );
}

// * * * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * //

void Foam::phaseChangeModel::limitInterfacialArea()
{
    if (residualIA_ != 0.0)
    {
        if (residualIACells_.size() > 0)
        {
            forAll(residualIACells_, i)
            {
                const label& celli(residualIACells_[i]);
                scalar& iAi(iA_[celli]);
                iAi = max(iAi, residualIA_);
            }
        }
        else
        {
            forAll(mesh_.cells(), celli)
            {
                scalar& iAi(iA_[celli]);
                iAi = max(iAi, residualIA_);
            }
        }
    }
}

void Foam::phaseChangeModel::limitMassTransfer()
{
    //- Prevent removing mass from phases that are not present in a cell
    dmdt_ = posPart(dmdt_)*pos(fluid1_) + negPart(dmdt_)*pos(fluid2_);

    /*
    Given that the continuity equations are:

    ddt(alpha1*rho1) + div(alphaRhoPhi1) = -dmdt
    ddt(alpha2*rho2) + div(alphaRhoPhi2) = dmdt

    and that they are solved explicitly (albeit with the whole MULES limiter
    step before the explicit solution), there should be a maximum allowable
    limit on the value of the mass transfer term. This can be computed from 
    the continuity equations assuming the maximum allowable ddt within a 
    single time-step (i.e. that ddt that would make my phase fraction change
    from its current value to 0 in one time step). Due to the way the 
    equations are implemented, a maximum on dmdt is derived from the 
    continuity equation of phase 1 while a minimum is derived from that of
    phase 2. Continuity errors are (for now) not considered.
    */
    if (pimple_.dict().found("massTransferSafetyFactor"))
    {
        scalar SF(pimple_.dict().get<scalar>("massTransferSafetyFactor"));

        scalar dt(mesh_.time().deltaT().value());
        const scalarField& V(mesh_.V());
        scalar totV(0);
        forAll(V, i)
        {
            totV += V[i];
        }

        scalarField div1 = fvc::div(fluid1_.alphaRhoPhi())().primitiveField();
        scalarField div2 = fvc::div(fluid2_.alphaRhoPhi())().primitiveField();

        forAll(dmdt_, i)
        {
            scalar& dmdt(dmdt_[i]);
            if (dmdt > 0.0)
            {
                const scalar& Vi(V[i]);
                scalar maxDmdt
                (
                    max(SF*fluid1_[i]*fluid1_.rho()[i]/dt-div1[i], 0.0)
                );
                scalar dmdt0(dmdt);
                dmdt = min(dmdt, maxDmdt);
                if (pimple_.finalIter())
                    dmLostToLimiter_ += (dmdt0-dmdt)*Vi*dt;
            }
            else if (dmdt < 0.0)
            {
                const scalar& Vi(V[i]);
                scalar minDmdt
                (
                    min(-SF*fluid2_[i]*fluid2_.rho()[i]/dt+div2[i], 0.0)
                );
                scalar dmdt0(dmdt);
                dmdt = max(dmdt, minDmdt);
                if (pimple_.finalIter())
                    dmLostToLimiter_ += (dmdt0-dmdt)*Vi*dt;
            }
        }

        /*if (pimple_.finalIter())
        {
            Info<< "Cumulative dm lost to limiter = " 
                << (dmLostToLimiter_) 
                << " kg/m3" << endl;
        }*/
    }

    dmdt_.correctBoundaryConditions();
}

void Foam::phaseChangeModel::correctInterfacialTemperature()
{
    if (mesh_.relaxField("T.interface"))
        iT_.storePrevIter();
    if (pimple_.dict().found("maxTInterfaceDdt"))
    {
        scalar maxDTdt
        (
            pimple_.dict().get<scalar>
            (
                "maxTInterfaceDdt"
            )
        );
        
        scalar dt(mesh_.time().deltaTValue());
        saturationPtr_->correctField(iT_, "TSat");
        
        forAll(iT_, i)
        {
            const scalar& iT0(iT_.oldTime()[i]);
            scalar& iT(iT_[i]);
            scalar dTdt(min(max((iT-iT0)/dt, -maxDTdt), maxDTdt));
            iT = iT0 + dTdt*dt;
        }
    }
    else
    {
        saturationPtr_->correctField(iT_, "TSat");
    }

    iT_.relax();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::phaseChangeModel::correct()
{
    //- Limit interfacial area so boiling can start
    //  (very crude, it's the best I have for now)
    limitInterfacialArea();

    //- Update saturation temperature
    correctInterfacialTemperature();

    //- Update latent heat
    myOps::storePrevIterIfRelax(L_);
    latentHeatPtr_->correctField(L_);
    L_.relax();

    //- Calculate dmdtI <- depends on the actual run-time selected phaseChange 
    //  model
    correctInterfacialDmdt();

    //- Update total dmdt and relax
    dmdt_.storePrevIter();
    dmdt_ = dmdtI_ + dmdtW_;
    limitMassTransfer(); 
    dmdt_.relax();

    //- The rest of this function is to set heSources_ in an 
    //  energy-conservative way
    //- Refs and fields
    const volScalarField& he1(fluid1_.thermo().he());
    const volScalarField& he2(fluid2_.thermo().he());
    const volScalarField& T1(fluid1_.thermo().T());
    const volScalarField& T2(fluid2_.thermo().T());
    volScalarField& htc1(pair_.htc(fluid1_.name()));
    volScalarField& htc2(pair_.htc(fluid2_.name()));
    volScalarField he1I(fluid1_.thermo().he(p_, iT_));
    volScalarField he2I(fluid2_.thermo().he(p_, iT_));
    
    /*-----------------------------------------------------------------------*/

    volScalarField c1
    (
        IOobject
        (
            "heDmdtSourceCoeff."+fluid1_.name(),
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("", dimDensity/dimTime, 0),
        zeroGradientFvPatchScalarField::typeName
    );
    volScalarField c2
    (
        IOobject
        (
            "heDmdtSourceCoeff."+fluid2_.name(),
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("", dimDensity/dimTime, 0),
        zeroGradientFvPatchScalarField::typeName
    );
    
    //- Fraction of total dmdt that is due to wall boiling (can be negative,
    //  e.g. vapour at saturation but sub-cooled boiling implies that dmdtI
    //  and dmdt (=dmdtI+dmdtW) have opposing signs)
    scalarField fw(mesh_.cells().size(), 0.0);
    scalar dmdtMin(1e-2);
    forAll(mesh_.cells(), i)
    {
        const scalar& dmdti(dmdt_[i]);
        if (mag(dmdti) > dmdtMin)
        {
            fw[i] = dmdtW_[i]/dmdti;
        }
    }

    //- The htcs are set to 0 in the phase-change region and the mass 
    //  transfer enthalpy contribution are accounted for via a SuSp(c, he) term
    //  for both phases. The next for loop is used to calculate these terms.
    //  Nonetheless, for visual purposes, I still want the htcs to have their
    //  model-computed values in the boiling/condensing region just so that I
    //  can have an idea of what is going on in paraView. This storePrevIter is
    //  thus used for chaching purposes (the htcs are restored only in FFPair.C
    //  correct() after the interfacial non-mass-transfer 
    //  enthalpy contributions are added).
    htc1.storePrevIter();
    htc2.storePrevIter();
    forAll(mesh_.cells(), i)
    {
        const scalar& dmdti(dmdt_[i]);
        if (mag(dmdti) > dmdtMin)
        {
            const scalar& he1i(he1[i]);
            const scalar& he2i(he2[i]);
            const scalar& iTi(iT_[i]);
            const scalar& iAi(iA_[i]);
            const scalar& Li(L_[i]);
            const scalar& fwi(fw[i]);
            scalar& htc1i(htc1[i]);
            scalar& htc2i(htc2[i]);
            scalar& c1i(c1[i]);
            scalar& c2i(c2[i]);
            scalar dmdtIi((1.0-fwi)*dmdti);
            scalar dmdtWi(fwi*dmdti);
            scalar q1i(iAi*htc1i*(T1[i]-iTi));
            scalar q2i(iAi*htc2i*(T2[i]-iTi));

            //- Interfacial contribution
            /*
                Nothing guarantees that the computed dmdtI
                already satisfies energy conservation in the sense that
                dmdtIi = (q1i+q2i)/Li, as it might have: 1) been computed
                via other approaches (e.g. gas kinetic theory); 2) been
                under-relaxed after its calculation, even if it was done
                via energy conservative approaches (e.g. heat conduction
                limited). Thus, to be consistent, I need to adjust the
                interfacial heat fluxes by the ratio of the original
                dmdt had it been computed via a perfectly 
                energy-conservative approach (i.e. dmdtIi0) to the actual
                dmdti which is not assured to satisfy dmdtiI = (q1i+q2i)/Li.
                In this way I assure energy conservation when accounting 
                for dmdti in the energy equations, regardless of how dmdti
                was actually computed (note, there is nothing un-physical
                about computing dmdt with approaches that are not based
                on energy conservation, as long as the dmdtI is added to
                the energy equations in a way that is energy conservative).
                However, I am not sure that this energy-conservativity 
                necessarily lead to meaningful temperature profiles IF
                the dmdtI is not computed via approaches other than the
                heat-conduction limited one. Oh, well!
            */
            scalar dmdtIi0((q1i+q2i)/Li);
            dmdtIi0 = 
                (dmdtIi0 >= 0.0) ? 
                max(dmdtMin, dmdtIi0) : 
                min(-dmdtMin, dmdtIi0);
            scalar f(dmdtIi/dmdtIi0);
            scalar dmdtILi(dmdtIi*Li);
            c1i = (dmdtILi-f*q2i)/he1i;
            c2i = (dmdtILi-f*q1i)/he2i;

            //- Wall contribution
            if (dmdtWi > 0.0)
                c1i += (dmdtWi*Li)/he1i;
            else if (dmdtWi < 0.0)
                c2i += (dmdtWi*Li)/he2i;

            //- Set htcs to 0 so that the interfacial contribtuion is not
            //  counted twice (added in FFPair.C) in cells that are boiling
            htc1i = 0;
            htc2i = 0;
        }
    }

    /*-----------------------------------------------------------------------*/

    //- This is longer to write but faster than using posPart/negPart
    //  because I do not care about boundary values
    volScalarField dmdt12
    (
        IOobject
        (
            "posPartDmdt."+pair_.name(),
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("", dimDensity/dimTime, 0),
        zeroGradientFvPatchScalarField::typeName
    );
    volScalarField dmdt21
    (
        IOobject
        (
            "negPartDmdt."+pair_.name(),
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("", dimDensity/dimTime, 0),
        zeroGradientFvPatchScalarField::typeName
    );
    forAll(mesh_.cells(), i)
    {
        const scalar& dmdti(dmdt_[i]);
        if (dmdti > 0.0)
            dmdt12[i] = dmdti;
        else if (dmdti < 0.0)
            dmdt21[i] = -dmdti;
    }
    dmdt12.correctBoundaryConditions();
    dmdt21.correctBoundaryConditions();

    //- If phase 1 gains enthalpy from the disappearance of phase 2 (e.g. if 
    //  phase 1 is liquid, 2 is vapour and the vapour is condensing) then the
    //  enthalpy added to phase 1 is at saturation (i.e. the interfacial 
    //  enthalpy, he1I). If phase 1 loses enthalpy from its own disappearance 
    //  (e.g. if phase 1 is liquid, 2 is vapour and the liquid is boiling) then
    //  the enthalpy is removed from phase 1 at its current enthalpy (he1).
    //  This strategy is to prevent possible thermal run-aways. Note that 
    //  beacuse of this formulation, the enthalpy gain term can be treated
    //  implicitly but the loss term must be treated explcitly. The same logic
    //  is applied to the he source of phase 2. It works regardless of which
    //  phase is liquid and which phase is vapour
    *(heSources_[he1.name()]) = 
        fvm::SuSp(c1, he1)      //- Contrib. from interf. and wall mass transf.
    +   fvm::Sp(dmdt12, he1)    //- Contrib. from intrinsic phase change (imp.)
    -   dmdt21*he1I;            //- Contrib. from intrinsic phase change (exp.)
    *(heSources_[he2.name()]) = 
        fvm::SuSp(c2, he2)
    +   fvm::Sp(dmdt21, he2)
    -   dmdt12*he2I;
    
    //- Store prev iters for under-relaxation (at the next iteration)
    dmdtI_.storePrevIter();
    dmdtW_.storePrevIter();

    //- Reset wall contribution (it is computed by some 
    //  FSHeatTransferCoefficient models but due to how said models work, it
    //  cannot be reset from within them)
    forAll(mesh_.cells(), i)
    {
        dmdtW_[i] = 0;
    }
    dmdtW_.correctBoundaryConditions();

    /*
    forAll(mesh_.cells(), i)
    {
        scalar deltai
        (
            he1I[i]+L_[i]-he2I[i]
        );
        
        scalar Xi
        (
            (
                (
                    fluid1_.normalized()[i]*fluid1_.rho()[i]*he1[i]
                +   fluid2_.normalized()[i]*fluid2_.rho()[i]*(he2[i]+deltai)
                )/
                (
                    fluid1_.normalized()[i]*fluid1_.rho()[i]
                +   fluid2_.normalized()[i]*fluid2_.rho()[i]
                ) - he1I[i]
            )/(he2I[i]+deltai-he1I[i])
        );
        Info<< he1[i] << " " << he1I[i] << " " << (he2[i]+deltai) << " " 
            << (he2I[i]+deltai) << " " << L_[i] << " | X = " << Xi << endl;
    }*/
}

// ************************************************************************* //
