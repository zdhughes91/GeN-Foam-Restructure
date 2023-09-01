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

#include "onePhaseLegacy.H"
#include "addToRunTimeSelectionTable.H"
#include "regimeMapModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace thermalHydraulicsModels
{
    defineTypeNameAndDebug(onePhaseLegacy, 0);
    addToRunTimeSelectionTable
    (
        thermalHydraulicsModel,
        onePhaseLegacy,
        thermalHydraulicsModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermalHydraulicsModels::onePhaseLegacy::onePhaseLegacy
(
    Time& time,
    fvMesh& mesh,
    customPimpleControl& pimple,
    fv::options& fvOptions
)
:
    thermalHydraulicsModel
    (
        time,
        mesh,
        pimple,
        fvOptions
    ),
    //- The structure model needs to be created before the fluids because
    //  the turbulence models created by the fluids might require a reference
    //  to the structure (obtained via objectRegistry lookup in the specific
    //  turbulence model class)
    structure_
    (
        this->subDict("structureProperties"),
        mesh,
        this->powerDensityNeutronics_
    ),
    fluid_
    (
        (this->found("fluidProperties"))
    ?   this->subDict("fluidProperties") : *this,
        mesh,
        word(""),   //- This is the phase name, setting it to "" signals a
                    //  onePhase solver to the rest of the FFS library
        this->powerDensityNeutronicsToLiquid_,
        false       //- No need to read or write the fluid phaseFraction in
                    //  onePhase, it is tied to the structure phaseFraction
    ),
    FSPair_(fluid_, structure_, *this),
    fixedRho_(fluid_.thermo().rho()),
    rhok_
    (
        volScalarField
        (
            IOobject
            (
                "rhok",
                mesh_.time().timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("", dimless, 1.0),
            zeroGradientFvPatchScalarField::typeName
        )
    )
{
    fixedRho_.correctBoundaryConditions();
    rhok_.correctBoundaryConditions();

    //- checking if the incompressible treatment is activated. In that case
    //- the solver becomes a porous version of buoyantBoussinesqPIMPLFoam.
    //- If this is activated, the only equation of state accepted is
    //- rhoConst.
    dictionary phaseDict = fluid_.dict();
    incompressibleTreatment_ = bool
    (
    	phaseDict.lookup
    	(
        	"incompressibleBuoyantBoussinesqTreatment"
    	)
    );

    //- Checking if heRhoThermo has the correct equationOfState for the
    //- current solver.
    dictionary thermoDict = fluid_.thermo().subDict("thermoType");
    if(incompressibleTreatment_
       and word(thermoDict.lookup("equationOfState")) != "rhoConst")
    {
        Foam::error e("The equation of state is not rhoConst (constant density)! "
                      "If 'incompressibleBuoyantBoussinesqTreatment' is activated,"
                      "one is only allowed to use rhoConst!");
        e.exit(100);
    }

    //- Reading reference temperature and thermal expansion coefficient
    //- for the Boussinesq approximation
    dictionary eosDict = fluid_.thermo()
                         .subDict("mixture")
                         .subDict("equationOfState");

    beta_.set
    (
        new dimensionedScalar
        (
            "beta",
            pow(dimTemperature,-1),
            eosDict.lookupOrDefault<scalar>("beta", 0.0)
        )
    );

    Tref_.set
    (
        new dimensionedScalar
        (
            "Tref",
            dimTemperature,
            eosDict.lookupOrDefault<scalar>("T0", 0.0)
        )
    );

    //- Create turbulence model
    fluid_.constructTurbulenceModel();

    //- Set phase fraction fields (constant in time), structure has priority
    fluid_.volScalarField::operator=(1.0-structure_);

    //- The normalized field is non-trivial (i.e. different than 1) only in the
    //  twoPhase solver. However, it is used by some models in the shared
    //  thermal-hydraulics library, so it should be set nonetheless! The most
    //  important quantity that relies on this is the Reynolds computed by
    //  the FSPair object
    fluid_.normalized() = fluid_/(1.0-structure_);

    //- Calculating rhok value for boussinesq approximation if incompressible flow.
    //- also, update thermo.tho() to make sure the neutronics solver has access to
    //- the density feedback.
    if(incompressibleTreatment_)
    {
        rhok_ = 1.0 - beta_()*(fluid_.thermo().T() - Tref_());
        rhok_.correctBoundaryConditions();
        fluid_.thermo().rho() = fixedRho_*rhok_;
        fluid_.thermo().rho().correctBoundaryConditions();
    }

    //- Set fluid characteristic dimension to structure hydraulic diameter.
    //  This is handled by the fluidGeometry class in the twoPhase solver
    //  (as the fluid characteristic dimension will depend on the regime)
    //  so here I need to do it manually. Since I assume that the structure is
    //  immutable, this is done only once
    fluid_.Dh() = structure_.Dh();

    //- Initialize fluid-intensive fluxes (i.e. that depend on the phase
    //  fraction, namely alphaPhi and alphaRhoPhi, which are the REAL
    //  volumetric flux in m3/s and the REAL mass flux in kg/s. By REAL I mean
    //  not superficial). This is done after the phaseFraction normalization
    //  step to ensure consistency. This step has an effect ONLY IF the
    //  alphaPhi, alphaRhoPhi fields were NOT found on disk
    fluid_.initAlphaPhis();

    //- The total volumetric flux is the REAL fluid volumetric flux. Might
    //  as well remove phi_ entirely as a field, I know... Maybe in the future
    phi_ = fluid_.alphaPhi();

    //- Compute initialFluidMass
    initialFluidMass_ = fvc::domainIntegrate(fluid_.thermo().rho()*fluid_);

    //- Setting up the initial Darcy velocity and flux
    UDarcy_.set
    (
        new volVectorField
        (
            IOobject
            (
                "UDarcy",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            fluid_.U()*fluid_
        )
    );

    //- Making sure that the division of boundary conditions does not result in
    //- a "calculated" boundary field. So we set the same boundary conditions as
    //- the ones given for the real velocity. Exceptions are wedge and empty BCs,
    //- because those are not touched.

    forAll(fluid_.U().boundaryField(), bcInd)
    {
        if(fluid_.U().boundaryField()[bcInd].type() == "empty" or
           fluid_.U().boundaryField()[bcInd].type() == "wedge")
        {
            Info << "Skipping boundary type assignment for UDarcy"
                 << "because it is empty/wedge." << endl;
            continue;
        }

        word bcType = fluid_.U().boundaryField()[bcInd].type();

        tmp<fvPatchField<vector>> originalPatch(UDarcy_().boundaryField()[bcInd]);
        
        UDarcy_().boundaryFieldRef().set
        (
            bcInd,
            fvPatchField<vector>::New
            (
                bcType,
                UDarcy_().mesh().boundary()[bcInd],
                UDarcy_()
            )
        );

        forAll(UDarcy_().boundaryFieldRef()[bcInd], faceI)
        {
            UDarcy_().boundaryFieldRef()[bcInd][faceI] = originalPatch()[faceI];
        }
    }
    UDarcy_().write();

    phiDarcy_.set
    (
        new surfaceScalarField
        (
            IOobject
            (
                "phiDarcy",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fvc::interpolate(fixedRho_)*fvc::interpolate(fluid_)*fluid_.phi()
        )
    );

    //- Initialize continuity errors. It's important to do it here or, if not
    //  solving for fluidMechanics, these would never get corrected
    correctContErr();

    Info << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

//- Solve according to flags
void Foam::thermalHydraulicsModels::onePhaseLegacy::correct
(
    scalar& residual,
    bool solveFluidDynamics,
    bool solveEnergy
)
{
    correctModels(solveFluidDynamics, solveEnergy);
    if (solveFluidDynamics)
    {
        correctFluidMechanics(residual);
    }

    if (solveEnergy)
    {
        correctEnergy(residual);
    }
    Info << endl;
}

void Foam::thermalHydraulicsModels::onePhaseLegacy::correctFluidMechanics
(
    scalar& residual
)
{
    #include "UEqn_1pl.H"
    #include "pEqn_1pl.H"
}

void Foam::thermalHydraulicsModels::onePhaseLegacy::correctEnergy(scalar& residual)
{
    #include "EEqn_1pl.H"
}

void Foam::thermalHydraulicsModels::onePhaseLegacy::correctModels
(
    bool solveFluidDynamics,
    bool solveEnergy
)
{
    this->correctRegimeMaps();
    FSPair_.correct(solveFluidDynamics, solveEnergy);
}

void Foam::thermalHydraulicsModels::onePhaseLegacy::correctCourant()
{
    CoNum_ = 0.0;
    meanCoNum_ = 0.0;

    scalarField sumPhi
    (
        fvc::surfaceSum
        (
            mag
            (
                fvc::interpolate(UDarcy_()) & mesh_.Sf()
            )
        )().primitiveField()/fluid_.primitiveField()
    );

    CoNum_ = 0.5*gMax(sumPhi/mesh_.V().field())*runTime_.deltaTValue();

    meanCoNum_ =
        0.5*(gSum(sumPhi)/gSum(mesh_.V().field()))*runTime_.deltaTValue();

    Info<< "Courant Number (avg max) = " << meanCoNum_
        << " " << CoNum_ << endl;
}

void Foam::thermalHydraulicsModels::onePhaseLegacy::correctContErr()
{
	if(incompressibleTreatment_)
	{
		fluid_.contErr() =
    	(
    	    fvc::div(phiDarcy_())
    	);
    	fluid_.contErr().correctBoundaryConditions();
	}
	else
	{
		volScalarField& rho(fluid_.thermo().rho());
		fluid_.contErr() =
    	(
    	        fvc::ddt(fluid_, rho)
    	    +   fvc::div(phiDarcy_())
    	    -   (fvOptions_(fluid_, fixedRho_) & rho)
    	);
    	fluid_.contErr().correctBoundaryConditions();
	}

}

// ************************************************************************* //
