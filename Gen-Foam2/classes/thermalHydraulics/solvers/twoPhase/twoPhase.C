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

#include "twoPhase.H"
#include <chrono> 
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace thermalHydraulicsModels
{
    defineTypeNameAndDebug(twoPhase, 0);
    addToRunTimeSelectionTable
    (
        thermalHydraulicsModel, 
        twoPhase, 
        thermalHydraulicsModels
    );
}
}

const Foam::Enum
<
    Foam::thermalHydraulicsModels::twoPhase::alphaEqnsSolver
>
Foam::thermalHydraulicsModels::twoPhase::alphaEqnsSolverNames_
(
    {
        { 
            alphaEqnsSolver::MULES, 
            "MULES" 
        },
        { 
            alphaEqnsSolver::implicitUpwind, 
            "implicitUpwind" 
        }
    }
);

const Foam::Enum
<
    Foam::thermalHydraulicsModels::twoPhase::partialEliminationMode
>
Foam::thermalHydraulicsModels::twoPhase::partialEliminationModeNames_
(
    {
        { 
            partialEliminationMode::none, 
            "none" 
        },
        { 
            partialEliminationMode::legacy, 
            "explicit" 
        },
        { 
            partialEliminationMode::implicit, 
            "implicit" 
        },
        {
            partialEliminationMode::implicitWithDmdt, 
            "implicitWithDmdt"
        }
    }
);

const Foam::Enum
<
    Foam::thermalHydraulicsModels::twoPhase::contErrCompensationMode
>
Foam::thermalHydraulicsModels::twoPhase::contErrCompensationModeNames_
(
    {
        { 
            contErrCompensationMode::none, 
            "none" 
        },
        { 
            contErrCompensationMode::Su, 
            "Su" 
        },
        { 
            contErrCompensationMode::Sp, 
            "Sp" 
        },
        { 
            contErrCompensationMode::SuSp, 
            "SuSp" 
        }
    }
);

const Foam::Enum
<
    Foam::thermalHydraulicsModels::twoPhase::heStabilizationMode
>
Foam::thermalHydraulicsModels::twoPhase::heStabilizationModeNames_
(
    {
        { 
            heStabilizationMode::cutoff, 
            "cutoff" 
        },
        { 
            heStabilizationMode::source, 
            "source" 
        }
    }
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermalHydraulicsModels::twoPhase::twoPhase
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
    //- The structure needs to be created before the fluids because
    //  the turbulence models created by the fluids might require a reference
    //  to the structure (obtained via objectRegistry lookup in the specific
    //  turbulence model class)
    structure_
    (
        this->subDict("structureProperties"),
        mesh,
        this->powerDensityNeutronics_ 
    ),
    fluid1_
    (
        this->subDict
        (
                word(this->lookup("fluid1"))
            +   "Properties"
        ),
        mesh,
        word(this->lookup("fluid1")),
        this->powerDensityNeutronicsToLiquid_
    ),
    fluid2_
    (
        this->subDict
        (
                word(this->lookup("fluid2"))
            +   "Properties"
        ),
        mesh,
        word(this->lookup("fluid2")),
        this->powerDensityNeutronicsToLiquid_
    ),
    movingAlpha_
    (
        "movingAlpha",
        1.0 - structure_
    ),
    FFPair_(fluid1_, fluid2_, *this),
    F1SPair_(fluid1_, structure_, *this),
    F2SPair_(fluid2_, structure_, *this),
    U_
    (
        IOobject
        (
            "U",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("", dimVelocity, vector::zero)
    ),
    rho_
    (
        IOobject
        (
            "rho",  // <- Cannot be named anything other than "rho" as it is
                    // looked-up by methods internal to some standard OpenFOAM
                    // ones I use. I recall this is the mixture density
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimDensity, 0.0)
    ),
    bothPhasesArePresent_(false),
    withinMarginToPhaseChange_(false),
    alphaEqnsSolver_
    (
        alphaEqnsSolverNames_.get
        (
            mesh_.solverDict("alpha").get<word>("solver")
        )
    ),
    partialEliminationMode_
    (
        partialEliminationModeNames_.get
        (
            pimple_.dict().lookupOrDefault<word>
            (
                "partialEliminationMode", 
                "none"
            )
        )
    ),
    UContErrCompensationMode_
    (
        contErrCompensationModeNames_.get
        (
            (pimple_.dict().found("UContinuityErrorCompensationMode")) ?
                pimple_.dict().get<word>
                (
                    "UContinuityErrorCompensationMode"
                ) :
                pimple_.dict().lookupOrDefault<word>
                (
                    "continuityErrorCompensationMode",
                    "SuSp"
                )
        )
    ),
    heContErrCompensationMode_
    (
        contErrCompensationModeNames_.get
        (
            (pimple_.dict().found("heContinuityErrorCompensationMode")) ?
                pimple_.dict().get<word>
                (
                    "heContinuityErrorCompensationMode"
                ) :
                pimple_.dict().lookupOrDefault<word>
                (
                    "continuityErrorCompensationMode",
                    "SuSp"
                )
        )
    ),
    heStabilizationMode_
    (
        heStabilizationModeNames_.get
        (
            pimple_.dict().lookupOrDefault<word>
            (
                "enthalpyStabilizationMode",
                "cutoff"
            )
        )
    ),
    oscillationLimiterFraction_
    (
        pimple.dict().lookupOrDefault<scalar>
        (
            "oscillationLimiterFraction", 
            0.0
        )
    )
{
    //- Init autoPtr-managed fields, namely flowQuality, XLM (i.e. 
    //  Lockhart-Martinelli parameter) and dispersion
    fluid1_.initTwoPhaseFields();
    fluid2_.initTwoPhaseFields();

    //- Construct twoPhaseDragFactor (needs to be done AFTER fluid twoPhase
    //  field init as some fields of the fluid class that can be required by
    //  some twoPhaseDragFactor sub-models are not initialized yet)
    if (this->subDict("physicsModels").found("twoPhaseDragMultiplierModel"))
    {
        if 
        (
            this->subDict
            (
                "physicsModels"
            ).subDict
            (
                "twoPhaseDragMultiplierModel"
            ).toc().size() != 0
        )
        {
            twoPhaseDragFactorPtr_.reset
            (
                new twoPhaseDragFactor
                (
                    F1SPair_,
                    F2SPair_,
                    *this
                )
            );
        }
    }
    //- This is specifically to avoid problem when solveAlpha is set to solve
    //  for only one phase, the one for which BCs and initial conditions are
    //  provided, yet this phase starts at 0 while the other phase BC and
    //  initial conditions were not provided. This is problematic as by 
    //  default a phase defaults to 0, so if the provided phase is 0 everywhere
    //  too, you have issues with the normalization later. This is veeeery
    //  specific but hey, why not!
    IOobject fluid1Header
    (
        "alpha."+fluid1_.name(),
        mesh_.time().timeName(),
        mesh_,
        IOobject::NO_READ
    );
    IOobject fluid2Header
    (
        "alpha."+fluid2_.name(),
        mesh_.time().timeName(),
        mesh_,
        IOobject::NO_READ
    );
    if 
    (
        !fluid1Header.typeHeaderOk<volScalarField>(true)
    and fluid2Header.typeHeaderOk<volScalarField>(true)
    )
    {
        fluid1_.volScalarField::operator=(geometricOneField()-fluid2_);
    }
    if 
    (
        !fluid2Header.typeHeaderOk<volScalarField>(true)
    and fluid1Header.typeHeaderOk<volScalarField>(true)
    )
    {
        fluid2_.volScalarField::operator=(geometricOneField()-fluid1_);
    }

    //- Normalize phase fraction fields, structure is left unchanged
    volScalarField corr(movingAlpha_/(fluid1_+fluid2_));
    fluid1_.volScalarField::operator*=(corr);
    fluid2_.volScalarField::operator*=(corr);
    fluid1_.correctBoundaryConditions();
    fluid2_.correctBoundaryConditions();

    //- Set the normalized phase fraction fields
    fluid1_.normalized() = fluid1_/movingAlpha_;
    fluid2_.normalized() = fluid2_/movingAlpha_;

    //- Set the bothPhasesArePresentFlag. This is only used within the
    //  adjustTimeStep function and for avoiding doing subcycles if
    //  only one phase is present
    bothPhasesArePresent_ = 
    (
        max((fluid1_*fluid2_)()).value() >= 1e-4 
    );

    //- Initialize fluid-intensive fluxes (i.e. that depend on the phase
    //  fraction, namely alphaPhi and alphaRhoPhi, which are the REAL 
    //  volumetric flux in m3/s and the REAL mass flux in kg/s. By REAL I mean
    //  not superficial). This is done after the phaseFraction normalization 
    //  step to ensure consistency. This step has an effect ONLY IF the
    //  alphaPhi, alphaRhoPhi fields were NOT found on disk
    fluid1_.initAlphaPhis();
    fluid2_.initAlphaPhis();

    //- Set total volumetric flux (real one, not superficial)
    phi_ = fluid1_.alphaPhi() + fluid2_.alphaPhi();

    //- Create turbulence models. This is done outside of fluid constructors as
    //  turbulence models might require references to fields that do not exist
    //  yet (e.g. the Reynolds number between fluid and structure, or between
    //  the two fluids, which are in the FFPair/FSPair classes, that cannot be
    //  created before the fluids)
    fluid1_.constructTurbulenceModel();
    fluid2_.constructTurbulenceModel();

    //- Construct fluid diameter models
    fluid1_.constructDiameterModel();
    fluid2_.constructDiameterModel();

    //- Compute initialFluidMass
    //initialFluidMass_ = fvc::domainIntegrate(fluid_.rho()*fluid_);

    Info << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

//- Solve according to flags
void Foam::thermalHydraulicsModels::twoPhase::correct
(
    scalar& residual,
    bool solveFluidDynamics, 
    bool solveEnergy
)
{   
    correctModels(solveFluidDynamics, solveEnergy);
    
    //auto start = std::chrono::steady_clock::now();
    if (solveFluidDynamics)
    {
        correctFluidMechanics(residual);
    }
    
    /*if (solveEnergy and pimple_.finalIter())
    {
        int N(pimple_.dict().getOrDefault<int>("nEnergyIter", 1));
        for(int c = 1; c < N+1; c++)
        {
            Info << "\nEnergy: iteration " << c << endl;
            //if (c != 1)
                //correctModels(solveFluidDynamics, solveEnergy);
            correctEnergy(residual);
        }
    }*/

    if (solveEnergy)
    {
        correctEnergy(residual);
    }
    //auto end = std::chrono::steady_clock::now();
    //std::chrono::duration<double> elapsed_seconds = end-start;
    //Info << "solveAll: " << elapsed_seconds.count() << "s\n";

    Info << endl;
}

void Foam::thermalHydraulicsModels::twoPhase::correctFluidMechanics
(
    scalar& residual
)
{
    #include "alphaEqns_2p.H"
    #include "UEqns_2p.H"
    if (momentumMode_ == momentumMode::faceCentered)
    {
        #include "pEqnf_2p.H"
    }
    else
    {
        #include "pEqn_2p.H"
    }

    //- Continuity error adjustment and infos
    correctContErrs();
    printContErrs();
    calcCumulContErrs();
}

void Foam::thermalHydraulicsModels::twoPhase::correctEnergy(scalar& residual)
{
    #include "EEqns_2p.H"
}

void Foam::thermalHydraulicsModels::twoPhase::correctModels
(
    bool solveFluidDynamics, 
    bool solveEnergy
)
{
    //- Two-phase Reynolds. Logically it does not belong anywere else as it
    //  involves both fluids and the structure (via the hydraulic diameter)
    /*
    ReTwoPhase_ = 
        max
        (
            mag
            (
                fluid1_.normalized()*fluid1_.rho()*fluid1_.U()
            +   fluid2_.normalized()*fluid2_.rho()*fluid2_.U()
            )*structure_.Dh()/
            (
                fluid1_.normalized()*fluid1_.thermo().mu()
            +   fluid2_.normalized()*fluid2_.thermo().mu()
            ),
            dimensionedScalar("", dimless, 10)
        );*/
    //auto start0 = std::chrono::steady_clock::now();
    this->correctRegimeMaps();
    //auto end0 = std::chrono::steady_clock::now();
    //std::chrono::duration<double> elapsed_seconds0 = end0-start0;
    //Info << "correctRegimes: " << elapsed_seconds0.count() << "s\n";
    //auto start = std::chrono::steady_clock::now();
    FFPair_.correct(solveFluidDynamics, solveEnergy);
    F1SPair_.correct(solveFluidDynamics, solveEnergy);
    F2SPair_.correct(solveFluidDynamics, solveEnergy);
    if (twoPhaseDragFactorPtr_.valid())
        twoPhaseDragFactorPtr_->correct();
    //auto end = std::chrono::steady_clock::now();
    //std::chrono::duration<double> elapsed_seconds = end-start;
    //Info << "correctModels: " << elapsed_seconds.count() << "s\n";
    /*
    forAll(mesh_.cells(), i)
    {
        Info<< i << " " << F1SPair_.Re()[i] << " " << fluid1_.dispersion()[i] 
            << " " << fluid2_.dispersion()[i] << " " << F2SPair_.f()[i] 
            << endl;
    }
    */
}

void Foam::thermalHydraulicsModels::twoPhase::correctCourant()
{
    CoNum_ = 0.0;
    meanCoNum_ = 0.0;

    scalarField sumPhi
    (
        fvc::surfaceSum(mag(phi_))().primitiveField()/
        movingAlpha_.primitiveField()
    );

    CoNum_ = 0.5*gMax(sumPhi/mesh_.V().field())*runTime_.deltaTValue();

    meanCoNum_ =
        0.5*(gSum(sumPhi)/gSum(mesh_.V().field()))*runTime_.deltaTValue();

    Info<< "Courant Number (avg max) = " << meanCoNum_
        << " " << CoNum_ << endl;

    scalar UrCoNum = 
        0.5*gMax
        (
            fvc::surfaceSum
            (
                mag(fluid1_.phi()-fluid2_.phi())
            )().primitiveField()/
            movingAlpha_/
            mesh_.V().field()
        )*runTime_.deltaTValue();

    Info<< "Max Ur Courant Number = " << UrCoNum << endl;

    CoNum_ = max(CoNum_, UrCoNum);
}


void Foam::thermalHydraulicsModels::twoPhase::adjustTimeStep()
{
    this->correctCourant();

    bool phase1(max(fluid1_).value() >= 1e-6);
    bool phase2(max(fluid2_).value() >= 1e-6);
    bothPhasesArePresent_ = (phase1 and phase2);

    bool adjustTimeStep =
        runTime_.controlDict().lookupOrDefault("adjustTimeStep", false);

    if (adjustTimeStep)
    {
        scalar maxCo =
            runTime_.controlDict().lookupOrDefault<scalar>("maxCo", 1.0);

        if (runTime_.controlDict().found("maxCoTwoPhase"))
        {
            scalar maxCoTwoPhase
            (
                runTime_.controlDict().get<scalar>("maxCoTwoPhase")
            );

            if (bothPhasesArePresent_) maxCo = maxCoTwoPhase;
            else
            {
                if 
                (
                    runTime_.controlDict().found("marginToPhaseChange")
                //  and phaseChange_.valid()
                )
                {
                    scalar marginToPhaseChange
                    (
                        runTime_.controlDict().get<scalar>
                        (
                            "marginToPhaseChange"
                        )
                    );

                    scalar DT1
                    (
                        min(mag(fluid1_.T()-FFPair_.iT())().primitiveField())
                    );
                    scalar DT2
                    (
                        min(mag(fluid2_.T()-FFPair_.iT())().primitiveField())
                    );
                    if 
                    (
                        (
                            phase1 and !phase2 and DT1 < marginToPhaseChange
                        ) or
                        (
                            phase2 and !phase1 and DT2 < marginToPhaseChange
                        )
                    )
                    {
                        maxCo = maxCoTwoPhase;
                    }
                }
            }
        }

        scalar maxDeltaT =
            runTime_.controlDict().lookupOrDefault<scalar>("maxDeltaT", GREAT);

        scalar minDeltaT =
            runTime_.controlDict().lookupOrDefault<scalar>("minDeltaT", 1e-69);

        scalar f =
            std::abs
            (
                runTime_.controlDict().lookupOrDefault<scalar>
                (
                    "maxDeltaTMaxRelInc", 
                    0.1
                )
            );
        scalar maxDeltaTFact = maxCo/(CoNum_ + SMALL);
        scalar deltaTFact = 
            min
            (
                min
                (
                    maxDeltaTFact, 
                    1.0 + f*maxDeltaTFact
                ), 
                1.0 + f
            );

        runTime_.setDeltaT
        (
            max
            (
                min
                (
                    deltaTFact*runTime_.deltaTValue(),
                    maxDeltaT
                ),
                minDeltaT
            )
        );
    }
}


void Foam::thermalHydraulicsModels::twoPhase::correctContErrs()
{
    volScalarField& cE1(fluid1_.contErr());
    volScalarField& cE2(fluid2_.contErr());
    volScalarField& rho1(fluid1_.rho());
    volScalarField& rho2(fluid2_.rho());

    cE1 = 
    (
        fvc::ddt(fluid1_, rho1)
    +   fvc::div(fluid1_.alphaRhoPhi())
    -   (fvOptions_(fluid1_, rho1) & rho1)
    );
    cE2 = 
    (
        fvc::ddt(fluid2_, rho2) 
    +   fvc::div(fluid2_.alphaRhoPhi())
    -   (fvOptions_(fluid2_, rho2) & rho2)
    );
    if (FFPair_.phaseChange())
    {
        cE1 += FFPair_.dmdt();
        cE2 -= FFPair_.dmdt();
    }
    
    cE1.correctBoundaryConditions();
    cE2.correctBoundaryConditions();
}


void Foam::thermalHydraulicsModels::twoPhase::printContErrs()
{
    volScalarField contErrRel1(fluid1_.contErr()/fluid1_.rho());
    volScalarField contErrRel2(fluid2_.contErr()/fluid2_.rho());

    Info<< "Instantaneous relative continuity errors (" 
        << fluid1_.name() << " " << fluid2_.name() << ") = "
        << contErrRel1.weightedAverage(mesh_.V()).value() << " "
        << contErrRel2.weightedAverage(mesh_.V()).value() << " "
        << "1/s" << endl;
}

void Foam::thermalHydraulicsModels::twoPhase::calcCumulContErrs()
{
    if (pimple_.finalIter())
    {
        scalar& cumulContErr1(fluid1_.cumulContErr());
        scalar& cumulContErr2(fluid2_.cumulContErr());

        const volScalarField& cE1(fluid1_.contErr());
        const volScalarField& cE2(fluid2_.contErr());

        const scalarField& V(mesh_.V());
        const scalar& dt(mesh_.time().deltaTValue());

        scalar totV(0);
        scalar deltaCumulContErr1(0);
        scalar deltaCumulContErr2(0);

        forAll(V, i)
        {
            const scalar& Vi(V[i]);
            deltaCumulContErr1 += cE1[i]*Vi*dt;
            deltaCumulContErr2 += cE2[i]*Vi*dt;
            totV += Vi;
        }
        reduce(deltaCumulContErr1, sumOp<scalar>());
        reduce(deltaCumulContErr2, sumOp<scalar>());
        reduce(totV, sumOp<scalar>());

        cumulContErr1 += deltaCumulContErr1;
        cumulContErr2 += deltaCumulContErr2;

        Info<< "Cumulative continuity errors ("
            << fluid1_.name() << " " << fluid2_.name() << ") = " 
            << (cumulContErr1/totV) << " " << (cumulContErr2/totV)
            << " kg/m3" << endl;
    }
}


// ************************************************************************* //