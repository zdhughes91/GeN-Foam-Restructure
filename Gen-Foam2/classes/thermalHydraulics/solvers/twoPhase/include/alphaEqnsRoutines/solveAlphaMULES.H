//- Field references
const surfaceScalarField& phi1(fluid1_.phi());
const surfaceScalarField& phi2(fluid2_.phi());
word phirScheme("div(phir,alpha)");
scalar minAlphaSp
(
    alphaDict.lookupOrDefault<scalar>("minAlphaSp", 1e-6)
);

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////    

//- START OF LAMBDA FUNCTION DEFINITIONS

//- Solve MULES only for solveFluid, while otherFluid is computed
//  to ensure continuity. This is an alternative to always solving both
//  phase fractions with MULES. Needs testing to see which one performs
//  better in which scenarios
auto solveAlpha = [&]
(
    fluid& solveFluid, 
    fluid& otherFluid, 
    scalar rDeltaT
)
{
    surfaceScalarField alphaPhic 
    (
        fvc::flux(phi_, solveFluid, "div(phi,alpha)")
    +   fvc::flux
        (
            -fvc::flux
            (
                -(solveFluid.phi()-otherFluid.phi()), 
                movingAlpha_-solveFluid, 
                phirScheme
            ),
            solveFluid,
            phirScheme
        )
    +   fvc::flux
        (
            -fvc::flux(-solveFluid.phi(), structure_, phirScheme),
            solveFluid,
            phirScheme
        )
    );
    correctInflowOutflow(solveFluid, solveFluid.phi(), alphaPhic);
    volScalarField::Internal Sp
    (
        IOobject
        (
            "Sp",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("Sp", dimless/dimTime, 0)
    );
    volScalarField::Internal Su
    (
        
        IOobject
        (
            "Su",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("Su", dimless/dimTime, 0)
    );
    const volScalarField& dgdt(solveFluid.dgdt());
    forAll(dgdt, celli)
    {
        if (minAlphaSp != 0)
        {
            if (dgdt[celli] < 0.0)
            {
                Sp[celli] += 
                    dgdt[celli]/max(solveFluid[celli], minAlphaSp);
            }
            else
            {
                Su[celli] += dgdt[celli];
            }
        }
        else
        {
            Su[celli] += dgdt[celli];
        }
    }
    myMULES::explicitSolve
    (
        rDeltaT,
        geometricOneField(),
        solveFluid,
        phi_,
        alphaPhic,
        Sp,
        Su,
        movingAlpha_.primitiveField(),
        zeroField()
    );

    solveFluid.alphaPhi() = alphaPhic;
    otherFluid.alphaPhi() = phi_-solveFluid.alphaPhi();
    correctInflowOutflow
    (
        otherFluid, 
        otherFluid.phi(), 
        otherFluid.alphaPhi()
    );
    solveFluid.volScalarField::operator=
    (
        min(max(solveFluid, 0.0), movingAlpha_)
    );
    solveFluid.correctBoundaryConditions();
    otherFluid.volScalarField::operator=
    (
        min(max(movingAlpha_-solveFluid, 0.0), movingAlpha_)
    );
};

///////////////////////////////////////////////////////////////////////////////

//- This is a lambda function (standard C++ stuff) that solves the volume
//  fraction equations with the MULES algorithm (explicit) with a given
//  timestep (rDeltaT is its inverse). This function by itself does nothing
//  until it's called, later in this file. Also, the [&] specification
//  allows this function to have access to all of the variables defined in
//  the scope in which the lambda functions is being defined
auto solveAlphas = [&](scalar rDeltaT)
{
    //- Initialize the advective fluxes to be used in the explicit 
    //  continuity equation solution. This is mathematically/physically 
    //  equivalent to writing 
    //  fvc::flux(fluid1_.phi(), alpha1, "div(phi,alpha)"), but written 
    //  this way it enhances numerical stability
    surfaceScalarField alphaPhic1 
    (
        fvc::flux(phi_, alpha1, "div(phi,alpha)")
    +   fvc::flux
        (
            -fvc::flux(-(phi1-phi2), alpha2, phirScheme),
            alpha1,
            phirScheme
        )
    +   fvc::flux
        (
            -fvc::flux(-phi1, alphas, phirScheme),
            alpha1,
            phirScheme
        )
    );
    correctInflowOutflow(alpha1, phi1, alphaPhic1);
    //- Use MULES algorithm to adjust alphaPhic so that it is not overly
    //  diffusive (i.e. purely upwind) yet still bounded. How does it work?
    //  The general idea comes from Flux Correct Transport. The actual
    //  implementation is specific to Weller's work
    myMULES::limit
    (
        rDeltaT,
        geometricOneField(),
        alpha1,
        phi_,
        alphaPhic1,
        zeroField(),
        zeroField(),
        min(movingAlpha_.primitiveField(), 1.0)(),
        zeroField(),
        true
    );
    volScalarField::Internal Sp1
    (
        IOobject
        (
            "Sp1",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("Sp1", dimless/dimTime, 0)
    );
    volScalarField::Internal Su1
    (
        
        IOobject
        (
            "Su1",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("Su1", dimless/dimTime, 0)
    );
    //- See comments in fluid.H concerning the dgdt_ field to understand
    //  its physical meaning
    const volScalarField& dgdt1(fluid1_.dgdt());
    forAll(dgdt1, celli)
    {
        if (minAlphaSp != 0)
        {
            if (dgdt1[celli] < 0.0)
            {
                Sp1[celli] += 
                    dgdt1[celli]/max(alpha1[celli], minAlphaSp);
            }
            else
            {
                Su1[celli] += dgdt1[celli];
            }
        }
        else
        {
            Su1[celli] += dgdt1[celli];
        }
    }
    
    //- Do the same exact things that have been done above, yet do it for
    //  fluid 2
    surfaceScalarField alphaPhic2
    (
        fvc::flux(phi_, alpha2, "div(phi,alpha)")
    +   fvc::flux
        (
            -fvc::flux(-(phi2-phi1), alpha1, phirScheme),
            alpha2,
            phirScheme
        )
    +   fvc::flux
        (
            -fvc::flux(-phi2, alphas, phirScheme),
            alpha2,
            phirScheme
        )
    );
    correctInflowOutflow(alpha2, phi2, alphaPhic2);
    myMULES::limit
    (
        rDeltaT,
        geometricOneField(),
        alpha2,
        phi_,
        alphaPhic2,
        zeroField(),
        zeroField(),
        min(movingAlpha_.primitiveField(), 1.0)(),
        zeroField(),
        true
    );
    volScalarField::Internal Sp2
    (
        IOobject
        (
            "Sp2",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("Sp2", dimless/dimTime, 0)
    );
    volScalarField::Internal Su2
    (
        IOobject
        (
            "Su2",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("Su2", dimless/dimTime, 0)
    );
    const volScalarField& dgdt2(fluid2_.dgdt());
    forAll(dgdt2, celli)
    {
        if (minAlphaSp != 0)
        {
            if (dgdt2[celli] < 0.0)
            {
                Sp2[celli] += 
                    dgdt2[celli]/max(alpha2[celli], minAlphaSp);
            }
            else
            {
                Su2[celli] += dgdt2[celli];
            }
        }
        else
        {
            Su2[celli] += dgdt2[celli];
        }
    }
    
    //- Limit to preserve phi magnitude so that the corrected
    //  alphaPhis still sum up to phi_
    {
        surfaceScalarField alpha1fu
        (
            upwind<scalar>(mesh_, phi_).interpolate(alpha1)
        );
        surfaceScalarField alpha2fu
        (
            upwind<scalar>(mesh_, phi_).interpolate(alpha2)
        );
        surfaceScalarField deltaPhi
        (
            (
                    alphaPhic1
                +   alphaPhic2
                -   phi_
            )/(alpha1fu+alpha2fu)
        );
        alphaPhic1 -= deltaPhi*alpha1fu;
        alphaPhic2 -= deltaPhi*alpha2fu;
    }

    //- Solve fluid1 
    correctInflowOutflow(alpha1, phi1, alphaPhic1);
    myMULES::explicitSolve
    (
        rDeltaT,
        geometricOneField(),
        alpha1,
        alphaPhic1,
        Sp1,
        Su1
    );
    alpha1 = min(max(alpha1, 0.0), 1.0);
    fluid1_.alphaPhi() = alphaPhic1;

    //- Solve fluid2
    correctInflowOutflow(alpha2, phi2, alphaPhic2);
    myMULES::explicitSolve
    (
        rDeltaT,
        geometricOneField(),
        alpha2,
        alphaPhic2,
        Sp2,
        Su2
    );
    alpha2 = min(max(alpha2, 0.0), 1.0);
    fluid2_.alphaPhi() = alphaPhic2;
    
    //- Limit alphas so that they still add up to movingAlpha_ and 
    //  stabilize alphas if their sum goes to 0 (helps with single-phase 
    //  scenarios with Courant numbers > 1 when not doing subCycles).
    //  Please note that, depending on the choice of boundaryConditions for
    //  the phaseFractions, the correctBoundaryConditions might ever so
    //  slightly break the normalization by small amounts
    volScalarField zeroSum(neg0(alpha1+alpha2));
    alpha1 += zeroSum*alpha1.oldTime();
    alpha2 += zeroSum*alpha2.oldTime();
    volScalarField corr(movingAlpha_/(alpha1+alpha2));
    alpha1 *= corr;
    alpha2 *= corr;
    alpha1.correctBoundaryConditions();
    alpha2.correctBoundaryConditions();
}; 

///////////////////////////////////////////////////////////////////////////////

//- Actual function that either solves for one or both phaseFractions
//  via MULES.
auto solve = [&](word solvePhase, scalar rDeltaT)
{
    if (solvePhase == "both")
    {
        solveAlphas(rDeltaT);
    }
    else if (solvePhase == fluid1_.name())
    {
        solveAlpha
        (
            fluid1_,
            fluid2_,
            rDeltaT
        );
    }
    else if (solvePhase == fluid2_.name())
    {
        solveAlpha
        (
            fluid2_,
            fluid1_,
            rDeltaT
        );
    }
    else
        FatalErrorInFunction 
        << "Phase named " << solvePhase << " not found. Available phases " 
        << "are " << fluid1_.name() << " and " << fluid2_.name()
        << exit(FatalError);
};

//- END OF LAMBDA FUNCTION DEFINITIONS

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//- All right, let's talk sub-cycles. The phase fraction equations are the
//  only fully-explicit set of equations, which means that ultimately all
//  my Courant number limitations will come from here, and here only. How
//  can we deal with having larger Courant numbers? The idea is stupidly
//  simple: solve the phase fraction equations over a time step N times
//  smaller for N times, such that the Courant number that is "felt/seen"
//  by the phase fraction equations is smaller by N times. This spares us
//  from having to globally reduce the time step for all of the physics of
//  the solver. This technique is typically referred to as "sub-cycling".
//  While it exists in other multiphase solvers, their implementation is
//  different for resons I can't seem to comprehend, and their
//  implementations do not appear to work with my way of using MULES. 
//  Nonetheless, what I did here works and is as simple as it gets. 
//  Furthermore, I introduced a simple way to adjust the number of
//  sub-cycles at runTime based on a user provided maximum Co num for
//  the alphaEqns (user provided but it should never exceed 0.5)
word solvePhase
(
    alphaDict.lookupOrDefault<word>("solvePhase", "both")
);

bool adjustSubCycles
(
    alphaDict.lookupOrDefault<bool>("adjustSubCycles", false)
);

label maxSubCycles
(
    alphaDict.lookupOrDefault<label>("maxSubCycles", 50)
);

scalar alphaMaxCoNum
(
    (adjustSubCycles) ? alphaDict.get<scalar>("alphaMaxCo") : 0.49
);

label nSubCycles
(
    (adjustSubCycles and bothPhasesArePresent_) ?
    min(ceil((CoNum_+SMALL)/alphaMaxCoNum), maxSubCycles) :
    alphaDict.lookupOrDefault<label>("nSubCycles", 1)
);

Info<< "MULES:  Solving for " 
    << (
            (solvePhase == "both") ?
            word(alpha1.name()+", "+alpha2.name()) :
            word("alpha."+solvePhase)
        ) 
    << ", No SubCycles " << nSubCycles << " (alphaCo = " 
    << (CoNum_/nSubCycles) << ")" << endl;

//- Rescale the inverse of the time-step by nSubCycles (it's the essence of
//  sub-cycling)
scalar rDeltaT(nSubCycles/mesh_.time().deltaTValue());

//- If not doing sub-cycling, solve the phase fraction equations once with 
//  the global time step and call it a day
if (nSubCycles == 1)
{
    solve(solvePhase, rDeltaT);
}
else
{   
    //- Cache old time fields. This is required as the explicitSolve and
    //  limit functions in MULES require access to the old time step
    //  alpha field. But beware, MULES should "see" the old time step
    //  field with respect to the local, shorter deltaT. This is done
    //  later, in the loop. However, when I get out of MULES, I should
    //  reset the old time step fields with respect to the global time
    //  step, not the sub-cycled one. Thus, cache these fields
    volScalarField alpha10(alpha1.oldTime());
    volScalarField alpha20(alpha2.oldTime());
    
    //- Init sub-cycled volumetric flux placeholder fields
    surfaceScalarField alphaPhiSum1
    (
        IOobject
        (
            "phiSum" + alpha1.name(),
            runTime_.timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("", dimVolume/dimTime, 0)
    );
    surfaceScalarField alphaPhiSum2
    (
        IOobject
        (
            "phiSum" + alpha2.name(),
            runTime_.timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("", dimVolume/dimTime, 0)
    );
    
    //- Using the alphaSubCycle class (like ALL the other multiphase 
    //  solvers do) makes everything go wrong. Why? I have no fucking clue.
    //  Nonetheless what I did here works, so...
    //  fuck it all, fuck the world, fuck everything that you stand for
    //  (cit. Slipknot)
    for(label subCycle = 0; subCycle < nSubCycles; subCycle++)
    {
        //- Compute new alphas, alphaPhis with reduced timeStep
        solve(solvePhase, rDeltaT);

        //- Set the computed alphas as oldTime fields (required to make it
        //  work)
        alpha1.oldTime() = alpha1;
        alpha2.oldTime() = alpha2;

        //- Integrate volumetric fluxes over sub-cycles (final ones will be
        //  an average)
        alphaPhiSum1 += fluid1_.alphaPhi();
        alphaPhiSum2 += fluid2_.alphaPhi();
    }

    //- Set volumetric fluxes as an average of the subCycled fluxes
    fluid1_.alphaPhi() = alphaPhiSum1/nSubCycles;
    fluid2_.alphaPhi() = alphaPhiSum2/nSubCycles;

    //- Reset old time fields and time indices to what they were before
    //  sub-cycling
    label timeIndex(runTime_.timeIndex());
    alpha1.oldTime() = alpha10;
    alpha1.timeIndex() = timeIndex;
    alpha1.oldTime().timeIndex() = timeIndex;
    alpha2.oldTime() = alpha20;
    alpha2.timeIndex() = timeIndex;
    alpha2.oldTime().timeIndex() = timeIndex;
}

//- Correct mass fluxes
fluid1_.alphaRhoPhi() = 
        fvc::interpolate(fluid1_.rho())*fluid1_.alphaPhi();
fluid2_.alphaRhoPhi() = 
    fvc::interpolate(fluid2_.rho())*fluid2_.alphaPhi();
    