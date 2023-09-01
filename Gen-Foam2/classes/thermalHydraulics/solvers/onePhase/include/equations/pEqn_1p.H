#include "updateAlphaPhiAndInitControls_1p.H"

//- Construct 1/A field and interpolate to faces
volScalarField rAU
(
    "rAU",
    1.0/(UEqn.A())
);
surfaceScalarField alpharAUf
(
    fvc::interpolate(alpha*rAU)
);

bool convergedPISO(false);

//- Pressure corrector loop
while (pimple_.correct(convergedPISO))
{
    p_rgh_ = p_ - rho*gh_;

    //- Correct fixed-flux phi BCs to be consistent with the velocity BCs
    #include "correctBoundaryFlux_1p.H"

    //- Construct H/A field
    volVectorField HbyA
    (
        "HbyA",
        rAU*UEqn.H()
    );

    //- Face-interpolated momentum source terms divided by diagonal coeffs.
    //  Currently includes only gravity/buoyancy and the optional pump momentum
    //  source. Any new explicit source terms in the momentum equation (that
    //  were NOT included in the momentum equation) should be added here
    #include "calcPhiSbyA_1p.H"

    //- Construct phase flux
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        alphaf*
        (
            fvc::flux(HbyA)
        +   phiSbyA
        )
    );

    #include "calcrAUfSetSnGradCalcComp_1p.H"

    //- Cache p prior to solve for density update
    volScalarField p_rgh_0(p_rgh_);

    //- Iterate over the pressure equation to correct for non-orthogonality
    while (pimple_.correctNonOrthogonal())
    {
        #include "solveP_rgh_1p.H"

        if (pimple_.finalNonOrthogonalIter())
        {
            phi_ = phiHbyA + pEqnIncomp.flux();

            fluid_.phi() = phi_/alphaf;

            fluid_.dgdt() = -pEqnComp & p_rgh_;

            p_rgh_.relax();
            
            if (momentumMode_ == momentumMode::cellCenteredFaceReconstruction)
            {
                #include "reconstructU_1p.H"
            }
            else
            {
                U =
                    HbyA 
                +   fvc::reconstruct
                    (
                        alpharAUf*pEqnIncomp.flux()/rAUf + phiSbyA
                    );
            }

            U.correctBoundaryConditions();
            fvOptions_.correct(U);
            fluid_.magU() = mag(U);
        }
    }
    
    //- Includes closed-volume compressibility effects
    #include "updatePP_rghRho_1p.H"
}
