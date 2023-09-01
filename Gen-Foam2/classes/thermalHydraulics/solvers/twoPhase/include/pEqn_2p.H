#include "refsAndControls_2p.H"

if (!implicitPartialElimination)
{
    //- Construct 1/A fields and limit to avoid problems for alpha = 0
    volScalarField rAU1
    (
        "rAU."+fluid1_.name(),
        1.0/
        (
            UEqn1.A()
        +   max(residualAlpha1 - alpha1, 0.0)*
            rho1/mesh_.time().deltaT()
        )
    );
    volScalarField rAU2
    (
        "rAU."+fluid2_.name(),
        1.0/
        (
            UEqn2.A()
        +   max(residualAlpha2 - alpha2, 0.0)*
            rho2/mesh_.time().deltaT()
        )
    );

    //- Fields required either way, whether doing explicit Partial Elimination
    //  or not
    volScalarField rAUKd1(rAU1*Kd12);
    volScalarField rAUKd2(rAU2*Kd12);
    surfaceScalarField rAUKdf1(fvc::interpolate(rAUKd1));
    surfaceScalarField rAUKdf2(fvc::interpolate(rAUKd2));

    //- Interpolate 1/A to faces and limit
    surfaceScalarField alpharAUf1
    (
        fvc::interpolate(alpha1*rAU1)
    );
    surfaceScalarField alpharAUf2
    (
        fvc::interpolate(alpha2*rAU2)
    );

    bool convergedPISO(false);

    while (pimple_.correct(convergedPISO))
    {
        rho_ = 
            (
                fluid1_*fluid1_.rho() 
            +   fluid2_*fluid2_.rho()
            )/movingAlpha_;
        p_rgh_ = p_ - rho_*gh_;

        //- Correct boundary fluxes
        #include "correctBoundaryFlux_2p.H"

        //- Construct H/A fields
        volVectorField HbyA1
        (
            "HbyA."+fluid1_.name(),
            rAU1*UEqn1.H()
        );
        volVectorField HbyA2
        (
            "HbyA."+fluid2_.name(),
            rAU2*UEqn2.H()
        );

        //- Face-interpolated momentum source terms.
        //  Currently includes only gravity/buoyancy and the optional pump 
        //  momentum source. Any new explicit source terms in the momentum  
        //  equation (that were NOT included in the momentum equation) should
        //  be added here 
        #include "calcPhiS_2p.H"

        surfaceScalarField phiSbyA1(fvc::interpolate(rAU1)*phiS1);
        surfaceScalarField phiSbyA2(fvc::interpolate(rAU2)*phiS2);

        //- Construct phase fluxes
        surfaceScalarField phiHbyA1
        (
            "phiHbyA."+fluid1_.name(),
            fvc::flux(HbyA1)
        +   phiSbyA1
        );
        surfaceScalarField phiHbyA2
        (
            "phiHbyA."+fluid2_.name(),
            fvc::flux(HbyA2)
        +   phiSbyA2
        );
        //- Add back explicit interphase drag contributions if not doing
        //  Partial Elimination
        if (noPartialElimination)
        {
            HbyA1 += rAUKd1*U2;
            HbyA2 += rAUKd2*U1;
            phiHbyA1 += rAUKdf1*phi2;
            phiHbyA2 += rAUKdf2*phi1;
        }
        
        //- Construct mixture flux
        surfaceScalarField phiHbyA
        (
            "phiHbyA",
            alphaf1*phiHbyA1
        +   alphaf2*phiHbyA2
        );
        if (explicitPartialElimination)
        {
            phiHbyA += 
                alphaf1*rAUKdf1*phi2
            +   alphaf2*rAUKdf2*phi1;
        }

        //- Interpolate mixture 1/A field
        surfaceScalarField rAUf
        (
            "rAUf",
            mag(alphaf1*alpharAUf1 + alphaf2*alpharAUf2)
        );

        #include "setSnGradP_rgh_2p.H"

        #include "calcComp_2p.H"

        //- Cache p prior to solve for density update
        volScalarField p_rgh_0(p_rgh_);

        //- Iterate over the pressure equation to correct for non-orthogonality
        while (pimple_.correctNonOrthogonal())
        {
            #include "solveP_rgh_2p.H"

            if (pimple_.finalNonOrthogonalIter())
            {
                surfaceScalarField pEqnFlux(pEqnIncomp.flux());

                phi_ = phiHbyA + pEqnFlux;

                surfaceScalarField mSfGradp("mSfGradp", pEqnFlux/rAUf);
                
                phi1 = phiHbyA1 + alpharAUf1*mSfGradp;
                phi2 = phiHbyA2 + alpharAUf2*mSfGradp;

                if (explicitPartialElimination)
                {
                    surfaceScalarField phir
                    (
                        ((phi1 + rAUKdf1*phi2) - (phi2 + rAUKdf2*phi1))/
                        (1.0 - rAUKdf1*rAUKdf2)
                    );

                    phi1 = (phi_ + alphaf2*phir)/movingAlphaf;
                    phi2 = (phi_ - alphaf1*phir)/movingAlphaf;
                }

                p_rgh_.relax();

                fluid1_.dgdt() = -pEqnComp1 & p_rgh_;
                fluid2_.dgdt() = -pEqnComp2 & p_rgh_;

                mSfGradp = pEqnIncomp.flux()/rAUf;

                if 
                (
                    momentumMode_ 
                ==  momentumMode::cellCenteredFaceReconstruction
                )
                {
                    #include "reconstructU_2p.H"
                }
                else
                {
                    volVectorField Us1
                    ( 
                        HbyA1 
                    +   fvc::reconstruct(alpharAUf1*mSfGradp + phiSbyA1)
                    );
                    volVectorField Us2
                    ( 
                        HbyA2 
                    +   fvc::reconstruct(alpharAUf2*mSfGradp + phiSbyA2)
                    );

                    if (explicitPartialElimination)
                    {
                        volVectorField U
                        (
                            alpha1*(Us1 + rAUKd1*U2) + alpha2*(Us2 + rAUKd2*U1)
                        );

                        volVectorField Ur
                        (
                            ((1.0 - rAUKd2)*Us1 - (1.0 - rAUKd1)*Us2)/
                            (1.0 - rAUKd1*rAUKd2)
                        );

                        U1 = (U + alpha2*Ur)/movingAlpha_;
                        U2 = (U - alpha1*Ur)/movingAlpha_;
                    }
                    else
                    {
                        U1 = Us1;
                        U2 = Us2;
                    }

                    U1.correctBoundaryConditions();
                    fvOptions_.correct(U1);
                    fluid1_.magU() = mag(U1);

                    U2.correctBoundaryConditions();
                    fvOptions_.correct(U2);
                    fluid2_.magU() = mag(U2);

                    //- Update mixture velocity
                    U_ = fluid1_*U1 + fluid2_*U2;
                }
            }
        }
        
        #include "updatePP_rghRho_2p.H"
    }
}
else    //- i.e. if doing implicitPartialElimination
{
    //- Diagonal coefficients with time derivative contribution
    volScalarField A1
    (
        UEqn1.A() 
    +   max(residualAlpha1 - alpha1, 0.0)*rho1/mesh_.time().deltaT()
    );
    volScalarField A2
    (
        UEqn2.A() 
    +   max(residualAlpha2 - alpha2, 0.0)*rho2/mesh_.time().deltaT()
    );

    //- Drag factor by diagonal coeffs required throughout the algorithm by
    //  Partial Elimination
    volScalarField KdbyA1(Kd12/A1);
    volScalarField KdbyA2(Kd12/A2);

    //- Construct 1/A fields with additional terms coming from Partial 
    //  Elimination
    volScalarField rAU1
    (
        "rAU."+fluid1_.name(),
        1.0/(A1 + Kd12*(1.0+A1/A2))
    );
    volScalarField rAU2
    (
        "rAU."+fluid2_.name(),
        1.0/(A2 + Kd12*(1.0+A2/A1))
    );
    surfaceScalarField alpharAUf1(fvc::interpolate(alpha1*rAU1));
    surfaceScalarField alpharAUf2(fvc::interpolate(alpha2*rAU2));
    surfaceScalarField alpharAUf1Star
    (
        fvc::interpolate((alpha1+KdbyA2*(alpha1+alpha2))*rAU1)
    );
    surfaceScalarField alpharAUf2Star
    (
        fvc::interpolate((alpha2+KdbyA1*(alpha1+alpha2))*rAU2)
    );

    bool convergedPISO(false);

    while (pimple_.correct(convergedPISO))
    {
        rho_ = 
            (
                fluid1_*fluid1_.rho() 
            +   fluid2_*fluid2_.rho()
            )/movingAlpha_;
        p_rgh_ = p_ - rho_*gh_;

        //- Correct boundary fluxes
        #include "correctBoundaryFlux_2p.H"

        //- Construct H/A fields and limit to avoid problems for alpha = 0
        volVectorField H1(UEqn1.H());
        volVectorField H2(UEqn2.H());
        volVectorField HbyA1
        (
            "HbyA."+fluid1_.name(),
            rAU1*((1.0+KdbyA2)*H1 + KdbyA2*H2)
        );
        volVectorField HbyA2
        (
            "HbyA."+fluid2_.name(),
            rAU2*((1.0+KdbyA1)*H2 + KdbyA1*H1)
        );
        
        //- Face-interpolated momentum source terms
        //  Currently includes only gravity/buoyancy and the optional pump 
        //  momentum source. Any new explicit source terms in the momentum  
        //  equation (that were NOT included in the momentum equation) should
        //  be added here    
        #include "calcPhiS_2p.H"

        surfaceScalarField phiSbyA1
        (
            fvc::interpolate(rAU1*(1.0+KdbyA2))*phiS1
        +   fvc::interpolate(rAU1*KdbyA2)*phiS2
        ); 
        surfaceScalarField phiSbyA2
        (
            fvc::interpolate(rAU2*(1.0+KdbyA1))*phiS2
        +   fvc::interpolate(rAU2*KdbyA1)*phiS1
        );

        //- Construct phase fluxes
        surfaceScalarField phiHbyA1
        (
            "phiHbyA."+fluid1_.name(),
            fvc::flux(HbyA1)
        +   phiSbyA1
        );
        surfaceScalarField phiHbyA2
        (
            "phiHbyA."+fluid2_.name(),
            fvc::flux(HbyA2)
        +   phiSbyA2
        );
        
        //- Construct mixture flux
        surfaceScalarField phiHbyA
        (
            "phiHbyA",
            alphaf1*phiHbyA1
        +   alphaf2*phiHbyA2
        );

        //- Construct mixture 1/A field
        surfaceScalarField rAUf
        (
            "rAUf",
            mag(alphaf1*alpharAUf1Star + alphaf2*alpharAUf2Star)
        );

        #include "setSnGradP_rgh_2p.H"

        #include "calcComp_2p.H"

        //- Cache p prior to solve for density update
        volScalarField p_rgh_0(p_rgh_);

        //- Iterate over the pressure equation to correct for non-orthogonality
        while (pimple_.correctNonOrthogonal())
        {
            #include "solveP_rgh_2p.H"

            if (pimple_.finalNonOrthogonalIter())
            {
                surfaceScalarField pEqnFlux(pEqnIncomp.flux());

                phi_ = phiHbyA + pEqnFlux;

                surfaceScalarField mSfGradp("mSfGradp", pEqnFlux/rAUf);
                
                phi1 = phiHbyA1 + alpharAUf1Star*mSfGradp;
                phi2 = phiHbyA2 + alpharAUf2Star*mSfGradp;

                p_rgh_.relax();

                fluid1_.dgdt() = -pEqnComp1 & p_rgh_;
                fluid2_.dgdt() = -pEqnComp2 & p_rgh_;

                mSfGradp = pEqnIncomp.flux()/rAUf;

                if 
                (
                    momentumMode_ 
                ==  momentumMode::cellCenteredFaceReconstruction
                )
                {
                    #include "reconstructU_2p.H"
                }
                else
                {
                    U1 = 
                        HbyA1 
                    +   fvc::reconstruct
                        (
                            alpharAUf1Star*mSfGradp 
                        +   phiSbyA1
                        );
                    U2 = 
                        HbyA2 
                    +   fvc::reconstruct
                        (
                            alpharAUf2Star*mSfGradp 
                        +   phiSbyA2
                        );

                    U1.correctBoundaryConditions();
                    fvOptions_.correct(U1);
                    fluid1_.magU() = mag(U1);
                    
                    U2.correctBoundaryConditions();
                    fvOptions_.correct(U2);
                    fluid2_.magU() = mag(U2);

                    //- Update mixture velocity
                    U_ = fluid1_*U1 + fluid2_*U2;
                }
            }
        }
        
        #include "updatePP_rghRho_2p.H"
    }
}
