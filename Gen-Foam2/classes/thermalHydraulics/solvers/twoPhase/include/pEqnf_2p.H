#include "refsAndControls_2p.H"

//- The faceMomentum algorithm relies on manually adding the time derivative
//  to the pressure equation, so that it can be defined on the faces (i.e.
//  deal with it in terms of fluxes). Thus, the time derivative coefficients
//  are defined here
surfaceScalarField alphaRhof1byDt
(
    "alphaRhof1",
    fvc::interpolate(max(alpha1, residualAlpha1)*rho1)/
    mesh_.time().deltaT()
);
surfaceScalarField alphaRhof2byDt
(
    "alphaRhof2",
    fvc::interpolate(max(alpha2, residualAlpha2)*rho2)/
    mesh_.time().deltaT()
);
surfaceScalarField alphaRhof10byDt
(
    "alphaRhof10",
    fvc::interpolate(max(alpha1.oldTime(), residualAlpha1)*rho1.oldTime())/
    mesh_.time().deltaT()
);
surfaceScalarField alphaRhof20byDt
(
    "alphaRhof20",
    fvc::interpolate(max(alpha2.oldTime(), residualAlpha2)*rho2.oldTime())/
    mesh_.time().deltaT()
);

surfaceScalarField Kdf12(fvc::interpolate(Kd12));

if(!implicitPartialElimination)
{
    //- Construct 1/A fields
    surfaceScalarField rAUf1
    (
        "rAUf."+fluid1_.name(),
        1.0/
        (
            fvc::interpolate(UEqn1.A())
        +   alphaRhof1byDt
        )
    );
    surfaceScalarField rAUf2
    (
        "rAUf."+fluid2_.name(),
        1.0/
        (
            fvc::interpolate(UEqn2.A())
        +   alphaRhof2byDt
        )
    );
    surfaceScalarField alpharAUf1(alphaf1*rAUf1);
    surfaceScalarField alpharAUf2(alphaf2*rAUf2);

    bool convergedPISO(false);

    while (pimple_.correct(convergedPISO))
    {
        rho_ = 
            (
                fluid1_*fluid1_.rho() 
            +   fluid2_*fluid2_.rho()
            )/movingAlpha_;
        p_rgh_ = p_ - rho_*gh_;

        //- Correct fixed-flux BCs to be consistent with the velocity BCs
        #include "correctBoundaryFlux_2p.H"

        //- Face-interpolated momentum source terms
        //  Currently includes only gravity/buoyancy and the optional pump 
        //  momentum source. Any new explicit source terms in the momentum  
        //  equation (that were NOT included in the momentum equation) should
        //  be added here 
        #include "calcPhiS_2p.H"

        //- Add explicit time derivative contribution
        phiS1 += alphaRhof10byDt*phi1.oldTime();
        phiS2 += alphaRhof20byDt*phi2.oldTime();

        //- Construct phase fluxes
        surfaceScalarField phiHbyA1
        (
            IOobject::groupName("phiHbyA", fluid1_.name()),
            rAUf1*(fvc::flux(UEqn1.H())+phiS1)
        );
        surfaceScalarField phiHbyA2
        (
            IOobject::groupName("phiHbyA", fluid2_.name()),
            rAUf2*(fvc::flux(UEqn2.H())+phiS2)
        );
        if (noPartialElimination)
        {
            phiHbyA1 += rAUf1*Kdf12*phi2;
            phiHbyA2 += rAUf2*Kdf12*phi1;
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
                alphaf1*rAUf1*Kdf12*phi2
            +   alphaf2*rAUf2*Kdf12*phi1;
        } 

        //- Construct mixture 1/A field
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
                        (
                            (phi2 + rAUf2*Kdf12*phi1) 
                        -   (phi1 + rAUf1*Kdf12*phi2)
                        )/(1.0 - rAUf1*rAUf2*sqr(Kdf12))
                    );
                    
                    phi1 = (phi_ - alphaf2*phir)/movingAlphaf;
                    phi2 = (phi_ + alphaf1*phir)/movingAlphaf;
                }

                p_rgh_.relax();

                fluid1_.dgdt() = -pEqnComp1 & p_rgh_;
                fluid2_.dgdt() = -pEqnComp2 & p_rgh_;

                #include "reconstructU_2p.H"
            }   
        }

        #include "updatePP_rghRho_2p.H"
    }
}
else //- i.e. if doing implicitPartialElimination
{
    //- Diagonal coefficients with time derivative contribution
    surfaceScalarField Af1(fvc::interpolate(UEqn1.A()) + alphaRhof1byDt);
    surfaceScalarField Af2(fvc::interpolate(UEqn2.A()) + alphaRhof2byDt);

    //- Construct 1/A fields with additional terms coming from Partial 
    //  Elimination
    autoPtr<surfaceScalarField> KdbyAf1Ptr;
    autoPtr<surfaceScalarField> KdbyAf2Ptr;
    autoPtr<surfaceScalarField> rAUf1Ptr;
    autoPtr<surfaceScalarField> rAUf2Ptr;
    autoPtr<surfaceScalarField> phiSu1Ptr;
    autoPtr<surfaceScalarField> phiSu2Ptr;
    
    KdbyAf1Ptr.reset(new surfaceScalarField(Kdf12/Af1));
    KdbyAf2Ptr.reset(new surfaceScalarField(Kdf12/Af2));
    
    //- These are the actual 1/A fields
    rAUf1Ptr.reset
    (
        new surfaceScalarField(1.0/(Af1 + Kdf12*(1.0+Af1/Af2)))
    );
    rAUf2Ptr.reset
    (
        new surfaceScalarField(1.0/(Af2 + Kdf12*(1.0+Af2/Af1)))
    );
    surfaceScalarField& KdbyAf1(KdbyAf1Ptr());
    surfaceScalarField& KdbyAf2(KdbyAf2Ptr());
    surfaceScalarField& rAUf1(rAUf1Ptr());
    surfaceScalarField& rAUf2(rAUf2Ptr());

    surfaceScalarField alpharAUf1(alphaf1*rAUf1);
    surfaceScalarField alpharAUf2(alphaf2*rAUf2);
    surfaceScalarField alpharAUf1Star
    (
        (alphaf1+KdbyAf2*(alphaf1+alphaf2))*rAUf1
    );
    surfaceScalarField alpharAUf2Star
    (
        (alphaf2+KdbyAf1*(alphaf1+alphaf2))*rAUf2
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

        //- Correct fixed-flux BCs to be consistent with the velocity BCs
        #include "correctBoundaryFlux_2p.H"

        //- Face-interpolated momentum source terms.
        //  Currently includes only gravity/buoyancy and the optional pump 
        //  momentum source. Any new explicit source terms in the momentum  
        //  equation (that were NOT included in the momentum equation) should
        //  be added here 
        #include "calcPhiS_2p.H"

        //- Add explicit time derivative contribution
        phiS1 += alphaRhof10byDt*phi1.oldTime();
        phiS2 += alphaRhof20byDt*phi2.oldTime();

        //- Construct phase fluxes
        surfaceScalarField phiH1(fvc::flux(UEqn1.H())+phiS1);
        surfaceScalarField phiH2(fvc::flux(UEqn2.H())+phiS2);
        surfaceScalarField phiHbyA1
        (
            IOobject::groupName("phiHbyA", fluid1_.name()),
            rAUf1*
            (
                (1.0+KdbyAf2)*phiH1
            +   KdbyAf2*phiH2
            )
        );
        surfaceScalarField phiHbyA2
        (
            IOobject::groupName("phiHbyA", fluid2_.name()),
            rAUf2*
            (
                (1.0+KdbyAf1)*phiH2
            +   KdbyAf1*phiH1
            )
        );
        if (phiSu1Ptr.valid())
        {
            phiHbyA1 -= rAUf1*phiSu1Ptr();
        }
        if (phiSu2Ptr.valid())
        {
            phiHbyA2 -= rAUf2*phiSu2Ptr();
        }

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
                surfaceScalarField mSfGradp("mSfGradp", pEqnFlux/rAUf);

                phi_ = phiHbyA + pEqnFlux;
                phi1 = phiHbyA1 + alpharAUf1Star*mSfGradp;
                phi2 = phiHbyA2 + alpharAUf2Star*mSfGradp;

                p_rgh_.relax();

                fluid1_.dgdt() = -pEqnComp1 & p_rgh_;
                fluid2_.dgdt() = -pEqnComp2 & p_rgh_;
                
                #include "reconstructU_2p.H"
            }
        }

        #include "updatePP_rghRho_2p.H"
    }
}
