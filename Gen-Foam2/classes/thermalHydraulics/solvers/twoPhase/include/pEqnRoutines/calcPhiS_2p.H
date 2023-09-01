//- Construct flux due to explicit sources, namely buoyancy/gravity effects
//  and the structure momentum source, if present. Passing "true" to 
//  fluid_.rho() returns the rho calculated according to the Boussinesq 
//  thermophysical model IF the Boussinesq equationOfState is being used.
//  In fact, not passing anything to fluid_.rho() return the constant rho0
//  as a field (which is the rho ref called in UEqns.H). Again, this all
//  happens under the hood in the fluid class, and when NOT using 
//  Boussinesq, fluid_.rho(), fluid_.rho(false), fluid_.rho(true) all
//  return the same rho field, according to the chosen equationOfState
volScalarField& rho1b(fluid1_.rho(true));
volScalarField& rho2b(fluid2_.rho(true));
surfaceScalarField ghSnGradRho
(
    "ghSnGradRho",
    ghf_*
    fvc::snGrad((fluid1_*rho1b+fluid2_*rho2b)/movingAlpha_)*
    mesh_.magSf()
);
surfaceScalarField deltaRhof12
(
    fvc::interpolate(rho1b) - fvc::interpolate(rho2b)
);
surfaceScalarField phiS1
(
    //- Gravity/buoyancy
    -   alphaf1*
        (
            ghSnGradRho
        -   alphaf2/movingAlphaf*
            deltaRhof12*(g_ & mesh_.Sf())
        )
);
surfaceScalarField phiS2
(
    //- Gravity/buoyancy
    -   alphaf2*
        (
            ghSnGradRho
        -   alphaf1/movingAlphaf*
        -   deltaRhof12*(g_ & mesh_.Sf())
        )
);
if (structure_.momentumSourcePtr().valid())
{
    surfaceScalarField fluxMS(fvc::flux(structure_.momentumSource()));
    phiS1 += alphaf1*fluxMS;
    phiS2 += alphaf2*fluxMS;
}

/* NOTE 
    
    To add a new source term S in the momentum equation for phase i, e.g.:
    
    ddt(alpha_i, rho_i, U_i) + ... = + S_i
    
    it is sufficient to add it to phiS_i as:
    
    fvc::flux(S_i)

    if alpha figures in S_i so that it can be written as:

    S_i = alpha_i*Q_i

    it is better to add the following to phiSource

    alphaf_i*fvc:flux(Q_i)

    rather than what was said before
*/
