//- Construct flux due to explicit sources, namely buoyancy/gravity effects
//  and the structure momentum source, if present. Passing "true" to 
//  fluid_.rho() returns the rho calculated according to the Boussinesq 
//  thermophysical model IF the Boussinesq equationOfState is being used.
//  In fact, not passing anything to fluid_.rho() return the constant rho0
//  as a field (which is the rho ref called in UEqns.H). Again, this all
//  happens under the hood in the fluid class, and when NOT using 
//  Boussinesq, fluid_.rho(), fluid_.rho(false), fluid_.rho(true) all
//  return the same rho field, according to the chosen equationOfState
surfaceScalarField phiSbyA
(
    -alpharAUf*ghf_*fvc::snGrad(fluid_.rho(true))*mesh_.magSf()
);
if (structure_.momentumSourcePtr().valid())
{
   phiSbyA += alpharAUf*fvc::flux(structure_.momentumSource());
}

/* NOTE 
    
    To add a new source term S in the momentum equation, e.g.:
    
    ddt(alpha,rho,U) + ... = + S
    
    it is sufficient to add it to phiSbyA as:
    
    rAUf*fvc::flux(S)

    if alpha figures in S so that it can be written as:

    S = alpha*Q

    it is better to add the following to phiSource

    alpharAUf*fvc:flux(Q)

    rather than what was said before
*/
