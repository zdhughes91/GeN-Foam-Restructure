//- Pressure laplacian coefficient
surfaceScalarField rAUf
(
    "rAUf",
    mag(alphaf*alpharAUf)
);

//- Update the fixedFluxPressure BCs to ensure flux consistency. Please
//  make sure that wherever you use this BC, the U is NOT zeroGradient,
//  as the boundary U (and thus phi) are not yet updated at this point.
//  You get un-physical results in those circumstances. Instead, use
//  zeroGradient for pressure were you use zeroGradient for U
setSnGrad<fixedFluxPressureFvPatchScalarField>
(
    p_rgh_.boundaryFieldRef(),
    (
        phiHbyA.boundaryField()
      - (
            alphaf.boundaryField()*fluid_.phi().boundaryField()
        )
    )/(mesh_.magSf().boundaryField()*rAUf.boundaryField())
);

//- Construct the compressibility parts of the pressure equation
tmp<fvScalarMatrix> pEqnComp;
pEqnComp =
        (
            fluid_.contErr() // +dmdt in contErr
          - fvc::Sp(fvc::ddt(alpha) + fvc::div(fluid_.alphaPhi()), rho)
        )/rho
    +   (alpha*fluid_.thermo().psi()/rho)*correction(fvm::ddt(p_rgh_));
