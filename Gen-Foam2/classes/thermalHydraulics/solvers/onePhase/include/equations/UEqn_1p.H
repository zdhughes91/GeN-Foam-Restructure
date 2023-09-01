//- The velocity predictor is missing in this implementation for
//  consistency with the two-phase treatment. It could be implemented
//  to improve convergence of specifc cases

//- Field refs
volVectorField& U(fluid_.U());
const volScalarField& alpha(fluid_);
volScalarField& rho(fluid_.rho());
fvVectorMatrix UEqn(U, dimMass*dimLength/dimTime/dimTime);
const volTensorField& Kds(FSPair_.Kd());

fluid_.turbulence().correct();

{
    volScalarField nuEff(fluid_.turbulence().nuEff()());

    UEqn =
    (
        fvm::div(fluid_.alphaRhoPhi(), U, "div(alphaRhoPhi,U)")
    +   fvm::SuSp(-fluid_.contErr(), U)
        //- Diffusion term (i.e. fluid_.turbulence().divDevRhoReff(U))
        //  but re-expressed so to: 1) either allow to use a tensor
        //  diffusivity (i.e. nuEff*structure_.tortuosity()) if the 
        //  structure is tortuous; 2) if not tortuous, still allow me to
        //  specify the div interpolation scheme (i.e. "div(alphaRhoPhiNu,U))
    +   (
            (structure_.tortuosityPtr().valid()) ? 
            (
            -   fvm::laplacian
                (
                    alpha*rho*nuEff*structure_.tortuosity(), 
                    U
                )
            -   fvc::div
                (
                    alpha*rho*nuEff*structure_.tortuosity()
                &   dev2(T(fvc::grad(U))), 
                    "div(alphaRhoPhiNu,U)"
                )
            ) :
            (
            -   fvm::laplacian(alpha*rho*nuEff, U)
            -   fvc::div
                (
                    alpha*rho*nuEff*dev2(T(fvc::grad(U))), 
                    "div(alphaRhoPhiNu,U)"
                )
            )
        )
    +   fvm::Sp((1.0/3.0)*tr(Kds), U) + (dev(Kds)&U)
    ==
        fvOptions_(alpha, rho, U)
    );

    if (momentumMode_ != momentumMode::faceCentered)
    {
        UEqn += fvm::ddt(alpha, rho, U);
    }

    UEqn.relax();
    fvOptions_.constrain(UEqn);
    fvOptions_.correct(U);
}
