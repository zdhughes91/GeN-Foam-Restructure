//- Solve MULES only for solveFluid, while otherFluid is computed
//  to ensure continuity. This is an alternative to always solving both
//  phase fractions with MULES. Needs testing to see which one performs
//  better in which scenarios
auto solveAlpha = [&]
(
    fluid& solveFluid, 
    fluid& otherFluid
)
{
    volScalarField& solveAlpha(solveFluid);
    volScalarField& otherAlpha(otherFluid);

    //- The source for the phase fraction equation is the dgdt term
    //  computed in pEqn for each fluid. 
    //  The phase fraction equation should be:
    //  
    //  ddt(alpha_i) + div(phi_i, alpha_i) = dgdt
    //
    //  However, I'll cheat and treat dgdt as if it was made of:
    //
    //  dgdt = alpha*Sp + Su
    //
    //  and I will treat the Sp term implicitly ONLY when it contributes
    //  positively to the discretized matrix diagonal, i.e. when 
    //  -alpha*Sp (moved it to the left hand side) is positive, i.e. only
    //  when Sp is negative. Thus Sp = dgdt/alpha when dgdt < 0. As it 
    //  necessitates a division by alpha, it is limited to a minimum alpha,
    //  below which, even if dgdt is nevagative, it is treated explicitly
    //  in Su. This approach is also used in the MULES variant of the
    //  alphaEqns. The minimum alpha is minAlphaSp. If set to 0, dgdt is
    //  always treated entirely explicitly.
    scalar minAlphaSp
    (
        alphaDict.lookupOrDefault<scalar>("minAlphaSp", 1e-4)
    );
    minAlphaSp = max(minAlphaSp, 1e-69);
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
        const scalar& alpha(solveFluid[celli]);
        if (alpha >= minAlphaSp and dgdt[celli] < 0.0)
        {
            Sp[celli] += dgdt[celli]/alpha;
        }
        else
        {
            Su[celli] += dgdt[celli];
        }
    }

    //- Just upwind the shit out of everything and enjoy first-order
    //  boundedness yolo
    upwind<scalar> upwindScheme(mesh_, phi_);
    Foam::fv::gaussConvectionScheme<scalar> gaussUpwindScheme
    (
        mesh_, 
        phi_, 
        upwindScheme
    );

    fvScalarMatrix alphaEqn 
    (
        fvm::ddt(solveAlpha)
    +   gaussUpwindScheme.fvmDiv(solveFluid.phi(), solveAlpha)
    ==  
        fvm::Sp(Sp, solveAlpha)
    +   Su
    );
    alphaEqn.relax();
    myOps::storePrevIterIfRelax(solveAlpha);
    alphaEqn.solve();
    solveAlpha.relax();
    solveAlpha.correctBoundaryConditions();

    solveAlpha = min(max(solveFluid, 0.0), movingAlpha_);
    otherAlpha = 1.0-alphas-solveAlpha;

    surfaceScalarField solveAlphaf(upwindScheme.interpolate(solveAlpha));
    solveFluid.alphaPhi() = solveAlphaf*solveFluid.phi();
    correctInflowOutflow
    (
        solveFluid, 
        solveFluid.phi(), 
        solveFluid.alphaPhi()
    );
    otherFluid.alphaPhi() = phi_-solveFluid.alphaPhi();
    correctInflowOutflow
    (
        otherFluid, 
        otherFluid.phi(), 
        otherFluid.alphaPhi()
    );

    //- For consistency, interpolate the rhos in an upwind manner too
    solveFluid.alphaRhoPhi() = 
        upwindScheme.interpolate(solveFluid.rho())*solveFluid.alphaPhi();
    otherFluid.alphaRhoPhi() = 
        upwindScheme.interpolate(otherFluid.rho())*otherFluid.alphaPhi();
};

word solvePhase
(
    alphaDict.get<word>("solvePhase")
);
if (solvePhase == fluid1_.name())
{
    solveAlpha
    (
        fluid1_,
        fluid2_
    );
}
else if (solvePhase == fluid2_.name())
{
    solveAlpha
    (
        fluid2_,
        fluid1_
    );
}
else
    FatalErrorInFunction 
    << "Phase named " << solvePhase << " not found. Available phases " 
    << "are " << fluid1_.name() << " and " << fluid2_.name()
    << exit(FatalError);
