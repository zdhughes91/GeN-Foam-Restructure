{
    const dictionary& alphaDict(mesh_.solverDict("alpha"));
    volScalarField& alpha1(fluid1_);
    volScalarField& alpha2(fluid2_);
    volScalarField& alphas(structure_);

    //- This is a lambda function (standard C++ stuff) that adjusts the
    //  boundary field of a given volumetric flux alphaPhi so that it is
    //  consistent with the BCs for the fields (alpha and U, thus phi) from
    //  which it is derived. It is basically a physically-motivated
    //  correctBoundaryConditions() for alphaPhi
    auto correctInflowOutflow = []
    (
        const volScalarField& alpha,
        const surfaceScalarField& phi,
        surfaceScalarField& alphaPhi
    )
    {
        const volScalarField::Boundary& alphaBf = alpha.boundaryField();
        const surfaceScalarField::Boundary& phiBf = phi.boundaryField();
        surfaceScalarField::Boundary& alphaPhiBf = alphaPhi.boundaryFieldRef();

        forAll(alphaPhiBf, patchi)
        {
            fvsPatchScalarField& alphaPhip = alphaPhiBf[patchi];

            if (!alphaPhip.coupled())
            {
                alphaPhip = phiBf[patchi]*alphaBf[patchi];
            }
        }
    };

    if (alphaEqnsSolver_ == alphaEqnsSolver::implicitUpwind)
    {
        #include "solveAlphaImplicitUpwind.H"
    }
    else 
    {
        #include "solveAlphaMULES.H"
    }
    //- Correct continuity errors
    correctContErrs();

    //- Update normalized alpha fields (only useful for regime map 
    //  parametrization, if selected as regime map parameter)
    fluid1_.normalized() = fluid1_/movingAlpha_;
    fluid2_.normalized() = fluid2_/movingAlpha_;

    //- Infos
    myOps::fieldInfo(alpha1);
    myOps::fieldInfo(alpha2);
    volScalarField alphaSum(alpha1 + alpha2 + alphas);
    alphaSum.rename("sum(alpha)");
    myOps::fieldInfo(alphaSum);
}
