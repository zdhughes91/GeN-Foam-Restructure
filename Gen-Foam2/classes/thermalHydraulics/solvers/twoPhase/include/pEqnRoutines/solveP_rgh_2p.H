fvScalarMatrix pEqnIncomp
(
    fvc::div(phiHbyA)
-   fvm::laplacian(rAUf, p_rgh_)
);

fvScalarMatrix pEqn(pEqnComp1() + pEqnComp2() + pEqnIncomp);

pEqn.setReference(pRefCell_, pRefValue_, forcePRef_);

pEqn.relax();

bool converged =
( 
    pEqn.solve
    (
        mesh_.solver(p_rgh_.select(pimple_.finalInnerIter()))
    ).max().nIterations()
==  0.0 
);

if (pimple_.corrNonOrtho() == 1)
    convergedPISO = converged;

