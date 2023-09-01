surfaceScalarField& phi1(fluid1_.phi());
surfaceScalarField& phi2(fluid2_.phi());
phi1.storePrevIter();
phi2.storePrevIter();
dimensionedScalar residualAlpha1(fluid1_.residualAlpha());
dimensionedScalar residualAlpha2(fluid2_.residualAlpha());

autoPtr<surfaceScalarField> alphaf1Ptr;
autoPtr<surfaceScalarField> alphaf2Ptr;
if (alphaEqnsSolver_ == alphaEqnsSolver::implicitUpwind)
{
    upwind<scalar> upwindScheme(mesh_, phi_);
    alphaf1Ptr.reset(new surfaceScalarField(upwindScheme.interpolate(alpha1)));
    alphaf2Ptr.reset(new surfaceScalarField(upwindScheme.interpolate(alpha2)));
}
else
{
    alphaf1Ptr.reset(new surfaceScalarField(fvc::interpolate(alpha1)));
    alphaf2Ptr.reset(new surfaceScalarField(fvc::interpolate(alpha2)));
}
surfaceScalarField& alphaf1(alphaf1Ptr());
surfaceScalarField& alphaf2(alphaf2Ptr());
surfaceScalarField movingAlphaf(alphaf1+alphaf2);

bool noPartialElimination
(
    partialEliminationMode_ == partialEliminationMode::none
);
bool explicitPartialElimination
(
    partialEliminationMode_ == partialEliminationMode::legacy
);
bool implicitPartialElimination
(
    partialEliminationMode_ == partialEliminationMode::implicit
);
