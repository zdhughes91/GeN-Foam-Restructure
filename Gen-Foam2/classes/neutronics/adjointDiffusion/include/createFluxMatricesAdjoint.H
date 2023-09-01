
PtrList<fvMatrix<scalar> > neutroEqList(xs_.energyGroups());

forAll(neutroEqList, energyI)
{

    Dalbedo_ = D[energyI];
    fluxStarAlbedo_ *= 0.0;
    neutroEqList.set
    (
        energyI,
        new fvMatrix<scalar>
        (
             fvm::ddt((1-eigenvalueNeutronics_)*IV[energyI]/discFactor[energyI],fluxStar_[energyI])
             - fvm::laplacian(D[energyI]/discFactor[energyI],fluxStar_[energyI])
             - fvm::Sp((nuSigmaEff[energyI]/keff_*(1.0-BetaTot)*chiPrompt[energyI]*(1-eigenvalueNeutronics_)-sigmaDisapp[energyI]/discFactor[energyI]),fluxStar_[energyI])
        )
    );
}


