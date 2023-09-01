
PtrList<fvMatrix<scalar> > neutroEqList(xs_.energyGroups());

forAll(neutroEqList, energyI)
{

    Dalbedo_ = D[energyI];
    fluxStarAlbedo_ = fluxStar2_[energyI];
    neutroEqList.set
    (
        energyI,
        new fvMatrix<scalar>
        (
             fvm::ddt((1-eigenvalueNeutronics_)*IV[energyI]/discFactor[energyI],fluxStar_[energyI])
             - fvm::laplacian(D[energyI]/discFactor[energyI],fluxStar_[energyI])
             + fvm::Sp(sigmaDisapp[energyI]/discFactor[energyI],fluxStar_[energyI])
        )
    );
}


