neutroResidual = 0.0;

if(aitkenAcceleration)
{
    #include "aitkenAcceleration.H"
}   

forAll(flux_,energyI)
{
    Dalbedo_ = D[energyI];
    fluxStarAlbedo_ *= 0.0;
    #include "calcScatteringSource.H"
    neutroSource_ -= nuSigmaEff[energyI]*flux_[energyI]*(1-eigenvalueNeutronics_);

    fvMatrix<scalar> neutroEq
    (
        neutroEqList[energyI]
      - neutroSource_*(1.0-BetaTot)*chiPrompt[energyI]/keff_
      - delayedNeutroSource_*chiDelayed[energyI]
      - scatteringSourceExtra_ 
    );


    scalar initRes = neutroEq.solve().max().initialResidual();
    neutroResidual = max(neutroResidual,initRes);

    //Info<< "neutro Residual = " << initRes << nl << endl;

    neutroEq.solve();

    flux_[energyI]=fluxStar_[energyI]/discFactor[energyI];
    flux_[energyI].correctBoundaryConditions();
    neutroSource_ += nuSigmaEff[energyI]*flux_[energyI]*(1-eigenvalueNeutronics_);

}




