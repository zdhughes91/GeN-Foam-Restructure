neutroResidual = 0.0;

forAll(flux_,energyI)
{
    Dalbedo_ = D[energyI];
    fluxStarAlbedo_ *= 0.0;
    #include "calcScatteringSourceAdjoint.H"
    adjointFlux_chiP_ -= chiPrompt[energyI]*flux_[energyI]*(1-eigenvalueNeutronics_);

    fvMatrix<scalar> neutroEq
    (
      neutroEqList[energyI]
      - adjointFlux_chiP_*(1.0-BetaTot)*nuSigmaEff[energyI]/keff_
      - 1/keff_*nuSigmaEff[energyI]*adjointPrec_beta_
      - scatteringSourceExtra_ 
    );


    scalar initRes = neutroEq.solve().max().initialResidual();
    neutroResidual = max(neutroResidual,initRes);

    Info<< "neutro Residual = " << initRes << nl << endl;

    neutroEq.solve();
   
    flux_[energyI]=fluxStar_[energyI]/discFactor[energyI];
    flux_[energyI].correctBoundaryConditions();
    adjointFlux_chiP_ += chiPrompt[energyI]*flux_[energyI]*(1-eigenvalueNeutronics_);
}


