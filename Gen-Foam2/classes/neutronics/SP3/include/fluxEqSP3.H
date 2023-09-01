neutroResidual = 0.0;

if(aitkenAcceleration)
{
    #include "aitkenAcceleration.H"
}   

forAll(flux_,energyI)
{

    Dalbedo_ = D[energyI];
    fluxStarAlbedo_ = fluxStar2_[energyI];

    #include "calcScatteringSource.H"

    fvMatrix<scalar> neutroEq
    (
      neutroEqList[energyI]
      - neutroSource_*(1.0-BetaTot)*chiPrompt[energyI]/keff_
      - delayedNeutroSource_*chiDelayed[energyI]
      - scatteringSourceExtra_
      - (sigmaDisapp[energyI]*(fluxStar2_[energyI]/discFactor[energyI])*2.0) 
      - fvc::ddt((1-eigenvalueNeutronics_)*IV[energyI]/discFactor[energyI],fluxStar2_[energyI])*2.0 
    );

    scalar initRes = neutroEq.solve().max().initialResidual();
    neutroResidual = max(neutroResidual,initRes);

    Info<< "neutro Residual = " << initRes << nl << endl;

    neutroEq.solve();

    fluxStarAlbedo_ = fluxStar_[energyI];
    /*
    Dalbedo_ = D[energyI]/discFactor[energyI]*9.0/7.0; //like in PARCS
    fvMatrix<scalar> neutroEq2 //like in PARCS
    (
      fvm::ddt((1-eigenvalueNeutronics_)*IV[energyI]/discFactor[energyI]*3.0,fluxStar2_[energyI])
      - fvm::laplacian((D[energyI]/discFactor[energyI]*9.0/7.0),fluxStar2_[energyI])
      + fvm::Sp(((1/D[energyI]/3.0)*5.0/3.0+(sigmaDisapp[energyI]*4.0/3.0)),fluxStar2_[energyI]) 
      - (sigmaDisapp[energyI]*(fluxStar_[energyI]/discFactor[energyI]))*2.0/3.0
      + neutroSource_*(1.0-BetaTot)*chiPrompt[energyI]/keff_*2.0/3.0
      + delayedNeutroSource_*chiDelayed[energyI]*2.0/3.0
      + scatteringSourceExtra_*2.0/3.0
      - fvc::ddt((1-eigenvalueNeutronics_)*IV[energyI]/discFactor[energyI],fluxStar_[energyI])*2.0/3.0
     );
    */
    Dalbedo_ = 1.0/(sigmaDisapp[energyI]+sigmaFromTo[energyI][energyI])*3.0/7.0;
    fvMatrix<scalar>  neutroEq2
    (
      fvm::ddt((1-eigenvalueNeutronics_)*IV[energyI]/discFactor[energyI]*3.0,fluxStar2_[energyI])
      - fvm::laplacian((discFactor[energyI]/(sigmaDisapp[energyI]+sigmaFromTo[energyI][energyI])*3.0/7.0),fluxStar2_[energyI])
      + fvm::Sp(((sigmaDisapp[energyI]+sigmaFromTo[energyI][energyI])*5.0/3.0+(sigmaDisapp[energyI]*4.0/3.0)),fluxStar2_[energyI]) 
      - (sigmaDisapp[energyI]*(fluxStar_[energyI]/discFactor[energyI]))*2.0/3.0
      + neutroSource_*(1.0-BetaTot)*chiPrompt[energyI]/keff_*2.0/3.0
      + delayedNeutroSource_*chiDelayed[energyI]*2.0/3.0
      + scatteringSourceExtra_*2.0/3.0 
      - fvc::ddt((1-eigenvalueNeutronics_)*IV[energyI]/discFactor[energyI],fluxStar_[energyI])*2.0/3.0
     );

    initRes = neutroEq2.solve().max().initialResidual();
    neutroResidual = max(neutroResidual,initRes);
    neutroEq2.solve();

    flux_[energyI]=(fluxStar_[energyI]-fluxStar2_[energyI]*2.0)/discFactor[energyI];
    flux_[energyI].correctBoundaryConditions();




}




