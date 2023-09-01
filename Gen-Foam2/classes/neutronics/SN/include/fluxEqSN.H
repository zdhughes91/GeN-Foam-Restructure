neutroResidual = 0.0;

forAll(flux_,energyI)
{
    #include "calcSourceSN.H"
    scalar innerResidual = 0.0;
    int innerIter = 0;

    do
    {
        innerResidual = 0.0;
        forAll(angularFlux_[energyI], directionI)
        {
            #include "calcThisIntraSource.H"
            facePhi_ = facePhis_[directionI];
            scalar localResidual =  solve
            (
                //neutroEqList[energyI][directionI] // to be used if created before (createFluxMatricesSN)
                fvm::ddt((1-eigenvalueNeutronics_)*IV[energyI],angularFlux_[energyI][directionI])
                + fvm::div(facePhi_, angularFlux_[energyI][directionI],"div(facePhi_,angularFlux_)")
                //+ fvm::Sp(sigmaDisapp[energyI] + sigmaFromToP[0][energyI][energyI], angularFlux_[energyI][directionI])
                + fvm::Sp(sigmaDisapp[energyI] + sigmaFromToTempEJEI[0][energyI], angularFlux_[energyI][directionI])
                ==
                extraSources[directionI]
                + intraSources[directionI]
                + delayedNeutroSource_*chiDelayed[energyI]
            ).max().initialResidual();
            innerResidual = max(innerResidual, localResidual);
            neutroResidual = max(neutroResidual, localResidual);
            //if(eigenvalueNeutronics_) angularFlux_[energyI][directionI] = max(angularFlux_[energyI][directionI],(angularFlux_[energyI][directionI]*0.0));
            
        }
        #include "updateThisScalarFlux.H"
        innerIter++;
    }
    while((innerResidual > (neutroResidual*innToOutTolerance))&&(innerIter < maxInnerIter));
    intraSources.clear();
    extraSources.clear();

}
