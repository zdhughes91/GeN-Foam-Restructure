
intraSources[directionI] = flux_[energyI] * (nuSigmaEff[energyI] * chiPrompt[energyI] * (1 - BetaTot) / keff_ + sigmaFromToTempEJEI[0][energyI]);

if (xs_.legendreMoments() > 1)
{
    forAll(angularFlux_[energyI], directionJ)
    {
        volScalarField intraSourceScatteringCoeff(nuSigmaEff[energyI]*0.0); 
        for(int momentL = 1; momentL < xs_.legendreMoments(); momentL++)
        {
                scalar scatteringCoeff = legendreMatrices_[momentL][directionJ][directionI] * (2.0 * momentL + 1.0);
                intraSourceScatteringCoeff += scatteringCoeff * sigmaFromToTempEJEI[momentL][energyI] ;
        }
        intraSources[directionI] += angularFlux_[energyI][directionJ] * directionWeights_[directionJ] * intraSourceScatteringCoeff ;
    }
}


