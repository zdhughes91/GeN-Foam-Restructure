
forAll(fluxStar_, energyI)
{
    fluxStar_[energyI] *= pTarget_/pTot;
    fluxStar2_[energyI] *= pTarget_/pTot;
    flux_[energyI] = (fluxStar_[energyI]-fluxStar2_[energyI]*2.0)/discFactor[energyI];
}

