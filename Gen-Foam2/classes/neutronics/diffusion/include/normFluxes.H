forAll(fluxStar_, energyI)
{
    fluxStar_[energyI] *= pTarget_/pTot;
    flux_[energyI]=fluxStar_[energyI]/discFactor[energyI];
}
