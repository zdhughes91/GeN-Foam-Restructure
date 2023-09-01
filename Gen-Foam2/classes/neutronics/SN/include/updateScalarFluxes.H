forAll(flux_,energyI)
{
    flux_[energyI] *= 0.0;
    forAll(angularFlux_[energyI],directionI)
    { 
        flux_[energyI] += angularFlux_[energyI][directionI] * directionWeights_[directionI];
    }
}

oneGroupFlux_ *= 0.0;
forAll(flux_,energyI)
{
    oneGroupFlux_ += flux_[energyI];
}
