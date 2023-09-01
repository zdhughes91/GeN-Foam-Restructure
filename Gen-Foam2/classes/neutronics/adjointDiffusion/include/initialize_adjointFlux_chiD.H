

adjointFlux_chiD_ *= 0.0;
forAll(flux_,energyI)
{
	adjointFlux_chiD_ += chiDelayed[energyI]*flux_[energyI];
}

