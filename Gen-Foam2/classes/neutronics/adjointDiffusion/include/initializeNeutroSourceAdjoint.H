
adjointFlux_chiP_ *= 0.0;
forAll(flux_,energyI)
{
	adjointFlux_chiP_ += chiPrompt[energyI]*flux_[energyI];
}
