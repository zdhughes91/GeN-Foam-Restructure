neutroSource_ *= 0.0;
forAll(flux_, energyI)
{
	neutroSource_ += nuSigmaEff[energyI]*flux_[energyI];
}
