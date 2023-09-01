
forAll(flux_, energyI)
{
	flux_[energyI] *= pTarget_/pTot;
	forAll(angularFlux_[energyI], directionI)
	{
		angularFlux_[energyI][directionI] *= pTarget_/pTot;
	}
}
