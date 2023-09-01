scatteringSourceExtra_ *= 0.0;

forAll(flux_, energyJ)
{
    if(energyJ != energyI)
    {
        scatteringSourceExtra_ += 
        	sigmaFromTo[energyJ][energyI] * flux_[energyJ];
    }
}
