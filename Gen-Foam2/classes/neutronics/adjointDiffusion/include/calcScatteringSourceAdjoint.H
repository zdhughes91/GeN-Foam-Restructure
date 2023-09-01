
scatteringSourceExtra_ *= 0.0;

forAll(flux_,energyJ)
{
    if(energyJ != energyI)
    {
        scatteringSourceExtra_ += sigmaFromTo[energyI][energyJ] * flux_[energyJ];
    }
}

