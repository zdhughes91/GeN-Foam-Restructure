
List<PtrList<fvMatrix<scalar> > > neutroEqList(xs_.energyGroups());

forAll(neutroEqList, energyI)
{
    neutroEqList[energyI].setSize(discreteDirections_);
    forAll(neutroEqList[energyI], directionI)
    {
        neutroEqList[energyI].set
        (
            directionI,
            new fvMatrix<scalar>
            (
                fvm::ddt((1-eigenvalueNeutronics_)*IV[energyI],angularFlux_[energyI][directionI])
                        + fvm::div(facePhis_[directionI], angularFlux_[energyI][directionI],"div(facePhi_,angularFlux_)")
                        + fvm::Sp((sigmaDisapp[energyI]+sigmaFromToP[0][energyI][energyI]), angularFlux_[energyI][directionI])
                        
            )
        );
    }
}

