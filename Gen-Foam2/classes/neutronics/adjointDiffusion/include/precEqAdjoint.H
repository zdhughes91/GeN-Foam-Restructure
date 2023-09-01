forAll(prec_,precI)
{
    if(liquidFuel_)
    {
        fvScalarMatrix precEq
        (
            fvm::ddt(alphaPtr_()*(1-eigenvalueNeutronics_),precStar_[precI])
            + fvm::Sp(lambda[precI]*alphaPtr_(), precStar_[precI])
            - adjointFlux_chiD_*lambda[precI]
            + fvm::div(-phiPtr_(), precStar_[precI])
            - fvm::laplacian(diffCoeffPrecPtr_(), precStar_[precI])                
        );
    	precEq.solve();
        prec_[precI] = precStar_[precI]*alphaPtr_(); //precs per total volume
    	prec_[precI].correctBoundaryConditions();
    }
    else
    {
        fvScalarMatrix precEq
        (
            fvm::ddt((1-eigenvalueNeutronics_),prec_[precI])
            + fvm::Sp(lambda[precI], prec_[precI])
            - adjointFlux_chiD_*lambda[precI]
        );
        precEq.solve();
    }
}
