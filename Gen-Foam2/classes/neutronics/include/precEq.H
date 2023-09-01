forAll(prec_,precI)
{
    if(liquidFuel_)
    {
        fvScalarMatrix precEq
        (
            fvm::ddt(alphaPtr_()*(1-eigenvalueNeutronics_), precStar_[precI])
        +   fvm::Sp(lambda[precI]*alphaPtr_(), precStar_[precI])
        -   neutroSource_/keff_*Beta[precI]
        +   fvm::div(phiPtr_(), precStar_[precI], "div(phi_,precStar_)")
        -   fvm::laplacian(diffCoeffPrecPtr_(), precStar_[precI])                
        );
        precEq.solve();
        prec_[precI] = precStar_[precI] * alphaPtr_(); // Precs per total volume
        prec_[precI].correctBoundaryConditions();
        /*
        //- Correct but works on discontinuous fields
        fvScalarMatrix precEq
        (
            fvm::ddt((1-eigenvalueNeutronics_), prec_[precI])
            + fvm::Sp(lambda[precI], prec_[precI])
            - neutroSource_/keff_*Beta[precI]
            + fvm::div(phi_/fvc::interpolate(porosity_), prec_[precI])
            - fvm::laplacian(diffCoeffPrec_, prec_[precI])                
        );
        precEq.solve();
        */
    }
    else
    {
        fvScalarMatrix precEq
        (
            fvm::ddt((1-eigenvalueNeutronics_), prec_[precI])
        +   fvm::Sp(lambda[precI], prec_[precI])
        -   neutroSource_/keff_*Beta[precI]
        );
        precEq.solve();
    }
}
