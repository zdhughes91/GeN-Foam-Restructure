//- Update and limit the static pressure
p_ = max(p_rgh_ + rho*gh_, pMin_);

//- Correct p_ to account for compressibility effects in closed volumes.
//  Do this outside on the last correct() loop as otherwise the loop itself
//  gets as stable as the Italian government
if 
(
    p_rgh_.needReference() //- i.e. closed volume 
and (pimple_.corrPISO() == pimple_.nCorrPISO() or convergedPISO)
)
{    
    //- This commented bit was originally other compressible solvers, e.g.
    //  buoyantPimpleFoam. However, it arguably does not make much sense,
    //  as the pressure was already set in setReference, so there is no
    //  need for this correction
    /*
    if (!compressible)
    {
        p_ += dimensionedScalar
        (
            "p",
            p_.dimensions(),
            pRefValue_ - getRefCellValue(p_, pRefCell_)
        );
    }*/
    if (compressible and !fluid_.Boussinesq())
    {
        dimensionedScalar deltaMass
        (
            initialFluidMass_ - fvc::domainIntegrate(rho*fluid_)
        );
        
        Info<< "Compressibile relative mass defect = " 
        << (deltaMass.value()/initialFluidMass_).value() << endl;
        p_ += deltaMass/compressibility;
    }

    //- Re-limit after correction
    p_ = max(p_, pMin_);
}

//- Update densities from change in p_rgh
rho += fluid_.thermo().psi()*(p_rgh_ - p_rgh_0);

//- Correct p_rgh for consistency with p and the updated densities
p_rgh_ = p_ - rho*gh_;
p_rgh_.correctBoundaryConditions();
