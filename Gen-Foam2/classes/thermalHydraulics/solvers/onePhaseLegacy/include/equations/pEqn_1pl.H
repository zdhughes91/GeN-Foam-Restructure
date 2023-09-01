
//- Compute the total volume here
dimensionedScalar VTot("", dimVol, gSum(mesh_.V()));
surfaceScalarField alphaf(fvc::interpolate(alpha));

//- Check if the density depends on the pressure
dimensionedScalar compressibility(fvc::domainIntegrate(psi));
bool compressible = (compressibility.value()/VTot.value() >= 1e-10);

//- Construct 1/A field and interpolate to faces
volScalarField rAU("rAUDarcy", 1.0/UEqn.A());

//- Pressure corrector loop.
while (pimple_.correct())
{
    this->correctContErr();
    
    if(!incompressibleTreatment_)
    {
        fixedRho_ = fluid_.thermo().rho();
    }

    surfaceScalarField alphaRhorAUf(fvc::interpolate(alpha*fixedRho_*rAU));

    volVectorField HbyA("HbyA", UDarcy);
    HbyA = rAU*(UEqn.H());
    
    //- Note: rho is a lagged density for compressible flows. It will be 
    //- updated throughout this iteration scheme. The latest value of rho 
    //- is accessible through the thermo_ member variable. Also, rhok=1 for 
    //- compressible fluids.
    surfaceScalarField phig(-alphaRhorAUf*ghf_*fvc::snGrad(fixedRho_*rhok_)*mesh_.magSf());
    
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (
            fvc::flux(fixedRho_*HbyA)
            + alphaRhorAUf*fvc::ddtCorr(fixedRho_, UDarcy, phiDarcy)
        )
        + phig
    );

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p_rgh_, fixedRho_, UDarcy, phiHbyA, alphaRhorAUf);

    fvScalarMatrix p_rghDDtEqn
    (
        (incompressibleTreatment_)*(fvc::ddt(alpha,fixedRho_) + psi*correction(fvm::ddt(alpha,p_rgh_)))
        + fvc::div(phiHbyA)
    );

    //- Cache p prior to solve for density update
    volScalarField p_rgh_0(p_rgh_);

    //- Iterate over the pressure equation to correct for non-orthogonality
    while (pimple_.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            p_rghDDtEqn
            - fvm::laplacian(alphaRhorAUf, p_rgh_)
        );

        pEqn.setReference
        (
            pRefCell_, 
            !incompressibleTreatment_ ? getRefCellValue(p_rgh_, pRefCell_) : pRefValue_
        );

        solve
        (
            pEqn,
            mesh_.solver(p_rgh_.select(pimple_.finalInnerIter()))
        );

        if (pimple_.finalNonOrthogonalIter())
        {
            phiDarcy = phiHbyA + pEqn.flux();

            //- The following fields are updated because they provide the 
            //- coupling to the turbulence model and other sub-solvers.

            fluid_.alphaRhoPhi() = phiDarcy;
            fluid_.alphaPhi() = fluid_.alphaRhoPhi()/fvc::interpolate(fixedRho_);
            fluid_.phi() = fluid_.alphaPhi() /alphaf;

            p_rgh_.relax();

            UDarcy = HbyA + rAU*alpha*fvc::reconstruct((phig + pEqn.flux())/alphaRhorAUf);

            UDarcy.correctBoundaryConditions();
            
            U = UDarcy/alpha;

            U.correctBoundaryConditions();

            fvOptions_.correct(UDarcy);

            fluid_.magU() = mag(U);

            fluid_.correctAlphaRhoMagU();
        }
    }

    //- Update and limit the static pressure
    p_ = max(p_rgh_ + fixedRho_*rhok_*gh_, pMin_);

    if(!incompressibleTreatment_)
    {
        //- Update densities from change in p_rgh
        fluid_.thermo().rho() += psi*(p_rgh_ - p_rgh_0);

        solve(fvm::ddt(alpha,fixedRho_) + fvc::div(phiDarcy));
        
        if(p_rgh_.needReference()
           and pimple_.corrPISO() == pimple_.nCorrPISO() 
           and compressible)
        {
            dimensionedScalar deltaMass
            (
                initialFluidMass_ - fvc::domainIntegrate(alpha*fluid_.thermo().rho())
            );

            Info<< "Compressibile relative mass defect = " 
            << (deltaMass.value()/initialFluidMass_).value() << endl;
            p_ += deltaMass/compressibility;

            //- Re-limit after correction
            p_ = max(p_, pMin_);

            fixedRho_ = fluid_.thermo().rho();
        }
    }

    //- Correct p_rgh for consistency with p and the updated densities
    p_rgh_ = p_ - fixedRho_*rhok_*gh_;
    p_rgh_.correctBoundaryConditions();

    //- Computing continuity errors. If the incompressible treatment is chosen 
    //- the error is the residual in the continuity equation. Otherwise it is 
    //- the deviation between the predicted and computed fluid masses.
    if(!incompressibleTreatment_)
    {
        dimensionedScalar totalMass = fvc::domainIntegrate(alpha*fixedRho_);

        scalar sumLocalContErr =
        (
            fvc::domainIntegrate(alpha*mag(fixedRho_ - fluid_.thermo().rho()))/totalMass
        ).value();
    
        scalar globalContErr =
        (
            fvc::domainIntegrate(alpha*(fixedRho_ - fluid_.thermo().rho()))/totalMass
        ).value();
    
        cumulativeContErr_ += globalContErr;
    
        Info<< "time step continuity errors (" << mesh_.name() << ")"
            << ": sum local = " << sumLocalContErr
            << ", global = " << globalContErr
            << ", cumulative = " << cumulativeContErr_
            << endl;
    }
    else
    {
        volScalarField contErr(fvc::div(phiDarcy));

        scalar sumLocalContErr = runTime_.deltaTValue()*
                mag(contErr)().weightedAverage(mesh_.V()).value();
        
        scalar globalContErr = runTime_.deltaTValue()*
                contErr.weightedAverage(mesh_.V()).value();
        
        cumulativeContErr_ += globalContErr;

        Info<< "time step continuity errors : sum local = " << sumLocalContErr
        << ", global = " << globalContErr
        << ", cumulative = " << cumulativeContErr_
        << endl;
    }

}


