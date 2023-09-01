surfaceScalarField alphaf(fvc::interpolate(alpha));

//- Update fluxes
fluid_.alphaPhi() = alphaf*fluid_.phi();
fluid_.alphaRhoPhi() = 
    fvc::interpolate(rho)*fluid_.alphaPhi();

this->correctContErr();

//- Check if compressible
dimensionedScalar compressibility(fvc::domainIntegrate(fluid_.thermo().psi()));
dimensionedScalar VTot("", dimVol, gSum(mesh_.V()));
bool compressible(compressibility.value()/VTot.value() >= 1e-8);
