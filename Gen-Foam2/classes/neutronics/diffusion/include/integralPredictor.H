
{

    Info << "Accelerating neutronic solution: integral predictor" << endl ;

    oneGroupFlux_ *= 0.0;
    forAll(fluxStar_, energyI)
    {
           oneGroupFlux_ +=fluxStar_[energyI];
    }
    dimensionedScalar integralFlux = fvc::domainIntegrate(oneGroupFlux_);

    volScalarField oneGroupFluxOld("",oneGroupFlux_);
    oneGroupFluxOld *= 0.0;
    forAll(fluxStar_, energyI)
    {
           oneGroupFluxOld +=fluxStar_[energyI].oldTime();
    }

    dimensionedScalar integralFluxOld = fvc::domainIntegrate(oneGroupFluxOld);

    dimensionedScalar temporalTermInt("",(dimless/dimTime),scalar(0.0));
    dimensionedScalar temporalTerm("",(dimless/dimLength),scalar(0.0));

    dimensionedScalar leakageTermInt("",(dimless/dimTime),scalar(0.0));
    dimensionedScalar leakageTerm("",(dimless/dimLength),scalar(0.0));

    dimensionedScalar disappearenceTermInt("",(dimless/dimTime),scalar(0.0));
    dimensionedScalar disappearenceTerm("",(dimless/dimLength),scalar(0.0));

    dimensionedScalar sourceTermInt("",(dimless/dimTime),scalar(0.0));
    dimensionedScalar sourceTerm("",(dimless/dimLength),scalar(0.0));

    dimensionedScalar scatteringSourceTermInt("",(dimless/dimTime),scalar(0.0));
    dimensionedScalar scatteringSourceTerm("",(dimless/dimLength),scalar(0.0));

    for(label energyI=0; energyI<xs_.energyGroups() ; energyI++)
    {

        Dalbedo_ = D[energyI];

        temporalTermInt += fvc::domainIntegrate(IV[energyI]*fluxStar_[energyI]/discFactor[energyI]/mesh_.time().deltaT());

        leakageTermInt += fvc::domainIntegrate(fvc::laplacian(D[energyI]/discFactor[energyI], fluxStar_[energyI]));

        disappearenceTermInt += fvc::domainIntegrate(sigmaDisapp[energyI]*fluxStar_[energyI]/discFactor[energyI]);

        #include "calcScatteringSource.H"
        scatteringSourceTermInt += fvc::domainIntegrate(scatteringSourceExtra_);
    }

    sourceTermInt = fvc::domainIntegrate(neutroSource_/keff_*(1.0-BetaTot));

    dimensionedScalar precursorTermInt = fvc::domainIntegrate(delayedNeutroSource_);

    temporalTerm = temporalTermInt/integralFlux;

    leakageTerm = leakageTermInt/integralFlux;

    disappearenceTerm = disappearenceTermInt/integralFlux;

    sourceTerm = sourceTermInt/integralFlux;

    scatteringSourceTerm = scatteringSourceTermInt/integralFlux;

    dimensionedScalar precursorTerm = precursorTermInt;

    Info<< "Predict new integral flux " << endl;
    dimensionedScalar newIntegralFlux = (precursorTerm + integralFluxOld*temporalTerm)/(temporalTerm-leakageTerm+disappearenceTerm-sourceTerm-scatteringSourceTerm);

    Info<< "Old integral flux: " << integralFlux.value() << endl;
    Info<< "New predicted integral flux: " << newIntegralFlux.value() << endl;

    Info<< "Reconstruct fluxes " << endl;
    forAll(fluxStar_, energyI)
    {
        Dalbedo_ = D[energyI];
        fluxStarAlbedo_ *= 0.0;
        fluxStar_[energyI] *= newIntegralFlux/integralFluxOld;
        fluxStar_[energyI].correctBoundaryConditions();
        flux_[energyI]=fluxStar_[energyI]/discFactor[energyI];
        flux_[energyI].correctBoundaryConditions();
    }

    #include "initializeNeutroSource.H"
}

