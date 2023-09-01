Info << "Accelerating neutronic solution: implicit predictor (not tested for SP3)" << endl ;
{

    PtrList<volScalarField> leakages(xs_.energyGroups());
    forAll(leakages,energyI)
    {
        Dalbedo_ = D[energyI];
        fluxStarAlbedo_ = fluxStar2_[energyI];
        leakages.set
        (
            energyI,
            new volScalarField
            (
                IOobject
                (
                    "leakages"+Foam::name(energyI),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fvc::laplacian(D[energyI]/discFactor[energyI],fluxStar_[energyI]) // /fluxStar_[energyI]
            )
        );
    }

    PtrList<volScalarField> ddtFlux2(xs_.energyGroups());

    forAll(ddtFlux2,energyI)
    {
        ddtFlux2.set
        (
            energyI,
            new volScalarField
            (
                IOobject
                (
                    "ddtFlux2"+Foam::name(energyI),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fvc::ddt(IV[energyI]/discFactor[energyI],fluxStar2_[energyI])
            )
        );
    }


    forAll(fluxStar_[0].internalField(),i)
    {
        SquareMatrix<scalar> matrix(xs_.energyGroups(),xs_.energyGroups(),Foam::zero());
        List<scalar> source(xs_.energyGroups(),0.0);
        List<scalar> predictedFlux(xs_.energyGroups(),0.0);

        forAll(flux_,energyI)//initialized predicted flux
        {
            predictedFlux[energyI] = fluxStar_[energyI].internalField()[i] ;
        }

        forAll(flux_,energyI)
        {
            forAll(flux_,energyJ)
            {
                if(energyJ==energyI)
                {
                    matrix[energyI][energyJ]=
                                  IV[energyI].internalField()[i]/mesh_.time().deltaT().value()/discFactor[energyI].internalField()[i]
                                - nuSigmaEff[energyI].internalField()[i]/keff_*(1.0-BetaTot.internalField()[i])*chiPrompt[energyI].internalField()[i]/discFactor[energyI].internalField()[i]
                                + sigmaDisapp[energyI].internalField()[i]/discFactor[energyI].internalField()[i];
                                //- leakages[energyI].internalField()[i];
                }else
                {
                    matrix[energyI][energyJ]=
                                - sigmaFromTo[energyJ][energyI].internalField()[i]/discFactor[energyJ].internalField()[i]
                                - nuSigmaEff[energyJ].internalField()[i]/keff_*(1.0-BetaTot.internalField()[i])*chiPrompt[energyI].internalField()[i]/discFactor[energyJ].internalField()[i];
                }
            }

            source[energyI]=
                      IV[energyI].internalField()[i]*fluxStar_[energyI].oldTime().internalField()[i]/mesh_.time().deltaT().value()/discFactor[energyI].internalField()[i]
                    + leakages[energyI].internalField()[i]
                    + delayedNeutroSource_.internalField()[i]*chiDelayed[energyI].internalField()[i]
                    + 2*sigmaDisapp[energyI].internalField()[i]*fluxStar2_[energyI].internalField()[i]/discFactor[energyI].internalField()[i]
                    + 2*ddtFlux2[energyI].internalField()[i];
        }

        solve(predictedFlux,matrix,source);

        forAll(fluxStar_,energyI)
        {

            fluxStar_[energyI].primitiveFieldRef()[i] = predictedFlux[energyI];

        }

    }

    forAll(flux_,energyI)
    {
        Dalbedo_ = D[energyI];
        fluxStarAlbedo_ = fluxStar2_[energyI];;
        fluxStar_[energyI].correctBoundaryConditions();
        flux_[energyI]=(fluxStar_[energyI]-fluxStar2_[energyI]*2.0)/discFactor[energyI];
        flux_[energyI].correctBoundaryConditions();

    }

    #include "initializeNeutroSource.H"

}


