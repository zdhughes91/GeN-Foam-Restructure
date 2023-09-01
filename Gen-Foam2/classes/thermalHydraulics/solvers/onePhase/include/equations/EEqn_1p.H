{
    //- Field references
    volScalarField& alpha(fluid_);
    volScalarField& he(fluid_.thermo().he());
    volScalarField& T(fluid_.thermo().T());
    const volVectorField& U(fluid_.U());
    const volScalarField& contErr(fluid_.contErr());
    volScalarField& rho(fluid_.rho());

    volScalarField alphaEff
    (
        fluid_.thermo().alphaEff(fluid_.turbulence().alphat())()
    );

    //- Energy equations
    fvScalarMatrix EEqn
    (
            fvm::ddt(alpha, rho, he) 
        +   fvm::div(fluid_.alphaRhoPhi(), he, "div(alphaRhoPhi,he)")
        +   fvm::SuSp(-contErr, he)
        -   (
                (structure_.tortuosityPtr().valid()) ?
                fvm::laplacian
                (
                    fvc::interpolate(alpha)*
                    fvc::interpolate(structure_.tortuosity()*alphaEff),
                    he
                ) :
                fvm::laplacian
                (
                    fvc::interpolate(alpha)*
                    fvc::interpolate(alphaEff),
                    he
                )
            )
        ==  
            alpha*rho*(U&g_) 
        +   fvOptions_(alpha, rho, he)
    );

    if (fluid_.thermo().dpdt())
        EEqn -= alpha*fvc::ddt(p_);

    //- Add kinetic energy contributions
    if (runTime_.controlDict().lookupOrDefault("includeKineticEnergy", true))
    {
        //- Kinetic energy per unit volume
        volScalarField K(0.5*U&U);
        EEqn +=
        (
            fvc::ddt(alpha, rho, K) 
        +   fvc::div(fluid_.alphaRhoPhi(), K, "div(alphaRhoPhi,K)")
        -   contErr*K
        );
    }

    //- Add structure heat sources (i.e. powerModels and HX)
    //const volScalarField& Hs(*htcs_["FSPair"]);
    const volScalarField& Hs(FSPair_.htc());
    /*
    volScalarField Hs
    (
        IOobject
        (
            "Hs",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("", dimPower/dimArea/dimTemperature, 0),
        zeroGradientFvPatchTensorField::typeName
    );
    */
    EEqn -= structure_.linearizedSemiImplicitHeatSource(fluid_, Hs)();

    //- Add additional sources deposited directly in the fluids, if any
    EEqn -= alpha * fluid_.powerDensity();

    //- Cache before solution
    T.storePrevIter();

    //- Relax and solve equations, re-limit enthalpies, update fluid
    //  temperatures
    EEqn.relax();
    fvOptions_.constrain(EEqn);
    EEqn.solve();
    fvOptions_.correct(he);
    he.relax();
    fluid_.thermo().correct();

    //- Correct surface temperature of the power-generating subStructures 
    structure_.correct
    (
        Hs*T.prevIter(),
        Hs
    );

    //- Correct fluid thermo placeholder fields and Prandtl;
    fluid_.kappa() = fluid_.thermo().kappa();
    fluid_.Cp() = fluid_.thermo().Cp();
    fluid_.mu() = fluid_.thermo().mu();
    fluid_.Pr() = fluid_.Cp()*fluid_.mu()/fluid_.kappa();

    Info<< "T (avg min max) = "
        << T.weightedAverage(mesh_.V()).value()
        << " " << min(T).value()
        << " " << max(T).value()
        << " K" << endl;
}
