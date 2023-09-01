{
    //- Field references
    volScalarField& alpha(fluid_);
    volScalarField& he(fluid_.thermo().he());
    volScalarField& T(fluid_.thermo().T());
    const volVectorField& U(fluid_.U());
    const volVectorField& UDarcy(UDarcy_());
    const surfaceScalarField& phiDarcy(phiDarcy_());
    const volScalarField& contErr(fluid_.contErr());
    volScalarField rho(fluid_.thermo().rho());
    tmp<volScalarField> cp(fluid_.thermo().Cp());

    //- Kinetic energy per unit volume
    volScalarField K(0.5*U&U);

    volTensorField alphaEff
    (
        structure_.tortuosityPtr().valid() ?
        structure_.tortuosity()
        *fluid_.thermo().alphaEff(fluid_.turbulence().alphat())()
        :
        tensor::I
        *fluid_.thermo().alphaEff(fluid_.turbulence().alphat())()
    );

    //- Energy equations
    fvScalarMatrix EEqn
    (
            fvm::ddt(alpha, rho, he) + fvm::div(phiDarcy, he, "div(phiDarcy,he)")
        -   fvm::SuSp(contErr, he)
        +   fvc::ddt(alpha, rho, K) + fvc::div(phiDarcy, K, "div(phiDarcy,K)")
        -   contErr*K
        +   (
                he.name() == "e"
                ?   fvc::div
                    (
                        fvc::absolute(phiDarcy/fvc::interpolate(rho), UDarcy),
                        p_,
                        "div(phiv,p)"
                    )
                :   -alpha*fvc::ddt(p_)
            )
        -   fvm::laplacian
            (
                fvc::interpolate(alpha*alphaEff),
                he
            )
        ==
        fvOptions_(alpha, rho, he)
    );

    //- Structure heat sources
    //- Correct surface temperature of the power-generating subStructures
    //  of the structures and add source/sink terms to the fluids energy
    //  eqns
    const volScalarField& Hs(FSPair_.htc());

    structure_.correct
    (
        Hs*T,
        Hs
    );
    EEqn -= structure_.linearizedSemiImplicitHeatSource(fluid_, Hs)();

    //- Add additional sources deposited directly in the fluids, if any
    EEqn -= alpha * fluid_.powerDensity();


    //- Relax and solve equations, re-limit enthalpies, update fluid
    //  temperatures
    EEqn.relax();
    fvOptions_.constrain(EEqn);
    EEqn.solve();
    fvOptions_.correct(he);
    fluid_.thermo().correct();

    //- If the incompressible treatment is activated we have to update
    //- therm.rho() to make sure that the neutronics equations have acess
    //- to the density feedback.
    if(incompressibleTreatment_)
    {
        rhok_ = 1.0 - beta_()*(fluid_.thermo().T()-Tref_());
        rhok_.correctBoundaryConditions();
        fluid_.thermo().rho() = fixedRho_*rhok_;
        fluid_.thermo().rho().correctBoundaryConditions();
    }

    Info<< "T (avg min max) = "
        << T.weightedAverage(mesh_.V()).value()
        << " " << min(T).value()
        << " " << max(T).value()
        << " K" << endl;

    //- Energy conservation
    /*
    fluid_.heErr() =
        (
                fvc::ddt(alpha, rho, he)
            +   fvc::div(fluid_.alphaRhoPhi(), he, "div(alphaRhoPhi,he)")
            -   contErr*he
            -   fvc::laplacian
                (
                    fvc::interpolate(alpha)
                   *fvc::interpolate(alphaEff),
                    he
                )
            +   fvc::ddt(alpha, rho, K)
            +   fvc::div(fluid_.alphaRhoPhi(), K, "div(alphaRhoPhi,K)")
            -   contErr*K
            +   alpha*dpdt_
            -   structure_.explicitHeatSource(fluid_, Hs)
            -   alpha*rho*(U&g_)
        )/(alpha*rho*he);
    fluid_.heErr().correctBoundaryConditions();

    Info<< "Fluid energy error (avg min max) = "
        << fluid_.heErr().weightedAverage(mesh_.V()).value()
        << " " << min(fluid_.heErr()).value()
        << " " << max(fluid_.heErr()).value()
        << " 1/s" << endl;
    */
}
