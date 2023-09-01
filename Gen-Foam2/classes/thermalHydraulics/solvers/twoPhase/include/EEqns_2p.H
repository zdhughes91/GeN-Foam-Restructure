{
    //- Field references
    volScalarField& he1(fluid1_.thermo().he());
    volScalarField& he2(fluid2_.thermo().he());
    const volScalarField& T1(fluid1_.thermo().T());
    const volScalarField& T2(fluid2_.thermo().T());
    const volVectorField& U1(fluid1_.U());
    const volVectorField& U2(fluid2_.U());
    const volScalarField& contErr1(fluid1_.contErr());
    const volScalarField& contErr2(fluid2_.contErr());
    const volScalarField& Cp1(fluid1_.Cp());
    const volScalarField& Cp2(fluid2_.Cp());
    const volScalarField& rho1(fluid1_.rho());
    const volScalarField& rho2(fluid2_.rho());
    const volScalarField& iT12(FFPair_.iT());
    const volScalarField& alpha1(fluid1_);
    const volScalarField& alpha2(fluid2_);
    volScalarField he1i(fluid1_.thermo().he(p_, iT12));
    volScalarField he2i(fluid2_.thermo().he(p_, iT12));

    //- Check if structure power off criteria are met and turn off power if so
    structure_.checkPowerOff();

    //-
    volScalarField alphaEff1
    (
        fluid1_.thermo().alphaEff(fluid1_.turbulence().alphat())()
    );
    volScalarField alphaEff2
    (
        fluid2_.thermo().alphaEff(fluid2_.turbulence().alphat())()
    );

    autoPtr<volScalarField> dpdtPtr;
    if (fluid1_.thermo().dpdt() or fluid2_.thermo().dpdt())
    {
        //- Pressure power term limited to 10 bar/s
        dimensionedScalar maxDpdt
        (
            "", 
            dimPressure/dimTime, 
            pimple_.dict().lookupOrDefault<scalar>("maxDpdt", 1e6)
        );
        dpdtPtr.reset
        (
            new volScalarField
            (
                max
                (
                    min
                    (
                        fvc::ddt(p_), maxDpdt
                    ),
                    -maxDpdt
                )
            )
        );
    }
    
    //- Construct enthalpy equations
    fvScalarMatrix E1Eqn
    (
            fvm::ddt(alpha1, rho1, he1) 
        +   fvm::div(fluid1_.alphaRhoPhi(), he1, "div(alphaRhoPhi,he)")
        -   (
                (structure_.tortuosityPtr().valid()) ?
                fvm::laplacian
                (
                    fvc::interpolate(alpha1)*
                    fvc::interpolate(structure_.tortuosity()*alphaEff1),
                    he1
                ) :
                fvm::laplacian
                (
                    fvc::interpolate(alpha1)*
                    fvc::interpolate(alphaEff1),
                    he1
                )
            )
        ==  
            alpha1*rho1*(U1&g_) 
        +   fvOptions_(alpha1, rho1, he1)
    );
    fvScalarMatrix E2Eqn
    (
            fvm::ddt(alpha2, rho2, he2)
        +   fvm::div(fluid2_.alphaRhoPhi(), he2, "div(alphaRhoPhi,he)")
        -   (
                (structure_.tortuosityPtr().valid()) ?
                fvm::laplacian
                (
                    fvc::interpolate(alpha2)*
                    fvc::interpolate(structure_.tortuosity()*alphaEff2),
                    he2
                ) :
                fvm::laplacian
                (
                    fvc::interpolate(alpha2)*
                    fvc::interpolate(alphaEff2),
                    he2
                )
            )
        ==  
            alpha2*rho2*(U2&g_) 
        +   fvOptions_(alpha2, rho2, he2)
    );

    if (fluid1_.thermo().dpdt())
        E1Eqn -= alpha1*dpdtPtr();
    if (fluid2_.thermo().dpdt())
        E2Eqn -= alpha2*dpdtPtr();

    scalar c
    (
        pimple_.dict().lookupOrDefault<scalar>
        (
            "continuityErrorScaleFactor", 
            1.0
        )
    );

    switch (heContErrCompensationMode_)
    {
        case contErrCompensationMode::none :
        {
            break;
        }
        case contErrCompensationMode::Su :
        {
            E1Eqn -= c*fluid1_.contErr()*he1;
            E2Eqn -= c*fluid2_.contErr()*he2;
            break;
        }
        case contErrCompensationMode::Sp :
        {
            E1Eqn -= fvm::Sp(c*fluid1_.contErr(), he1);
            E2Eqn -= fvm::Sp(c*fluid2_.contErr(), he2);
            break;
        }
        case contErrCompensationMode::SuSp :
        {
            E1Eqn += fvm::SuSp(-c*fluid1_.contErr(), he1);
            E2Eqn += fvm::SuSp(-c*fluid2_.contErr(), he2);
            break;
        }
    }

    //- Add kinetic energy contributions
    if (runTime_.controlDict().lookupOrDefault("includeKineticEnergy", false))
    {
        volScalarField K1(0.5*U1&U1);
        E1Eqn +=
        (
            fvc::ddt(alpha1, rho1, K1) 
        +   fvc::div(fluid1_.alphaRhoPhi(), K1, "div(alphaRhoPhi,K)")
        -   contErr1*K1
        );
        volScalarField K2(0.5*U2&U2);
        E2Eqn +=
        (
            fvc::ddt(alpha2, rho2, K2) 
        +   fvc::div(fluid2_.alphaRhoPhi(), K2, "div(alphaRhoPhi,K)")
        -   contErr2*K2
        );
    }

    //- Print dmdt info
    if (FFPair_.phaseChange())
    {
        const volScalarField& dmdt(FFPair_.dmdt());
        Info<< "dmdt " << fluid1_.name() << "->" << fluid2_.name() 
            << " (avg min max) = " << dmdt.weightedAverage(mesh_.V()).value()
            << " " << min(dmdt).value()
            << " " << max(dmdt).value()
            << " kg/m3/s" << endl;
    }
    
    //- Add interfacial heat transfer contributions due to heat and mass
    //  transfer
    E1Eqn += FFPair_.heSource(E1Eqn.psi().name());
    E2Eqn += FFPair_.heSource(E2Eqn.psi().name());
    
    //- Structure heat sources
    //- Correct surface temperature of the power-generating subStructures 
    //  of the structures and add source/sink terms to the fluids energy 
    //  eqns. Also, limit fluid-structure HTCs only to regions where the fluid
    //  phase fraction is above thermoResidualAlpha
    volScalarField& H1s(F1SPair_.htc());
    //H1s *= pos1;
    volScalarField& H2s(F2SPair_.htc());
    //H2s *= pos2;
    
    E1Eqn -= structure_.linearizedSemiImplicitHeatSource(fluid1_, H1s)();
    E2Eqn -= structure_.linearizedSemiImplicitHeatSource(fluid2_, H2s)();

    //- Add additional sources deposited directly in the fluids, only if
    //  they are liquid
    if (fluid1_.isLiquid())
    {
        E1Eqn -= alpha1 * fluid1_.powerDensity();
    }
    if (fluid2_.isLiquid())
    {
        E2Eqn -= alpha2 * fluid2_.powerDensity();
    }

    //- Add stabilization sources so that the temperatures will tend to their
    //  interfacial values if the phase fraction is below thermoResidualAlpha
    if (heStabilizationMode_ == heStabilizationMode::source)
    {
        /* This is the old source stabilizer
        scalar HStab(1e10);
        if (fluid1_.thermoResidualAlpha().value() != 0.0)
        {
            volScalarField C
            (
                IOobject("", mesh_.time().timeName(), mesh_), mesh_, 
                dimensionedScalar("", dimPower/dimVolume/dimTemperature, 0.0),
                zeroGradientFvPatchScalarField::typeName
            );
            const scalar& t(fluid1_.thermoResidualAlpha().value());
            forAll(C, i)
            {
                scalar deltaAlpha(fluid1_.normalized()[i] - t);
                if (deltaAlpha < 0.0)
                    C[i] = (-neg1[i]*deltaAlpha/t)*HStab;
            }
            C.correctBoundaryConditions();
            E1Eqn -= (C*(he1/Cp1 + iT12-T1) - fvm::Sp(C/Cp1, he1));
        }
        if (fluid2_.thermoResidualAlpha().value() != 0.0)
        {
            volScalarField C
            (
                IOobject("", mesh_.time().timeName(), mesh_), mesh_, 
                dimensionedScalar("", dimPower/dimVolume/dimTemperature, 0.0),
                zeroGradientFvPatchScalarField::typeName
            );
            const scalar& t(fluid2_.thermoResidualAlpha().value());
            forAll(C, i)
            {
                scalar deltaAlpha(fluid2_.normalized()[i] - t);
                if (deltaAlpha < 0.0)
                    C[i] = (-neg2[i]*deltaAlpha/t)*HStab;
            }
            C.correctBoundaryConditions();
            E2Eqn -= (C*(he2/Cp2 + iT12-T2) - fvm::Sp(C/Cp2, he2));
        }*/
        
        //- Relaxation time
        scalar tau(1e-12);
        if (fluid1_.thermoResidualAlpha().value() != 0.0)
        {
            volScalarField C
            (
                IOobject("", mesh_.time().timeName(), mesh_), mesh_, 
                dimensionedScalar("", dimDensity/dimTime, 0.0),
                zeroGradientFvPatchScalarField::typeName
            );
            const scalar& t(fluid1_.thermoResidualAlpha().value());
            forAll(C, i)
            {
                scalar deltaAlpha(fluid1_.normalized()[i] - t);
                if (deltaAlpha < 0.0)
                    C[i] = t*rho1[i]/tau; //(-neg1[i]*deltaAlpha/t)*HStab;
            }
            C.correctBoundaryConditions();
            E1Eqn -= C*he1i-fvm::Sp(C, he1); 
                //(C*(he1/Cp1 + iT12-T1) - fvm::Sp(C/Cp1, he1));
        }
        if (fluid2_.thermoResidualAlpha().value() != 0.0)
        {
            volScalarField C 
            (
                IOobject("", mesh_.time().timeName(), mesh_), mesh_, 
                dimensionedScalar("", dimDensity/dimTime, 0.0),
                zeroGradientFvPatchScalarField::typeName
            );
            const scalar& t(fluid2_.thermoResidualAlpha().value());
            forAll(C, i)
            {
                scalar deltaAlpha(fluid2_.normalized()[i] - t);
                if (deltaAlpha < 0.0)
                    C[i] = t*rho2[i]/tau; //(-neg2[i]*deltaAlpha/t)*HStab;
            }
            C.correctBoundaryConditions();
            E2Eqn -= C*he2i-fvm::Sp(C, he2); 
                //(C*(he2/Cp2 + iT12-T2) - fvm::Sp(C/Cp2, he2));
        }
    }
    
    //- Relax and solve equations
    E1Eqn.relax();
    fvOptions_.constrain(E1Eqn);
    E1Eqn.solve();
    fvOptions_.correct(he1);
    he1.relax();
    E2Eqn.relax();
    fvOptions_.constrain(E2Eqn);
    E2Eqn.solve();
    fvOptions_.correct(he2);
    he2.relax(); 

    //- Stabilize enthalpies by cutoff if source mode not used
    if (heStabilizationMode_ == heStabilizationMode::cutoff)
    {
        const volScalarField& pos1(fluid1_.aboveThermoResidualAlpha());
        const volScalarField& pos2(fluid2_.aboveThermoResidualAlpha());
        const volScalarField& neg1(fluid1_.belowThermoResidualAlpha());
        const volScalarField& neg2(fluid2_.belowThermoResidualAlpha());
        he1 = 
            max
            (
                pos1*he1 + neg1*he1i, 
                dimensionedScalar("", dimEnergy/dimMass, 100)
            );
        he2 = 
            max
            (
                pos2*he2 + neg2*he2i, 
                dimensionedScalar("", dimEnergy/dimMass, 100)
            );
    }
    
    //- Correct temperatures based on new enthalpy and update other
    //  temperature-dpendent quantities
    fluid1_.thermo().correct();
    fluid2_.thermo().correct();

    //- Update structure temperature
    structure_.correct
    (
        H1s*T1 + H2s*T2, 
        H1s + H2s
    );

    //- Correct placeholder fields
    fluid1_.Cp() = fluid1_.thermo().Cp();
    fluid1_.kappa() = fluid1_.thermo().kappa();
    fluid1_.mu() = fluid1_.thermo().mu();
    fluid2_.Cp() = fluid2_.thermo().Cp();
    fluid2_.kappa() = fluid2_.thermo().kappa();
    fluid2_.mu() = fluid2_.thermo().mu();
    fluid1_.Pr() = Cp1*fluid1_.thermo().mu()/fluid1_.thermo().kappa();
    fluid2_.Pr() = Cp2*fluid2_.thermo().mu()/fluid2_.thermo().kappa();
    
    myOps::fieldInfo(iT12);
    myOps::fieldInfo(T1);
    myOps::fieldInfo(T2);
}
