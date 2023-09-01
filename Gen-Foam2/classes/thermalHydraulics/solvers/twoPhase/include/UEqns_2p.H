//- Field refs
volVectorField& U1(fluid1_.U());
volVectorField& U2(fluid2_.U());
const volScalarField& alpha1(fluid1_);
const volScalarField& alpha2(fluid2_);
volScalarField& rho1(fluid1_.rho());
volScalarField& rho2(fluid2_.rho());
fvVectorMatrix UEqn1(U1, dimMass*dimLength/dimTime/dimTime);
fvVectorMatrix UEqn2(U2, dimMass*dimLength/dimTime/dimTime);

//- Correct turbulence. What about this copy-setting fuss? Well, the standard
//  OpenFOAM turbulence models (e.g. kEpsilon) do not limit the phase fraction
//  before constructing the epsilon and k equations. This means that if a phase
//  fraction happens to be 0 somewhere, it will end up causing some null
//  entries in the k or epsilon matrices (...).
//  The solutions are: 1) copy some of the standard kEpslion models and make
//  them part of FFSEulerFoam and make the changes directly in the turbulence
//  classes; 2) limit the phase fractions before solving the turbulence models
//  and then reset them to the un-limited values cached before. Option 1)
//  is a NO-NO as I don't want to break compatibility with standard OpenFOAM
//  stuff and I want to avoid code duplication too (as if "code duplication" 
//  is a concern in a code that is 22k+ lines without even counting the 
//  underlying OpenFOAM library)
volScalarField fluid10(fluid1_);
volScalarField fluid20(fluid2_);
fluid1_.volScalarField::operator=(max(fluid1_, fluid1_.residualAlpha()));
fluid2_.volScalarField::operator=(max(fluid2_, fluid2_.residualAlpha()));
fluid1_.turbulence().correct();
fluid2_.turbulence().correct();

fluid1_.volScalarField::operator=(fluid10);
fluid2_.volScalarField::operator=(fluid20);

//- Fluid-fluid drag
const volScalarField& Kd12(FFPair_.Kd());

//- Drag between fluid1 and structure
const volTensorField& Kd1s(F1SPair_.Kd());
const volTensorField& Kd2s(F2SPair_.Kd());

{
    volScalarField nuEff1(fluid1_.turbulence().nuEff()());
    volScalarField nuEff2(fluid2_.turbulence().nuEff()());
    
    UEqn1 =
    (
        fvm::div(fluid1_.alphaRhoPhi(), U1, "div(alphaRhoPhi,U)")
        
        //- Diffusion term (i.e. fluid1_.turbulence().divDevRhoReff(U1))
        //  but re-expressed so to: 1) either allow to use a tensor
        //  diffusivity (i.e. nuEff1*structure_.tortuosity()) if the 
        //  structure is tortuous; 2) if not tortuous, still allow me to
        //  specify the div interpolation scheme (i.e. "div(alphaRhoPhiNu,U))
    +   (
            (structure_.tortuosityPtr().valid()) ? 
            (
            -   fvm::laplacian
                (
                    alpha1*rho1*nuEff1*structure_.tortuosity(), 
                    U1
                )
            -   fvc::div
                (
                    alpha1*rho1*nuEff1*structure_.tortuosity()
                &   dev2(T(fvc::grad(U1))), 
                    "div(alphaRhoPhiNu,U)"
                )
            ) :
            (
            -   fvm::laplacian(alpha1*rho1*nuEff1, U1)
            -   fvc::div
                (
                    alpha1*rho1*nuEff1*dev2(T(fvc::grad(U1))), 
                    "div(alphaRhoPhiNu,U)"
                )
            )
        )
        
        //- Implicit diagonal contribution of fluid-structure drag, 
        //  explicit off-diagonal contribution
    +   fvm::Sp((1.0/3.0)*tr(Kd1s), U1) + (dev(Kd1s)&U1)
    ==
        fvOptions_(alpha1, rho1, U1)
    );
    UEqn2 =
    (
        fvm::div(fluid2_.alphaRhoPhi(), U2, "div(alphaRhoPhi,U)")

        //- Diffusion term (i.e. fluid2_.turbulence().divDevRhoReff(U2))
        //  but re-expressed so to: 1) either allow to use a tensor
        //  diffusivity (i.e. nuEff2*structure_.tortuosity()) if the 
        //  structure is tortuous; 2) if not tortuous, still allow me to
        //  specify the div interpolation scheme (i.e. "div(alphaRhoPhiNu,U))
    +   (
            (structure_.tortuosityPtr().valid()) ? 
            (
            -   fvm::laplacian
                (
                    alpha2*rho2*nuEff2*structure_.tortuosity(), 
                    U2
                )
            -   fvc::div
                (
                    alpha2*rho2*nuEff2*structure_.tortuosity()
                &   dev2(T(fvc::grad(U2))), 
                    "div(alphaRhoPhiNu,U)"
                )
            ) :
            (
            -   fvm::laplacian(alpha2*rho2*nuEff2, U2)
            -   fvc::div
                (
                    alpha2*rho2*nuEff2*dev2(T(fvc::grad(U2))), 
                    "div(alphaRhoPhiNu,U)"
                )
            )
        )

        //- Implicit diagonal contribution of fluid-structure drag,
        //  explicit off-diagonal contribution
    +   fvm::Sp((1.0/3.0)*tr(Kd2s), U2) + (dev(Kd2s)&U2) 
    ==
        fvOptions_(alpha2, rho2, U2)
    );

    //- Add mass transfer contributions
    if (FFPair_.phaseChange())
    {
        volScalarField dmdt12(posPart(FFPair_.dmdt()));
        volScalarField dmdt21(-negPart(FFPair_.dmdt())); 
        //UEqn1 += fvm::Sp(dmdt12, U1) - dmdt21*U2;
        //UEqn2 += fvm::Sp(dmdt21, U2) - dmdt12*U1;
    }

    //- Add virtual mass forces, if applicable
    FFPair_.addVirtualMassForce(UEqn1);
    FFPair_.addVirtualMassForce(UEqn2);

    //- The implicit and implicitWithDmdt are formulated differently
    //  even when there is no mass transfer. I should write the derivations
    //  somewhere down in the source code (eventually)
    if (partialEliminationMode_ != partialEliminationMode::implicit)
    {
        UEqn1 += fvm::Sp(Kd12, U1);
        UEqn2 += fvm::Sp(Kd12, U2);
    }

    //- The time derivative is manually added to the pressure equation
    //  contribution if doing faceMomentum. Otherwise, do it here
    if (momentumMode_ != momentumMode::faceCentered) 
    {
        UEqn1 += fvm::ddt(alpha1, rho1, U1);
        UEqn2 += fvm::ddt(alpha2, rho2, U2);
    }

    scalar c
    (
        pimple_.dict().lookupOrDefault<scalar>
        (
            "continuityErrorScaleFactor", 
            1.0
        )
    );

    switch (UContErrCompensationMode_)
    {
        case contErrCompensationMode::none :
        {
            break;
        }
        case contErrCompensationMode::Su :
        {
            UEqn1 -= c*fluid1_.contErr()*U1;
            UEqn2 -= c*fluid2_.contErr()*U2;
            break;
        }
        case contErrCompensationMode::Sp :
        {
            UEqn1 -= fvm::Sp(c*fluid1_.contErr(), U1);
            UEqn2 -= fvm::Sp(c*fluid2_.contErr(), U2);
            break;
        }
        case contErrCompensationMode::SuSp :
        {
            UEqn1 += fvm::SuSp(-c*fluid1_.contErr(), U1);
            UEqn2 += fvm::SuSp(-c*fluid2_.contErr(), U2);
            break;
        }
    }

    //UEqn1.relax();
    fvOptions_.constrain(UEqn1);
    fvOptions_.correct(U1);
    //UEqn2.relax();
    fvOptions_.constrain(UEqn2);
    fvOptions_.correct(U2);
}
