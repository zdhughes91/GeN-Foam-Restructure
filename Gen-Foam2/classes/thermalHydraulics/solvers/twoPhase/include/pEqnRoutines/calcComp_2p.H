//- Construct the compressibility parts of the pressure equation
tmp<fvScalarMatrix> p_rghCorr(correction(fvm::ddt(p_rgh_)));    // <- might do 
                                                                // a lookup for
                                                                // rho
tmp<fvScalarMatrix> pEqnComp1;
tmp<fvScalarMatrix> pEqnComp2;
if (FFPair_.phaseChange())
{
    pEqnComp1 =
        (
            fvc::ddt(alpha1, rho1)
        +   fvc::div(fluid1_.alphaRhoPhi())
        -   fvc::Sp
            (
                fvc::ddt(alpha1) 
            +   fvc::div(fluid1_.alphaPhi()), 
                rho1
            )
        )/rho1
    +   (alpha1*fluid1_.thermo().psi()/rho1)*p_rghCorr()
    +   FFPair_.dmdt()/rho1;

    pEqnComp2 =
        (
            fvc::ddt(alpha2, rho2)
        +   fvc::div(fluid2_.alphaRhoPhi()) 
        -   fvc::Sp
            (
                fvc::ddt(alpha2) 
            +   fvc::div(fluid2_.alphaPhi()), rho2
            )
        )/rho2
    +   (alpha2*fluid2_.thermo().psi()/rho2)*p_rghCorr()
    -   FFPair_.dmdt()/rho2;
}
else
{
    pEqnComp1 =
        (
            fvc::ddt(alpha1, rho1)
        +   fvc::div(fluid1_.alphaRhoPhi())
        -   fvc::Sp
            (
                fvc::ddt(alpha1) 
            +   fvc::div(fluid1_.alphaPhi()), 
                rho1
            )
        )/rho1
    +   (alpha1*fluid1_.thermo().psi()/rho1)*p_rghCorr();

    pEqnComp2 =
        (
            fvc::ddt(alpha2, rho2)
        +   fvc::div(fluid2_.alphaRhoPhi()) 
        -   fvc::Sp
            (
                fvc::ddt(alpha2) 
            +   fvc::div(fluid2_.alphaPhi()), rho2
            )
        )/rho2
    +   (alpha2*fluid2_.thermo().psi()/rho2)*p_rghCorr();
}
  