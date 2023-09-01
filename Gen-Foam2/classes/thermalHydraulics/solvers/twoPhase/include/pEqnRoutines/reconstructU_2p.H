phi1.relax();
phi2.relax();

//- The ifs are just to avoid useless field reconstructions
//  for integer values of the porousInterfaceSharpness_ factor
if (porousInterfaceSharpness_ != 0.0)
{
    if (porousInterfaceSharpness_ == 1.0)
    {
        U1 = 
            fvc::reconstruct(phi1*movingAlphaf)/
            movingAlpha_;
        U2 = 
            fvc::reconstruct(phi2*movingAlphaf)/
            movingAlpha_;
    }
    else
    {
        U1 = 
            porousInterfaceSharpness_*
            (
                fvc::reconstruct
                (
                    phi1*movingAlphaf
                )/
                movingAlpha_
            )
        +   (1-porousInterfaceSharpness_)*
            (   
                fvc::reconstruct(phi1)
            );
        U2 = 
            porousInterfaceSharpness_*
            (
                fvc::reconstruct
                (
                    phi2*movingAlphaf
                )/
                movingAlpha_
            )
        +   (1-porousInterfaceSharpness_)*
            (   
                fvc::reconstruct(phi2)
            );
        }
}
else
{
    U1 = fvc::reconstruct(phi1);
    U2 = fvc::reconstruct(phi2);
}

U1.correctBoundaryConditions();
fvOptions_.correct(U1);
fluid1_.magU() = mag(U1);
fluid1_.correctAlphaRhoMagU();

U2.correctBoundaryConditions();
fvOptions_.correct(U2);
fluid2_.magU() = mag(U2);
fluid2_.correctAlphaRhoMagU();

//- Update mixture velocity
U_ = fluid1_*U1+fluid2_*U2;
