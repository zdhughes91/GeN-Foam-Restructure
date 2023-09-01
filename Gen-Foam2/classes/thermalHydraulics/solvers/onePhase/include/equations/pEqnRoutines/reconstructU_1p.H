//- The ifs are just to avoid useless field reconstructions
//  for integer values of the porousInterfaceSharpness_ factor
if (porousInterfaceSharpness_ != 0.0)
{
    if (porousInterfaceSharpness_ == 1.0)
        U = fvc::reconstruct(phi_)/alpha;
    else
        U = 
            porousInterfaceSharpness_*
            (
                fvc::reconstruct(phi_)/alpha
            )
        +   (1-porousInterfaceSharpness_)*
            (   
                fvc::reconstruct(fluid_.phi())
            );
}
else
    U = fvc::reconstruct(fluid_.phi());

fluid_.correctAlphaRhoMagU();