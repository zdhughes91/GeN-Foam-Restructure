if (mesh_.time().outputTime())
{
    volSymmTensorField sigma
    (
        IOobject
        (
            "sigma",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        rho_*sigmaD_
    );

    sigma = sigma - I*(rho_*threeKalpha_*(TStruct_-TStructRef_));

    volScalarField sigmaEq
    (
        IOobject
        (
            "sigmaEq",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sqrt((3.0/2.0)*magSqr(dev(sigma)))
    );

    Info<< "Max sigmaEq = " << max(sigmaEq).value()
        << endl;
}
