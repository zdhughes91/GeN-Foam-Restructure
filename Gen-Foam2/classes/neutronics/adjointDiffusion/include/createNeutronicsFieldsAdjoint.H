if (liquidFuel_)
{
    UPtr_.reset
    (
        new volVectorField
        (
            IOobject
            (
                "U",
                mesh.time().timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedVector("", dimVelocity, vector::zero),
            zeroGradientFvPatchVectorField::typeName
        )
    );
    alphaPtr_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "alpha",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("", dimless, 1.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    alphatPtr_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "alphat",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("", dimMass/dimLength/dimTime, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    muPtr_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "mu",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("", dimMass/dimLength/dimTime, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    phiPtr_.reset
    (
        new surfaceScalarField
        (
            IOobject
            (
                "phi",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fvc::flux(UPtr_())
        )
    );
    diffCoeffPrecPtr_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "diffCoeffPrec",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("", dimArea/dimTime, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}

forAll(flux_,fluxI)
{
    flux_.set
    (
        fluxI,
        new volScalarField
        (
            IOobject
            (
                "adjointFlux"+Foam::name(fluxI),
                mesh_.time().timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            defaultFlux_
        )
    );
}

forAll(fluxStar_,fluxI)
{
        fluxStar_.set
        (
            fluxI,
            new volScalarField
            (
                IOobject
                (
                    "adjointFluxStar"+Foam::name(fluxI),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                defaultFlux_
            )
        );
}

oneGroupFlux_ *= 0.0;
forAll(flux_,energyI)
{
    oneGroupFlux_+=flux_[energyI];
}


forAll(precStar_,precI)
{
    precStar_.set
    (
        precI,
        new volScalarField
        (
            IOobject
            (
                "adjointPrecStar"+Foam::name(precI),
                mesh_.time().timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            defaultPrec_
        )
    );
}


forAll(prec_,precI)
{
    prec_.set
    (
        precI,
        new volScalarField
        (
            IOobject
            (
                "adjointPrec"+Foam::name(precI),
                mesh_.time().timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            defaultPrec_
        )
    );
}






