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

forAll(angularFlux_,energyI)
{
    angularFlux_[energyI].setSize(discreteDirections_);
    forAll(angularFlux_[energyI],directionI)
    {
        angularFlux_[energyI].set
        (
            directionI,
            new volScalarField
            (
                IOobject
                (
                    "angularFlux_"+Foam::name(energyI)+"_"+Foam::name(directionI),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                defaultFlux_
            )
        );
    }
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
                "flux"+Foam::name(fluxI),
                mesh_.time().timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            defaultFlux_
        )
    );
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
                "precStar"+Foam::name(precI),
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
                "prec"+Foam::name(precI),
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
    delayedNeutroSource_ += xs_.lambda()[precI]*prec_[precI];
}

forAll(facePhis_, directionI)
{
    facePhis_.set
    (
        directionI,
        new surfaceScalarField("",directionVersors_[directionI] & mesh_.Sf())
    );
}    
