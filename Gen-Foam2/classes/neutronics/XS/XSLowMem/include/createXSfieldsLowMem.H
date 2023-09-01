// Nuclear data fields
Info<< "Creating nuclear data fields\n" << endl;


forAll(IV_, i)
{
    IV_.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "IV"+Foam::name(i),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("", dimensionSet(0,-1,1,0,0,0,0), 1.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}

forAll(D_, i)
{
    D_.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "D"+Foam::name(i),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,//IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("", dimensionSet(0,1,0,0,0,0,0), 1.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}

forAll(nuSigmaEff_, i)
{
    nuSigmaEff_.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "nuSigmaEff"+Foam::name(i),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("", dimensionSet(0,-1,0,0,0,0,0), 1.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}

forAll(sigmaPow_, i)
{
    sigmaPow_.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "sigmaPow"+Foam::name(i),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("", dimensionSet(1,1,-2,0,0,0,0), 1.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}

forAll(sigmaDisapp_, i)
{
    sigmaDisapp_.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "sigmaDisapp"+Foam::name(i),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar("", dimensionSet(0,-1,0,0,0,0,0), 1.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}

forAll(chiPrompt_, i)
{
    chiPrompt_.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "chiPrompt"+Foam::name(i),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("", dimensionSet(0,0,0,0,0,0,0), 1.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}

forAll(chiDelayed_, i)
{
    chiDelayed_.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "chiDelayed"+Foam::name(i),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("", dimensionSet(0,0,0,0,0,0,0), 1.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}

forAll(Beta_, i)
{
    Beta_.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "Beta"+Foam::name(i),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("", dimensionSet(0,0,0,0,0,0,0), 1.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}

forAll(lambda_, i)
{
    lambda_.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "lambda"+Foam::name(i),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("", dimensionSet(0,0,-1,0,0,0,0), 1.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}

forAll(discFactor_,energyI)
{
    discFactor_.set
    (
        energyI,
        new volScalarField
        (
            IOobject
            (
                "discFactor"+Foam::name(energyI),
                mesh_.time().timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar("", dimensionSet(0,0,0,0,0,0,0), 1.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}

Info << "done" << endl;




