if (writeOutput)
{
    //- Now write stuff
    Info<< "Writing output files" << endl << endl;

    out.precision(10);
    out << runTime.time().value() << ";";

    out.precision(6);
    out << neutronics.keff() << ";";

    out.precision(5);
    out << neutronics.power() << ";";

    out.precision(5);
    forAll(neutronics.fluxes(), energyI)
    {
        out << fvc::domainIntegrate(neutronics.fluxes()[energyI]).value() << ";";
    }

    if (fluidMesh.foundObject<volScalarField>("powerDensity.nuclearFuelPin"))
    {        
        out << pinIO->get<scalar>("Tfmax") << ";";
        out << pinIO->get<scalar>("Tfavav") << ";";
        out << pinIO->get<scalar>("Tfmin") << ";";
        out << pinIO->get<scalar>("Tcmax") << ";";
        out << pinIO->get<scalar>("Tcavav") << ";";
        out << pinIO->get<scalar>("Tcmin") << ";";
    }
    else if (multiphysics.liquidFuel())
    {
        out << "0;";
        out << "0;";
        out << "0;";
        out << "0;";
        out << "0;";
        out << "0;";
    }
    else
    {
        out << "0;";
        out << "0;";
        out << "0;";
        out << "0;";
        out << "0;";
        out << "0;";
    }

    //- Coolant
    out << max(liquid->T()).value() << ";";
    out << liquid->T().weightedAverage(fluidMesh.V()).value() << ";";
    out << min(liquid->T()).value() << ";";

    out << endl;
}
