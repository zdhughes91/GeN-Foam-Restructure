bool writeOutput(runTime.controlDict().lookupOrDefault<bool>("writeDat", true));

OFstream out(runTime.path()+"/GeN-Foam.dat");

if (writeOutput)
{
    out << "time(s);";
    out << "keff(-);";
    out << "power(W);";
    forAll(neutronics.fluxes(), energyI)
    {
        out << "flux" << energyI;
        if (energyI != neutronics.fluxes().size()-1)
        {
            out << ";";
        }
        else
        {
            out << "(m-2s-1);";
        }
    }
    out << "TFuel_Max;Avg;Min;";
    out << "TCladding_Max;Avg;Min;";
    out << "TCoolant_Max;Avg;Min(K)";
    out << endl << endl;
    
}

//- Initialize ptr to nuclearFuelPin, if any
const powerModel* const pinIO
(
    (fluidMesh.foundObject<volScalarField>("powerDensity.nuclearFuelPin")) ?
    fluidMesh.lookupClass<const powerModel>()["powerModel"] : nullptr
);
