//- As of now, these serve no purpose
scalar FMResidual = 0.0;
scalar EResidual = 0.0;

bool solveFM(multiphysics.solveFluidMechanics());
bool solveE(multiphysics.solveEnergy());
bool solveN(multiphysics.solveNeutronics());
bool solveTM(multiphysics.solveThermoMechanics());

//- Correct flow regime map
if (solveFM or solveE)
{
    thermalHydraulics.correctModels(solveFM, solveE);
}

//- Solve fluid-mechanics
if (solveFM) 
{
    thermalHydraulics.correctFluidMechanics(FMResidual);
    if (!solveE and !solveN and !solveTM)
        Info << endl;
}

//- Solve energy
if (!multiphysics.finalIter())
{
    if (solveE) 
    {
        if (solveN)//update power denisty only if neutronics is solved for
        {
            thermalHydraulics.interpolateCouplingFields
            (
                fluidToNeutro,
                multiphysics.liquidFuel()
            );
        }
        thermalHydraulics.correctEnergy(EResidual);
        if (!solveN and !solveTM)
            Info << endl;
    }
}

//- Solve energy-neutronics-thermomechanics coupling on last outer iteration
else if (solveE or solveN or solveTM)
{
    scalar couplingResidual = 0.0;
    label couplingIter = 0;
    do
    {
        Info << "Coupling iteration " << couplingIter << endl;
     
        //- Reset as the thermoMechanics.correct(couplingResidual) and 
        //  neutronics.correct(couplingResidual) always max() it against their
        //  solution residual, meaning that with no reset, it will stay stuck
        //  at its max value (likely the one of the first coupling iteration)
        couplingResidual = 0.0;

        if (solveE)
        {
            if (solveN)//update power denisty only if neutronics is solved for
            {
                thermalHydraulics.interpolateCouplingFields
                (
                    fluidToNeutro,
                    multiphysics.liquidFuel()
                );
            }
            thermalHydraulics.correctEnergy(couplingResidual);
            if (!solveN and !solveTM)
                Info << endl;
        }
        // correctCouplingFields.H is needed only beacuse of the speacial treatment of the 
        // thermalHydraulic mesh, where we go through a set of fields on a "phantom" mesh 
        // without baffles that is needed to facilitate the mapping.
        // If it was not for that, interpolateCouplingFields of neutronics and thermalMechanics
        // could be change to directly search for coupling fields directly in thermalHydraulics
        if (solveTM or solveN) 
        {
            #include "correctCouplingFields.H"
        }
        if (solveTM)
        {
            thermoMechanics.interpolateCouplingFields(mechToFluid);
            thermoMechanics.correct(couplingResidual);
            neutronics.deformMesh(mechToNeutro,thermoMechanics.meshDisp()); 
        }       
        if (solveN)
        {
            neutronics.interpolateCouplingFields(neutroToFluid);
            neutronics.correct(couplingResidual,couplingIter);
        }
        couplingIter++;
        Info << endl;
    }   
    while
        ( 
            (couplingResidual > multiphysics.timeStepResidual())
        &&  (couplingIter < multiphysics.maxTimeStepIterations()) 
        &&  (multiphysics.tightlyCoupled())
        );
}
