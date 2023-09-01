bool adjustTimeStep
(
    runTime.controlDict().lookupOrDefault("adjustTimeStep", false)
);

//- First, the thermalHydraulic class adjusts the global time step based on
//  possible Co number limitations. The adjustTimeStep flag is read within the
//  fuction itself and runTime.deltaT is set (if applicable) within it
if (multiphysics.solveFluidMechanics() or multiphysics.solveEnergy())
{
    thermalHydraulics.adjustTimeStep();
}

//- This bit of code is executed only the first time setDeltaT.H is included
#ifndef FIRSTDELTAT
#define FIRSTDELTAT

scalar pTotOld(neutronics.power());

//- Adjusting the time step based on neutronics limitations on the first time
//  step is useless because information relative to the "real" pTotOld is
//  unavailable. Thus, for conservativeness, only choose the smallest time
//  step between that specified by the thermalHydraulics (if any) and the
//  initial one 
if (adjustTimeStep)
{
    runTime.setDeltaT
    (
        min
        (
            runTime.deltaTValue(), 
            runTime.controlDict().get<scalar>("deltaT")
        )
    );
    Info << "initialDeltaT = " << runTime.deltaTValue() << endl;
}

#else //- Otherwise, if FIRSTDELTAT is set, this is executed

//- Then, the neutronics adjusts it again, if needed. The coding of this
//  section isn't amazing, but it is due to the fact that the thermalHydraulic
//  class was developed independently of GeN-Foam and then shoved in it. That
//  is why the thermalHydraulics adjusts the time step via an ad-hoc made
//  fuction while the neutronics does it this way
if 
(
    adjustTimeStep 
and (runTime.timeIndex()-runTime.startTimeIndex()) > 1
)
{
    scalar maxDeltaT =
        runTime.controlDict().lookupOrDefault<scalar>("maxDeltaT", GREAT);

    scalar maxPowerVariation =
        runTime.controlDict().lookupOrDefault<scalar>
        (
            "maxPowerVariation", 
            0.025
        );

    if (multiphysics.solveNeutronics())
    {
        scalar pTot = neutronics.power();
        
        scalar powerVariation = mag((pTot - pTotOld) / (pTotOld + SMALL));
        
        pTotOld = pTot;

        scalar maxDeltaTNeutroFact = 
            mag(maxPowerVariation/(powerVariation + SMALL));

        scalar deltaTNeutroFact = 
            min(min(maxDeltaTNeutroFact, 1.0 + 0.1*maxDeltaTNeutroFact), 1.2);

        //- If solveFluidMechanics is true, the timeStep could have already
        //  been adjusted by the thermalHyraulics, in which case set the new
        //  timeStep only if the neutronics requires an even shorter time-step.
        runTime.setDeltaT
        (
            min
            (
                deltaTNeutroFact*runTime.deltaTValue(),
                (multiphysics.solveFluidMechanics()) ?
                runTime.deltaTValue() : maxDeltaT
            )
        );
    }

    Info << "deltaT = " << runTime.deltaTValue() << endl;
}
    
#endif
