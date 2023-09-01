const scalarField& alphal(liquid->primitiveField());
const scalarField& rhol(liquid->thermo().rho().primitiveField());

if (vapour != nullptr)//if 2 phases
{
    const scalarField& alphav(vapour->primitiveField());
    const scalarField& rhov(vapour->thermo().rho().primitiveField());
    
    bafflelessAlpha.primitiveFieldRef() = alphal+alphav;

    bafflelessRhoCool.primitiveFieldRef() = 
        (alphal*rhol+alphav*rhov)/bafflelessAlpha.primitiveField();
    
    //- Should I weigh coolant temperature by mass or volume, in two-phase
    //  scenarios? I don't know, I'll use the liquid temperature for now
    bafflelessTCool.primitiveFieldRef() = 
        liquid->thermo().T().primitiveField();

    //- Baffleless should be the superficial velocity, not the real one!
    //  Thus, multiply by the void fraction (no vapour here so it's just 
    //  alphal). 
    bafflelessU.primitiveFieldRef() = 
        liquid->U().primitiveField()*alphal 
    +   vapour->U().primitiveField()*alphav;
    bafflelessU.correctBoundaryConditions();    

    //- Weighted average of the two phases
    bafflelessAlphat.primitiveFieldRef() = 
        liquid->turbulence().alphat()().primitiveField()*liquid->normalized() 
    +   vapour->turbulence().alphat()().primitiveField()*vapour->normalized();
    bafflelessAlphat.correctBoundaryConditions();

    bafflelessMu.primitiveFieldRef() = 
        liquid->thermo().mu()().primitiveField()*liquid->normalized() 
    +   vapour->thermo().mu()().primitiveField()*vapour->normalized();
    bafflelessMu.correctBoundaryConditions();       
}
else
{
    bafflelessAlpha.primitiveFieldRef() = alphal;

    bafflelessRhoCool.primitiveFieldRef() = rhol;

    bafflelessTCool.primitiveFieldRef() = 
        liquid->thermo().T().primitiveField();

    bafflelessAlphat.primitiveFieldRef() = 
        liquid->turbulence().alphat()().primitiveField();
    bafflelessAlphat.correctBoundaryConditions();

    bafflelessMu.primitiveFieldRef() = 
        liquid->thermo().mu()().primitiveField();
    bafflelessMu.correctBoundaryConditions();

    //- Baffleless should be the superficial velocity, not the real one!
    //  Thus, multiply by the void fraction (no vapour here so it's just 
    //  alphal)
    bafflelessU.primitiveFieldRef() = liquid->U().primitiveField()*alphal;
    bafflelessU.correctBoundaryConditions();
}

bafflelessAlpha.correctBoundaryConditions();
bafflelessRhoCool.correctBoundaryConditions();
bafflelessTCool.correctBoundaryConditions();

//- The fuel temperature is the liquid temperature if doing liquidFuel
bafflelessTFuelAv.primitiveFieldRef() = 
    (multiphysics.liquidFuel()) ?   
    liquid->thermo().T().primitiveField() :
    TFuelAv.primitiveField();
bafflelessTFuelAv.correctBoundaryConditions();
bafflelessTCladAv.primitiveFieldRef() = TCladAv.primitiveField();
bafflelessTCladAv.correctBoundaryConditions();
bafflelessTStruct.primitiveFieldRef() = TStruct.primitiveField();
bafflelessTStruct.correctBoundaryConditions();

//- Update phi on the NB registry if doing pointKinetics with GEMs (which is 
//  the only condition in which the ptr is set)
if (phiFMNBReg != nullptr)
{
    *phiFMNBReg = phiFMReg;
}
