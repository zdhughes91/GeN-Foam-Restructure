//- Brief note on use of the const keyword with pointers:
/*
    Type* ptr;  // ptr CAN be reassigned and points to an object of type Type
                // that CAN be modified via ptr

    const Type* ptr; // ptr CAN be reassigned and points to an object of type
                // Type that CANNOT be modified via ptr

    Type* const ptr; // ptr CANNOT be reassigned and points to an object of 
                // type Type that CAN be modified via ptr

    const Type* const ptr // ptr CANNOT be reassigned and points to an object 
                // of type Type that CANNOT be modified via ptr
*/

//- Determine if thermalHydraulics is onePhase or twoPhase by counting fluid
//  objects
HashTable<fluid*> fluids(fluidMesh.lookupClass<fluid>());
int nFluids(fluids.size());

//- Non-const pointers to non-const fluid objects. The pointers are non-const
//  as they need to be set here. If only one fluid is present, no big deal. If
//  two are present, determine which is the vapour and which is the liquid.
//  This is based by comparing enthalpies of formation as specified in the
//  thermophysicalProperties dictionaries.
//  This is necessary to correctly determine the coupling fields, which could
//  be an avarage of 2 fields in the 2-phase scenarios
fluid* liquid(nullptr);
fluid* vapour(nullptr);
if (nFluids == 1)
{
    liquid = fluids[fluids.toc()[0]];
}
else
{   
    if (multiphysics.liquidFuel())
    {
        Info << "WARNING: " << endl
        << "In two-phase simulations with liquid fuel, all the "  << endl
        << "powerDensity in neutronics goes to anything that is "  << endl
        << "liquid. You are supposed to have one liquid and one gas"  << endl
        << "otherwise, power will be counted twice" << endl; 
    }
    
    fluid* const fluid0(fluids[fluids.toc()[0]]);
    fluid* const fluid1(fluids[fluids.toc()[1]]);

    if (fluid0->isLiquid() and fluid1->isGas())
    {
        liquid = fluid0;
        vapour = fluid1;
    }
    else if (fluid1->isLiquid() and fluid0->isGas())
    {
        liquid = fluid1;
        vapour = fluid0;
    }
    else    //- Default case for liquid-liquid or gas-gas systems (no clue if 
            //  GeN-Foam will ever be used for those. Pardon the hard-coded 
            //  names)
    {
        liquid = fluid0;
        vapour = fluid1;
    }
}
//- Silence "unused variable" warnings, I actually do use them and have no idea
//  why the compiler does not "see" them
(void) liquid;
(void) vapour;

//- Get ref to passive structure temperature, and fuel and cladding average
//  temperatures
const volScalarField& TStruct
(
    (*fluidMesh.findObject<volScalarField>("T.passiveStructure"))
);

const volScalarField& TFuelAv
(
    (*fluidMesh.findObject<volScalarField>("T.fuelAvForNeutronics"))
);

const volScalarField& TCladAv
(
    (*fluidMesh.findObject<volScalarField>("T.cladAvForNeutronics"))
);

//- Then, create the actual fields for the coupling. These are created over the
//  fluidMeshNB mesh, which is a baffle-less mesh that is created if the 
//  removeBaffles option is set in controlDict. This is generally required when
//  doing multi-physics on meshes that contain baffles, as baffles are known to
//  mess up mesh-to-mesh projections. If the removeBaffles options is not
//  specified, fluidMeshNB defaults to fluidMesh (see createFluidMesh.H)
volScalarField bafflelessTCool
(
    IOobject
    (
        "bafflelessTCool",
        runTime.timeName(),
        fluidMeshNB
    ),
    fluidMeshNB,
    dimensionedScalar("", dimTemperature, 0.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField bafflelessTFuelAv
(
    IOobject
    (
        "bafflelessTFuelAv",
        runTime.timeName(),
        fluidMeshNB
    ),
    fluidMeshNB,
    dimensionedScalar("",  dimTemperature, 0.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField bafflelessTCladAv
(
    IOobject
    (
        "bafflelessTCladAv",
        runTime.timeName(),
        fluidMeshNB
    ),
    fluidMeshNB,
    dimensionedScalar("",  dimTemperature, 0.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField bafflelessTStruct
(
    IOobject
    (
        "bafflelessTStruct",
        runTime.timeName(),
        fluidMeshNB
    ),
    fluidMeshNB,
    dimensionedScalar("", dimTemperature, 0.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField bafflelessAlpha
(
    IOobject
    (
        "bafflelessAlpha",
        runTime.timeName(),
        fluidMeshNB
    ),
    fluidMeshNB,
    dimensionedScalar("", dimless, 0.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField bafflelessRhoCool
(
    IOobject
    (
        "bafflelessRhoCool",
        runTime.timeName(),
        fluidMeshNB
    ),
    fluidMeshNB,
    dimensionedScalar("", dimDensity, SMALL),
    zeroGradientFvPatchScalarField::typeName
);

//- Turbulent thermal diffusivity of liquid. For some reason, thermal
//  diffusivity as returned by thermo() is always multiplied by the density and
//  is not in m2/s, rather, kg/m/s (== Pa*s). It's a standard OpenFOAM thing,
//  not a bug, yet still...
volScalarField bafflelessAlphat
(
    IOobject
    (
        "bafflelessAlphat",
        runTime.timeName(),
        fluidMeshNB
    ),
    fluidMeshNB,
    dimensionedScalar("", dimArea/dimTime*dimDensity, 0.0),
    zeroGradientFvPatchScalarField::typeName
);

//- Molecular viscosity of liquid
volScalarField bafflelessMu
(
    IOobject
    (
        "bafflelessMu",
        runTime.timeName(),
        fluidMeshNB
    ),
    fluidMeshNB,
    dimensionedScalar("", dimPressure*dimTime, 0.0),
    zeroGradientFvPatchScalarField::typeName
);

volVectorField bafflelessU
(
    IOobject
    (
        "bafflelessU",
        runTime.timeName(),
        fluidMeshNB
    ),
    fluidMeshNB,
    dimensionedVector("", dimVelocity, vector::zero),
    zeroGradientFvPatchScalarField::typeName
);

//- This fields are sort of an exception required when using pointKinetics.
//  PointKinetics by itself simply rescales all neutronic fields according to
//  the predicted power change. However, if no powerDensities exists on the
//  neutroMesh, then look for powerDensity on the fluidMesh. If it exists, copy
//  it to a bafflelessPowerDensities that will be interpolated to the
//  powerDensities on the neutroMesh of pointKinetics in the 
//  getCouplingFieldRefs function.
volScalarField* bafflelessPowerDensity = nullptr;
volScalarField* bafflelessPowerDensityToLiquid = nullptr;

//- This is only required when point kinetics model Gas Expansion Modules (i.e.
//  a safety feature specific to one SFR, the FFTF. It is not baffleless as this
//  field is not interpolated, it is required on the original fluidMesh by the
//  pointKinetics model, yet it needs to be inserted in the fluidMeshNB 
//  objectRegistry, or the pointKinetics model won't be able to look it up.
//  phiFMNBReg is the phi but checked in the the fluidMeshNB registry, while
//  phiFMReg is the phi in the fluidMesh registry (i.e., the original one).
//  both are however defined on the fluidMesh mesh anyway
surfaceScalarField* phiFMNBReg = nullptr;
const surfaceScalarField& phiFMReg
(
    (fluidMesh.lookupObject<surfaceScalarField>("phi"))
);

if (neutronics.type() == "pointKinetics")
{
    
    bafflelessPowerDensity = 
        new volScalarField
        (
            IOobject
            (
                "bafflelessPowerDensity",
                runTime.timeName(),
                fluidMeshNB
            ),
            fluidMeshNB,
            dimensionedScalar("", dimPower/dimDensity, 0.0),
            zeroGradientFvPatchScalarField::typeName
        );
    bafflelessPowerDensity->primitiveFieldRef() = 
        fluidMesh.findObject<volScalarField>("powerDensityNeutronics")->primitiveField();
    bafflelessPowerDensity->correctBoundaryConditions();

    bafflelessPowerDensityToLiquid = 
        new volScalarField
        (
            IOobject
            (
                "bafflelessPowerDensityToLiquid",
                runTime.timeName(),
                fluidMeshNB
            ),
            fluidMeshNB,
            dimensionedScalar("", dimPower/dimDensity, 0.0),
            zeroGradientFvPatchScalarField::typeName
        );
    bafflelessPowerDensityToLiquid->primitiveFieldRef() = 
        fluidMesh.findObject<volScalarField>("powerDensityNeutronicsToLiquid")->primitiveField();
    bafflelessPowerDensityToLiquid->correctBoundaryConditions();
    
    //- Phi only required by GEM modelling
    if (neutronics.get<bool>("GEM"))
    {
        phiFMNBReg =
            new surfaceScalarField
            (
                IOobject
                (
                    "phi",
                    runTime.timeName(),
                    fluidMeshNB
                ),
                fluidMesh,
                dimensionedScalar("", dimVolume/dimTime, 0.0)
            );
        *phiFMNBReg = phiFMReg;
    }
}

//- Update the above fields so that, if any neutronic classes rely on these
//  fields for further local field initialization (e.g. pointKinetics), correct
//  initial conditions are set
if (multiphysics.solveNeutronics() or multiphysics.solveThermoMechanics())
{
    #include "correctCouplingFields.H"
}

//- Set local neutronics and thermoMechanics references to fields they require
//  for coupling. These are all obtained via lookups. Additionally, the
//  neutronics->thermalHydraulic interpolator is passed as some classes might
//  require to project some fields right after construction and before the 
//  start of the simulation e.g. pointKinetics
if (multiphysics.solveEnergy())
{
    thermalHydraulics.getCouplingFieldRefs
    (
        neutroMesh,
        fluidToNeutro
    );
}
if (multiphysics.solveNeutronics())
{
    neutronics.getCouplingFieldRefs
    (
        fluidMeshNB,
        neutroToFluid
    );
}
if (multiphysics.solveThermoMechanics())
{
    thermoMechanics.getCouplingFieldRefs
    (
        fluidMeshNB,
        mechToFluid
    );
}
