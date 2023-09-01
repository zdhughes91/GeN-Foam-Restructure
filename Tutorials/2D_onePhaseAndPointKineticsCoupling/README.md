# 2D 1-phase & Point-Kinetics Coupling

## Description

This is a simplified test-case for the pointKinetics neutronics model.
It consists of a 2-D square domain consisting of a single cellZone containing
fuel. A coolant (properties ~ sodium) flows from the bottom of the domain to
the top. Cross-sections have been arbitrarily chosen so to obtain a critical
system in 2 groups. They are used in an initial steady-state diffusion
calculation to obtain the power shape to be used in the point kinetics
transient calculations. In case an initial diffusion calculation is not
performed, the point kinetic solver will assume a spatially uniform power.

Three transients are simulated:

-   0.2$ reactivity insertion with the only feedback being fuel temperature,
    with the feedback set a -0.3 pcm/K (arbitrary);

-   like the trasient above, but including the feedback from an assumed
    driveline expansion. The control rod reactivity map (constant/neutroRegion/
    nuclearData) is arbitrarily set to 100 pcm/cm (more than one
    order of magnitude higher than realistic values, used here only for feature
    demonstration purposes). The absolute driveline expansion coefficient
    is assumed to be 8e-5 m/K;

- like the first transient, but with time dependent insertion of reactivity
  from Boron (removal). This reactivity profile is specified in nuclearData

In both transienst, the time steps are adjusted so to keep a power change of
~ 1% per time step.


## How to run

Use `./Allclean` to remove previous results, if any. Use `./Allrun` to
automatically run the steady state, the transient without driveline expansion
and the transient with driveline expansion. All the results of interest are
summarized in the log files, namely `steadyState.log`,
`transientNoDriveline.log` and `transientWithDriveline.log`.

To plot the results in terms of time-dependence of the power, total reactivity
and fuel temperature of the two transient cases, simply run
```bash
python plot.py transientNoDriveline.log transientWithDriveline.log
```

It requires matplotlib to plot results (tested with Python 2.7 and 3.6).


## Short discussion

The prompt jump (i.e. the height of the first power plateau) agrees with
expected theoretical results. A 0.2$ reactivity insertion should lead to an
absolute jump magnitude of ~ 1/(1-0.8) = 1.25, which agrees with the variation
of power from 1e7 to 1.25e7 W. Please note that the theoretical prompt jump
results are however derived under some assumption (e.g. 1 delayed group), and
will not be exactly reproduced by the code for obvious reasons. After the
prompt jump, the power rise is governed by precursor increase. However, shortly
after, the time scales are compatbile with fuel and driveline heat-up scales,
so that negative reactivity effects kick in. The interplay between power, fuel
temperature and driveline expansion and change in precursor population is what
ultimately shapes the power profile. In the case with driveline feedbacks, the
power quickly reaches a steady state, while in the one without, it is not
reached within the transient time scale. However it will, as the total
reactivity approaches 0 as the fuel slowly heats up.

--------------------------------------------------------------------------------

## Notes on the input files

-   constant/neutroRegion/reactorState. This file is only used to specify the
    initial power, external reactivity and, if applicable, initial precursor
    concentrations (expressed in W for consistency with the pointKinetic
    equation being solved for reactor power rathern than neutron density).
    The intial precursor concentrations can be specified via the
    precursorPowers keyword. If not found, precursor concentration are
    initialized so to be in equilibrium with the starting conditions (i.e.
    a steady state is assumed).

-   constant/neutroRegion/nuclearData. The pointKinetics-related keywords that
    can be specified are
    /*
        //- Self explanatory
        promptGenerationTime 1e-06;

        //- Delayed neutron fraction by group. The number of groups is deduced
        //  by the length of this list
        Beta
        (
            7.2315e-05
            0.000609661
            0.000471181
            0.00118907
            0.000445487
            9.58515e-05
        );

        //- Delayed neutron decay constants by group. The number of groups is
        //  deduced by the length of this list and should be consistend with
        //  Beta
        lambda
        (
            0.0125371
            0.0300828
            0.109879
            0.325484
            1.3036
            9.51817
        );

        //- If true, an additional logarithmic feedback of fuel temperature is
        //  added to the total reactivity
        fastNeutrons    false;

        //- Doppler coefficient for the logarithmic fuel temperature feedback.
        //  Used only if fastNeutrons is true
        feedbackCoeffFastDoppler 0;

        //- Reactivity contribution in the form
        //  feedbackCoeffTFuel*(TFuel-TFuelRef))
        feedbackCoeffTFuel -3e-06;

        //- Reactivity contribution in the form
        //  feedbackCoeffTClad*(TClad-TCladRef))
        feedbackCoeffTClad 0;

        //- Reactivity contribution in the form
        //  feedbackCoeffTCool*(TCool-TCoolRef))
        feedbackCoeffTCool 0;

        //- Reactivity contribution in the form
        //  feedbackCoeffRhoCool*(RhoCool-RhoCoolRef))
        feedbackCoeffRhoCool 0;

        //- Reactivity contribution in the form
        //  feedbackCoeffTStruct*(TStruct-TStructRef))
        feedbackCoeffTStruct 0;

        //- Used to compute drivelineExpansion as
        //  (TDriveline-TDrivelineRef)*absoluteDrivelineExpansionCoeff
        //  Dimensionally, it is m/K
        absoluteDrivelineExpansionCoeff 8e-05;

        //- Used to related drivelineExpansion to reactivity contribution
        //  Values of the expansion that are in between two map points
        //  result in a reactivity contribution that is computed via
        //  linear interpolation from the two. Any number of points can be
        //  provided. The points can be provided either in ascending or
        //  descending order of the first number (i.e. the expansion).
        controlRodReactivityMap
        (
            ( 0.1 -0.01 ) //- Pair in the form ( expans(m) reacitivity(-) )
            ( 0 0 )
            ( -0.1 0.01 )
        );

        //- Used to specify a reference fuel temperature for the computation
        //  of the TFuel and FastDoppler reactivity contributions. If not
        //  provided, the initial time-step flux-averaged fuel temperature is
        //  used, effectively starting from a steady state. The average is
        //  performed in the cellZones specified by fuelFeedbackZones
        //  TFuelRef    0;

        //- Used to specify a reference cladding temperature for the
        //  computation of the TClad reactivity contribution. If not
        //  provided, the initial time-step flux-averaged cladding temperature
        //  is used, effectively starting from a steady state. The average is
        //  performed in the cellZones specified by fuelFeedbackZones
        //  TCladRef    0;

        //- Used to specify a reference coolant temperature for the
        //  computation of the TCool reactivity contribution. If not
        //  provided, the initial time-step flux-averaged coolant temperature
        //  is used, effectively starting from a steady state. The average is
        //  performed in the cellZones specified by coolFeedbackZones
        //  TCoolRef    0;

        //- Used to specify a reference coolant density for the
        //  computation of the rhoCool reactivity contribution. If not
        //  provided, the initial time-step flux-averaged coolant density
        //  is used, effectively starting from a steady state. The average is
        //  performed in the cellZones specified by coolFeedbackZones
        //  rhoCoolRef    0;

        //- Used to specify a reference structure temperature for the
        //  computation of the TStruct reactivity contribution.
        //  If not provided, the initial time-step flux-averaged structure
        //  temperature is used, effectively starting from a steady state. The
        //  average is performed in the cellZones specified by
        //  structFeedbackZones
        //  TStructRef    0;

        //- Used to specify a reference driveline temperature for the
        //  computation of the driveline expansion reactivity contribution.
        //  If not provided, the initial time-step flux-averaged structure
        //  temperature is used, effectively starting from a steady state.The
        //  average is performed in the cellZones specified by
        //  drivelineFeedbackZones. Please note that the structure temperature
        //  is used as representative of the driveline in said zones
        //  TDrivelineRef    0;

        //- List of cellZones over which the flux-averaging of the fuel, clad
        //  temperature is performed to compute TFuel, TClad. If not provided,
        //  this defaults to the enitre mesh
        fuelFeedbackZones
        (
            "innerCore"
            "outerCore"
            ...
        );

        //- List of cellZones over which the flux-averaging of the coolant
        //  temperature and density are performed to compute TCool, rhoCool.
        //  If not provided, this defaults to the enitre mesh
        coolFeedbackZones
        (
            "innerCore"
            "outerCore"
            ...
        );

        //- List of cellZones over which the flux-averaging of the structure
        //  temperature is performed to compute TStruct. If not provided, this
        //  defaults to the enitre mesh
        structFeedbackZones
        (
            "diagrid"
            ...
        );

        //- List of cellZones over which the flux-averaging of the structure
        //  temperature is performed to compute TDriveline. If not provided,
        //  this defaults to the enitre mesh
        fuelFeedbackZones
        (
            "skirt"
            "driveline"
            ...
        );
    */
