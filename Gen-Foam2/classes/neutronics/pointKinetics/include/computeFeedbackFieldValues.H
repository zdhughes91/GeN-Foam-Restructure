//- Create field for weighting. This is currently done on the oneGroupFlux,
//  rather than on the energy-dependent adjoint. It will be changed in the
//  future
volScalarField oneGroupFluxSqr(sqr(oneGroupFlux_));

scalar domainIntegratedFluxSqrFuel
(
    fvc::domainIntegrate(oneGroupFluxSqr*fuelFeedbackCellField_).value()
);
scalar domainIntegratedFluxSqrCool
(
    fvc::domainIntegrate(oneGroupFluxSqr*coolFeedbackCellField_).value()
);

//- Compute average perturbed values
scalar TFuelValue
(
    fvc::domainIntegrate
    (
        oneGroupFluxSqr*TFuel_*fuelFeedbackCellField_
    ).value()/
    domainIntegratedFluxSqrFuel
);
scalar TCladValue
(
    fvc::domainIntegrate
    (
        oneGroupFluxSqr*TClad_*fuelFeedbackCellField_
    ).value()/
    domainIntegratedFluxSqrFuel
);
scalar TCoolValue
(
    fvc::domainIntegrate
    (
        oneGroupFluxSqr*TCool_*coolFeedbackCellField_
    ).value()/
    domainIntegratedFluxSqrCool
);
scalar rhoCoolValue
(
    fvc::domainIntegrate
    (
        oneGroupFluxSqr*rhoCool_*coolFeedbackCellField_
    ).value()/
    domainIntegratedFluxSqrCool
);
scalar TStructValue
(
    fvc::domainIntegrate
    (
        oneGroupFluxSqr*TStruct_*structFeedbackCellField_
    ).value()/
    fvc::domainIntegrate(oneGroupFluxSqr*structFeedbackCellField_).value()
);
scalar TDrivelineValue
(
    fvc::domainIntegrate
    (
        oneGroupFluxSqr*TStruct_*drivelineFeedbackCellField_
    ).value()/
    fvc::domainIntegrate(oneGroupFluxSqr*drivelineFeedbackCellField_).value()
);
