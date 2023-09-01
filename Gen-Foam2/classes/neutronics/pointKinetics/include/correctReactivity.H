scalar DopplerReactivity =
    -coeffFastDoppler_*
    Foam::log(max(TFuelValue, SMALL)/max(TFuelRef_, SMALL));

scalar TFuelReactivity =
    coeffTFuel_*(TFuelValue - TFuelRef_);

scalar TCladReactivity =
    coeffTClad_*(TCladValue - TCladRef_);

scalar TCoolReactivity =
    coeffTCool_*(TCoolValue - TCoolRef_);

scalar rhoCoolReactivity =
    coeffRhoCool_*(rhoCoolValue - rhoCoolRef_);

scalar TStructReactivity =
    coeffTStruct_*(TStructValue - TStructRef_);

scalar drivelineExpValue
(
    coeffDrivelineExp_*(TDrivelineValue - TDrivelineRef_)
);
scalar drivelineReactivity(calcDrivelineReactivity(drivelineExpValue));

Pair<scalar> GEMLevelAndReactivity(calcGEMLevelAndReactivity());
scalar GEMSodiumLevel(GEMLevelAndReactivity.first());
scalar GEMReactivity(GEMLevelAndReactivity.second());

(void) GEMSodiumLevel;

//- External reactivity
scalar t(mesh_.time().timeOutputValue());
scalar externalReactivity = externalReactivityTimeProfile_.value(t);

//- Boron reactivity
scalar boronReactivity = boronReactivityTimeProfile_.value(t);
// scalar boronReactivity = 0.0;
// if(boronReactivityPtr_.valid())
// {
//     boronReactivity =
//         boronReactivityPtr_->value
//         (
//             mesh_.time().timeOutputValue()
//         );
// }

totalReactivity_ =
        precEquilibriumReactivity_
    +   externalReactivity
    +   DopplerReactivity
    +   TFuelReactivity
    +   TCladReactivity
    +   TCoolReactivity
    +   rhoCoolReactivity
    +   TStructReactivity
    +   drivelineReactivity
    +   GEMReactivity
    +   boronReactivity;
