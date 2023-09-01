//- Field names must reflect those defined in createCouplingFields.H
TFuelOrig_ = 
    src.findObject<volScalarField>("bafflelessTFuelAv");
TCladOrig_ = 
    src.findObject<volScalarField>("bafflelessTCladAv");
TCoolOrig_ = 
    src.findObject<volScalarField>("bafflelessTCool");
rhoCoolOrig_ = 
    src.findObject<volScalarField>("bafflelessRhoCool");
if (liquidFuel_)
{
    UOrig_ = 
        src.findObject<volVectorField>("bafflelessU");
    alphaOrig_ = 
        src.findObject<volScalarField>("bafflelessAlpha");
    alphatOrig_ = 
        src.findObject<volScalarField>("bafflelessAlphat");
    muOrig_ =
        src.findObject<volScalarField>("bafflelessMu");
}
else
{
    UOrig_ = nullptr;
    alphaOrig_ = nullptr;
    alphatOrig_ = nullptr;
    muOrig_ = nullptr;
}

//- Initialize mapped fields
this->interpolateCouplingFields(neutroToFluid);

//- Initialize XS correctly
xs_.correct(TFuel_, TClad_, rhoCool_, TCool_, disp_);
