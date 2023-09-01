Info<< "    totalPower = " << power_ << " W" << endl;
Info<< "    -> fission = " << fissionPower_ << " W" << endl;
Info<< "    -> decay   = " << decayPower_ << " W" << endl;
Info<< "    totalReactivity   = " << (totalReactivity_*1e5)<< " pcm" << endl;
Info<< "    -> extReactivity  = " << (externalReactivity*1e5)<< " pcm" << endl;
Info<< "    -> Doppler (fast) = " << (DopplerReactivity*1e5)<< " pcm" << endl;
Info<< "    -> TFuel          = " << (TFuelReactivity*1e5) << " pcm" << endl;
Info<< "    -> TClad          = " << (TCladReactivity*1e5) << " pcm" << endl;
Info<< "    -> TCool          = " << (TCoolReactivity*1e5) << " pcm" << endl;
Info<< "    -> rhoCool        = " << (rhoCoolReactivity*1e5) << " pcm" << endl;
Info<< "    -> TStruct        = " << (TStructReactivity*1e5) << " pcm" << endl;
Info<< "    -> driveline      = " << (drivelineReactivity*1e5) << " pcm"<<endl;
Info<< "    -> Boron      = " << (boronReactivity*1e5) << " pcm"<<endl;
if(liquidFuel_)
{
    Info<< "    Initial reactivity loss due to delayed neutron precursors movement = "
    << (precEquilibriumReactivity_ *1e+5) << " pcm"<<endl;
    Info<< "    Estimated Beta including DNP recirculation = "
    << (liquidFuelBeta_ *1e+5) << " pcm"<<endl;
}
if (this->get<bool>("GEM"))
	Info<< "    -> GEM            = " << (GEMReactivity*1e5) << " pcm" << endl;
Info<< "    TFuel = " << TFuelValue << " K" << endl;
Info<< "    TClad = " << TCladValue << " K" << endl;
Info<< "    TCool = " << TCoolValue << " K" << endl;
Info<< "    rhoCool = " << rhoCoolValue << " kg/m3" << endl;
Info<< "    TStruct = " << TStructValue << " K" << endl;
Info<< "    driveline (T, exp) = "
    << TDrivelineValue << " K, "
    << drivelineExpValue << " m"
    << endl;
if (this->get<bool>("GEM"))
	Info<< "    GEMSodiumLevel = " << GEMSodiumLevel << " m" << endl;
