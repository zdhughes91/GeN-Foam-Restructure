Info << "Solving neutronics" << endl ;

xs_.correct(TFuel_, TClad_, rhoCool_, TCool_, disp_);
#include "setXSFieldsSN.H"

const dictionary& neutronTransport = mesh_.solutionDict().subDict("neutronTransport");
//const bool integralPredictor = neutronTransport.lookupOrDefault("integralPredictor", false);
//const bool implicitPredictor = neutronTransport.lookupOrDefault("implicitPredictor", false);
//const bool aitkenAcceleration = neutronTransport.lookupOrDefault("aitkenAcceleration", false);
const bool aitkenAccelerationEigenvalue =
    neutronTransport.lookupOrDefault("aitkenAccelerationEigenvalue", false);

scalar neutronIterationResidual =
    neutronTransport.lookupOrDefault<scalar>("neutronIterationResidual", 0.00001);

int maxNeutronIterations =
    neutronTransport.lookupOrDefault<int>("maxNeutronIterations", 10);

scalar innToOutTolerance =
    neutronTransport.lookupOrDefault<scalar>("innToOutTolerance", 0.01);

int maxInnerIter =
    neutronTransport.lookupOrDefault<int>("maxInnerIter", 3);


//for aitken
//List<scalar> integralFlux0(xs_.energyGroups());
//List<scalar> integralFlux1(xs_.energyGroups());
//List<scalar> integralFlux2(xs_.energyGroups());

if(aitkenAccelerationEigenvalue)
{
    #include "aitkenAccelerationEigenvalue.H"
}


#include "updateScalarFluxes.H"
scalar pTot(power());
#include "initializeNeutroSource.H"
#include "precEq.H"
#include "initializeDelayedNeutroSource.H"
if (eigenvalueNeutronics_ && (mesh_.time().timeIndex()-mesh_.time().startTimeIndex())==1 )
{
    #include "normFluxesSN.H"
}
#include "updateScalarFluxes.H"

/*
if(!eigenvalueNeutronics_ && couplingIter == 0)
{

    if(integralPredictor)
    {
        #include "integralPredictor.H"
    }
    if(implicitPredictor)
    {
        #include "implicitPredictor.H"
    }

}
*/

int neutroIter = 0;

//#include "createFluxMatricesSN.H"// Can be used to speed up, but costs a lot in memory

scalar neutroResidual = 0.0;

do
{
    Info << "Neutron iteration no:  " << neutroIter << nl << endl;

    #include "fluxEqSN.H"
    //The part here below should be included for transients
    //#include "updateScalarFluxes.H"
    //#include "initializeNeutroSource.H"
    //#include "precEq.H"
    //#include "initializeDelayedNeutroSource.H"

    neutroIter++;

    Info<< "neutroResidual (couplingIter: " << couplingIter << ",neutroIter: "
    << neutroIter << "): " << neutroResidual << endl << " after Aitken no: " << aitkenIterNo_ << endl;

}while( (neutroResidual > neutronIterationResidual) && (neutroIter < maxNeutronIterations) && !eigenvalueNeutronics_);

#include "updateScalarFluxes.H"

powerDensity_*=0.0;
forAll(flux_,energyI)
{
    powerDensity_ += flux_[energyI] * sigmaPow[energyI];
}

pTot = power();

Info<< "Final neutroResidual at coupling/neutroIter iter no " << couplingIter << "/" << neutroIter << ": "  << neutroResidual << " - pow:" << pTot << endl;

//- Determine power densities
secondaryPowerDenisty_ = powerDensity_  / max(secondaryPowerVolumeFraction, SMALL) * fractionToSecondaryPower  ;
powerDensity_ /= max(fuelFraction, SMALL) * (1.0 - fractionToSecondaryPower);



if (eigenvalueNeutronics_)
{
    Info << "pTarget_ = " << pTarget_ << nl << endl;
    Info << "keff_ previous = " << keff_ << nl << endl;
    keff_ *= pTot/pTarget_;
    keff_ = min(keff_,5);
    keff_ = max(keff_,0.1);
    Info << "keff_ = " << keff_ << nl << endl;

    // Since a keff_ is a scalar, the value inside the dictionary is not updated at runTime automatically.
    // The line below resets the SCALAR IN THE DICTIONARY at runTime so that the dictionary is written with the updated value.
    reactorState_.set("keff", keff_);
    // reactorState_.regIOobject::writeObject
    // (
    //     IOstream::ASCII,
    //     IOstream::currentVersion,
    //     reactorState_.time().writeCompression(),
    //     true
    // );

    #include "normFluxesSN.H"

}



residual = max(residual,neutroResidual);
