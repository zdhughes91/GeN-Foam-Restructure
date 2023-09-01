xs_.correct(TFuel_, TClad_, rhoCool_, TCool_, disp_);
#include "setXSFields.H"

const dictionary& neutronTransport = mesh_.solutionDict().subDict("neutronTransport");

const bool integralPredictor = neutronTransport.lookupOrDefault("integralPredictor", false);
const bool implicitPredictor = neutronTransport.lookupOrDefault("implicitPredictor", false);
const bool aitkenAcceleration = neutronTransport.lookupOrDefault("aitkenAcceleration", false);
scalar neutronIterationResidual =
    neutronTransport.lookupOrDefault<scalar>("neutronIterationResidual", 0.00001);
int maxNeutronIterations =
    neutronTransport.lookupOrDefault<int>("maxNeutronIterations", 10);

//for aitken
List<scalar> integralFlux0(xs_.energyGroups());
List<scalar> integralFlux1(xs_.energyGroups());
List<scalar> integralFlux2(xs_.energyGroups());

scalar pTot(1.0);

forAll(flux_,energyI)
{
    flux_[energyI]=(fluxStar_[energyI]-fluxStar2_[energyI]*2.0)/discFactor[energyI];
}

if (eigenvalueNeutronics_ && (mesh_.time().timeIndex()-mesh_.time().startTimeIndex())==1 )
{
    pTot = power();
    #include "normFluxesSP3.H"
    pTot = power();
}

#include "initializeNeutroSource.H"
#include "precEq.H"

#include "initializeDelayedNeutroSource.H"

if(!eigenvalueNeutronics_ && couplingIter == 0)
{

    if(integralPredictor)
    {
        #include "integralPredictorSP3.H"
    }
    if(implicitPredictor)
    {
        #include "implicitPredictorSP3.H"
    }

}

int aitkenIterCounter = 0;
int aitkenIterNo = 0;
bool aitkenYesNo = true;

int neutroIter = 0;

#include "createFluxMatricesSP3.H"

scalar neutroResidual = 0.0;

do
{
    Info << "Neutron iteration no:  " << neutroIter << nl << endl;

    #include "fluxEqSP3.H"
    #include "precEq.H"
    #include "initializeDelayedNeutroSource.H"

    neutroIter++;

    Info<< "neutroResidual (couplingIter: " << couplingIter << ",neutroIter: "
    << neutroIter << "): " << neutroResidual << " after Aitken no: " << aitkenIterNo << endl;

}while( (neutroResidual > neutronIterationResidual) && (neutroIter < maxNeutronIterations) && !eigenvalueNeutronics_);

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

    #include "normFluxesSP3.H"
    xs_.adjustDiscFactors(fluxStar_);
}

residual = max(residual,neutroResidual);

oneGroupFlux_ *= 0.0;
forAll(flux_,energyI)
{
    oneGroupFlux_+=flux_[energyI];
}
