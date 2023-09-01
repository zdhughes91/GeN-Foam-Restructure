//- Get solution contro
const dictionary& neutronTransport =
mesh_.solutionDict().subDict("neutronTransport");

scalar neutronIterationResidual =
neutronTransport.lookupOrDefault<scalar>
(
    "neutronIterationResidual",
    0.00001
);
int maxNeutronIterations =
neutronTransport.lookupOrDefault<int>("maxNeutronIterations", 10);

//- Update reactivity coefficients
#include "computeFeedbackFieldValues.H"
#include "correctReactivity.H"


//- Update old time step values if time step has changed
if(timeIndex_ != mesh_.time().timeIndex())
{
    timeIndex_ = mesh_.time().timeIndex();
    powerOld_ = power_;
}

scalar dt(mesh_.time().deltaT().value());

//- prepare fields necessary to weight precursors
PtrList<scalar> weightedPrec(delayedGroups_);


//- Iterate precursors and power till sufficiently small prec residual
int neutroIter = 0;
scalar neutroResidual = 0.0;

//- Check if you need to initialize precurosr distribution
bool initialize = false;
if(initPrecursorsLiquidFuel_ &&
  (mesh_.time().timeIndex() == mesh_.time().startTimeIndex()+1) && neutroIter==0)
{
    initialize = true;
}

//Calc equilibrium precPK at the steady-state
if(initialize)
{
    forAll(precPK_,precI)
    {
        fvMatrix<scalar> precEq
        (
            fvm::Sp(dimensionedScalar("",dimless/dimTime,lambdas_[precI])*alphaPtr_(), precPKStar_[precI])
            - initOneGroupFluxN_*dimensionedScalar("",dimPower,power_)*betas_[precI]/
            dimensionedScalar("",dimTime,promptGenerationTime_)
             +   fvm::div(phiPtr_(), precPKStar_[precI], "div(phi_,precStar_)")
             -   fvm::laplacian(diffCoeffPrecPtr_(), precPKStar_[precI])
        );
        precEq.solve();
        precPKStar_[precI].storePrevIter();
    }
}

/*
//Calc the precEquilibriumReactivity
if(initialize)
{
    precEquilibriumReactivity_= 0.0;
    precEquilibriumReactivity_ += beta_;

    forAll(precPKStar_,precI)
    {
        //- Conver to precs per total volume
        precPK_[precI] = precPKStar_[precI] * alphaPtr_();
        precPK_[precI].correctBoundaryConditions();

        weightedPrec.set(precI,new scalar(fvc::domainIntegrate(
        initOneGroupFluxN_*precPK_[precI]).value()/domainIntegratedInitOneGroupFluxN_ ));

        precEquilibriumReactivity_ -= weightedPrec[precI] * lambdas_[precI] *
        promptGenerationTime_ / powerOld_;
    }
    #include "correctReactivity.H"
    Info << "precEquilibriumReactivity[pcm] =  " << precEquilibriumReactivity_ *1e+5 << nl << endl;
}
*/



//Calc beta for liquid fuel
liquidFuelBeta_ *= 0.0;
forAll(precPKStar_,precI)
{
    //- Convert to precs per total volume
    precPK_[precI] = precPKStar_[precI] * alphaPtr_();
    precPK_[precI].correctBoundaryConditions();

    weightedPrec.set(precI,new scalar(fvc::domainIntegrate(
    initOneGroupFluxN_*precPK_[precI]).value()/domainIntegratedInitOneGroupFluxN_ ));

    liquidFuelBeta_ += weightedPrec[precI] * lambdas_[precI] *
    promptGenerationTime_ / powerOld_;
}

//Calc the precEquilibriumReactivity
if(initialize)
{
    precEquilibriumReactivity_ = 0.0;
    precEquilibriumReactivity_ += beta_;
    precEquilibriumReactivity_ -= liquidFuelBeta_;
    #include "correctReactivity.H"
}

Info << "Beta[pcm] =  " << beta_ *1e+5 << nl << endl;
Info << "liquidFuelBeta[pcm] =  " << liquidFuelBeta_ *1e+5 << nl << endl;
Info << "precEquilibriumReactivity[pcm] =  " << precEquilibriumReactivity_ *1e+5 << nl << endl;


do
{
    Info << "Neutron iteration no:  " << neutroIter << nl << endl;

    //- Solve for precursors
    neutroResidual = 0.0;


    forAll(precPK_,precI)
    {
        //- Solve precursors
        fvMatrix<scalar> precEq
        (
            fvm::ddt(alphaPtr_(), precPKStar_[precI])
            +   fvm::Sp(dimensionedScalar("",dimless/dimTime,lambdas_[precI])*alphaPtr_(), precPKStar_[precI])
            -   initOneGroupFluxN_*dimensionedScalar("",dimPower,power_)*betas_[precI]/
            dimensionedScalar("",dimTime,promptGenerationTime_)
            +   fvm::div(phiPtr_(), precPKStar_[precI], "div(phi_,precStar_)")
            -   fvm::laplacian(diffCoeffPrecPtr_(), precPKStar_[precI])
        );

        scalar initRes = precEq.solve().max().initialResidual();
        neutroResidual = max(neutroResidual,initRes);

        //- Conver to precs per total volume
        precPK_[precI] = precPKStar_[precI] * alphaPtr_();
        precPK_[precI].correctBoundaryConditions();



        //- Calculate weighted precursors to be fed to power equation
        weightedPrec.set(precI,new scalar(fvc::domainIntegrate(
        initOneGroupFluxN_*precPK_[precI]).value()/domainIntegratedInitOneGroupFluxN_ ));

    }


    //- Calc power
    scalar sumPrecLamda = 0.0;
    forAll(precPK_,precI)
    {
        sumPrecLamda += weightedPrec[precI] * lambdas_[precI];
    }

    power_ = (powerOld_/dt+sumPrecLamda)/
    (
    1.0/dt - (totalReactivity_-beta_)/promptGenerationTime_
    );


    Info<<"Intermediate Power for the next iteration in the precursors equation[W] =  "<< power_ <<nl <<endl;
          
    neutroIter++;


    //- Some info on iteartions
    Info<< "neutroResidual (couplingIter: " << couplingIter
    << ", neutroIter: " << neutroIter << "): " << neutroResidual
    <<  endl;
}
while((neutroResidual > neutronIterationResidual)
      &&  (neutroIter < maxNeutronIterations)
    );

//- Re-scale fluxes, powerDensity
scalar fluxScaleFactor = power_/powerOld_;


for (int i = 0; i < energyGroups_; i++)
{
    fluxes_[i] *= fluxScaleFactor;
    fluxes_[i].correctBoundaryConditions();
}

oneGroupFlux_ *= fluxScaleFactor;
oneGroupFlux_.correctBoundaryConditions();

powerDensity_ *= fluxScaleFactor;
powerDensity_.correctBoundaryConditions();
secondaryPowerDenisty_ *= fluxScaleFactor;
secondaryPowerDenisty_.correctBoundaryConditions();

Info << endl << "pointKinetics: " << endl;
#include "pointKineticsInfo.H"
