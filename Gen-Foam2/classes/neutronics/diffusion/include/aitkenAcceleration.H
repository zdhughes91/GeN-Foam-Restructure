List<scalar> integralFlux0(xs_.energyGroups());
List<scalar> integralFlux1(xs_.energyGroups());
List<scalar> integralFlux2(xs_.energyGroups());

if(aitkenIterCounter==0) // only happening at timeStepIteration==1 
{
    forAll(fluxStar_, energyI)
    {
          integralFlux0[energyI]=(fvc::domainIntegrate(fluxStar_[energyI])).value();
    }   
}  
if(aitkenIterCounter==1)
{
    forAll(fluxStar_, energyI)
    {
          integralFlux1[energyI]=(fvc::domainIntegrate(fluxStar_[energyI])).value();
    }   
}
if(aitkenIterCounter==2)
{
    aitkenYesNo = true;
    forAll(fluxStar_, energyI)
    {
        integralFlux2[energyI]=(fvc::domainIntegrate(fluxStar_[energyI])).value();

        if(
        mag(integralFlux2[energyI]-integralFlux1[energyI])/(max(mag(integralFlux1[energyI]-integralFlux0[energyI]),SMALL)) >0.975 
        || (integralFlux2[energyI]-2*integralFlux1[energyI]+integralFlux0[energyI]==0)
        ) aitkenYesNo=false;
    }

    if(aitkenYesNo==true)
    {
        Info << "Intra-step Aitken acceleration no: " << aitkenIterNo << endl;
        forAll(fluxStar_, energyI)
        {
            integralFlux2[energyI]=(fvc::domainIntegrate(fluxStar_[energyI])).value();
            scalar aitkenScaling = 
            (
                ( 
                    integralFlux0[energyI]
                    - (integralFlux1[energyI]-integralFlux0[energyI])
                    * (integralFlux1[energyI]-integralFlux0[energyI])
                    / (integralFlux2[energyI]-2*integralFlux1[energyI]+integralFlux0[energyI])
                )/integralFlux2[energyI]
            );

            fluxStar_[energyI] = fluxStar_[energyI] * aitkenScaling ;
            
            flux_[energyI]=fluxStar_[energyI]/discFactor[energyI];

            integralFlux0[energyI]=(fvc::domainIntegrate(fluxStar_[energyI])).value();
        }
        #include "initializeNeutroSource.H"
        aitkenIterNo++ ;
        aitkenIterCounter = 0;
        Info<< "Done "  << endl;
    }
    else
    {
        Info << "Aitken not performed because of diverging series" << endl;

        forAll(integralFlux0, energyI)
        {
            integralFlux0[energyI] = integralFlux1[energyI];
            integralFlux1[energyI] = integralFlux2[energyI];
            
        }
        aitkenIterCounter-- ;
    }

}
aitkenIterCounter++ ;
