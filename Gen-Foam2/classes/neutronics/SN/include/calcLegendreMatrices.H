
double legendrePoly(int, double);

Info<< "Calculating Legendre matrices\n" << endl;

Info<< "Warning: considering Legendre scattering moments up to 7^th at maximum." << endl;

scalarSquareMatrix muMatrix(discreteDirections_);

scalarSquareMatrix weightMatrix(discreteDirections_);

forAll(directionVersors_, directionI)
{
    forAll(directionVersors_, directionJ)
    {
        muMatrix[directionI][directionJ] = directionVersors_[directionI] & directionVersors_[directionJ];
        weightMatrix[directionI][directionJ] = directionWeights_[directionJ];
    }
}

Info<< "Done\n" << endl;

PtrList<scalar> maxLegendreIntegral(xs_.legendreMoments());
PtrList<scalar> minLegendreIntegral(xs_.legendreMoments());
PtrList<scalar> maxAbsLegendreIntegral(xs_.legendreMoments());
PtrList<scalar> minAbsLegendreIntegral(xs_.legendreMoments());
scalar legendreIntegral;
scalar legendreAbsIntegral;

forAll(legendreMatrices_, legendreI)
{

    Info<< "Calculating Legendre matrices... Legendre moment: " << Foam::name(legendreI) << endl;
    maxLegendreIntegral.set(legendreI, new scalar (-100.0));
    minLegendreIntegral.set(legendreI, new scalar (100.0));
    maxAbsLegendreIntegral.set(legendreI, new scalar (-100.0));
    minAbsLegendreIntegral.set(legendreI, new scalar (100.0));
    Info<< "Also, checking consistency of the matrix" << endl;
    legendreMatrices_.set
    (
        legendreI,
        new scalarSquareMatrix(discreteDirections_)
    );
    forAll(directionVersors_, directionI)
    {
        legendreIntegral = 0;
        legendreAbsIntegral = 0;
        forAll(directionVersors_, directionJ)
        {
            legendreMatrices_[legendreI][directionI][directionJ] = legendrePoly(legendreI, muMatrix[directionI][directionJ]);
            legendreIntegral += legendreMatrices_[legendreI][directionI][directionJ] * directionWeights_[directionJ] * Foam::constant::mathematical::pi * 4;
            legendreAbsIntegral += fabs(legendreMatrices_[legendreI][directionI][directionJ]) * directionWeights_[directionJ] * Foam::constant::mathematical::pi * 4;
        }
        maxLegendreIntegral[legendreI] = max(maxLegendreIntegral[legendreI],legendreIntegral);
        minLegendreIntegral[legendreI] = min(minLegendreIntegral[legendreI],legendreIntegral);
        maxAbsLegendreIntegral[legendreI] = max(maxAbsLegendreIntegral[legendreI],legendreAbsIntegral);
        minAbsLegendreIntegral[legendreI] = min(minAbsLegendreIntegral[legendreI],legendreAbsIntegral);
    }
}

Info<< "Done" << nl << endl;

Info<< "Consistency report:" << nl << "max" << nl << "min" << nl << "maxAbs" << nl << "minAbs" << nl << endl;

forAll(legendreMatrices_, legendreI)
{
    Info<< "P" << Foam::name(legendreI) << " matrix:" << endl;
    Info<< Foam::name(minLegendreIntegral[legendreI]) << endl;
    Info<< Foam::name(maxLegendreIntegral[legendreI]) << endl;
    Info<< Foam::name(minAbsLegendreIntegral[legendreI]) << endl;
    Info<< Foam::name(maxAbsLegendreIntegral[legendreI]) << nl << endl;
}

Info<< "Abs integrals should be:  (0 to 7)" << nl << "12.5664" << nl << "6.28319" << nl << "4.8366" << nl << "4.08407" << nl << "3.6014" << nl << "3.25777" << nl << "2.997" << nl << "2.79031" << nl << endl;    


Info<< "Done\n" << endl;







