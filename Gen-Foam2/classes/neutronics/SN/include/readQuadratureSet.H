Info<< "Reading quadrature set\n" << endl;

label discreteDirectionsOct = discreteDirections_/8;

scalarField directionUxOct("directionUxOct", quadratureSet_, discreteDirectionsOct);
scalarField directionUyOct("directionUyOct", quadratureSet_, discreteDirectionsOct);
scalarField directionUzOct("directionUzOct", quadratureSet_, discreteDirectionsOct);
scalarField directionWeightsOct("directionWeightsOct", quadratureSet_, discreteDirectionsOct);

vectorField octants(8);
octants[0] = vector(1.00, 1.00, 1.00);
octants[1] = vector(-1.0, 1.00, 1.00);
octants[2] = vector(-1.0, -1.0, 1.00);
octants[3] = vector(1.00, -1.0, 1.00);
octants[4] = vector(1.00, 1.00, -1.0);
octants[5] = vector(-1.0, 1.00, -1.0);
octants[6] = vector(-1.0, -1.0, -1.0);
octants[7] = vector(1.00, -1.0, -1.0);

forAll(octants, octI)
{
    forAll(directionWeightsOct, directionI)
    {
        directionVersors_[directionI+octI*discreteDirectionsOct].x() = directionUxOct[directionI] * octants[octI].x();
        directionVersors_[directionI+octI*discreteDirectionsOct].y() = directionUyOct[directionI] * octants[octI].y();
        directionVersors_[directionI+octI*discreteDirectionsOct].z() = directionUzOct[directionI] * octants[octI].z();
        directionWeights_[directionI+octI*discreteDirectionsOct] = directionWeightsOct[directionI];
    }
}
    
Info<< "Done\n" << endl;

Info<< "Checking quadrature consistency\n" << endl;

scalar totalDirectionsWeight = 0.0;
scalar maxVersorLenght = 0.0;
scalar minVersorLenght = 2.0;

forAll(directionVersors_, directionI)
{
    totalDirectionsWeight += directionWeights_[directionI];
    maxVersorLenght = max(mag(directionVersors_[directionI]),maxVersorLenght);
    minVersorLenght = min(mag(directionVersors_[directionI]),minVersorLenght);
}
    

Info<< "Normalizing to 1 instead of 4*pi..." << endl;

directionWeights_ *= 1/totalDirectionsWeight;

forAll(directionVersors_, directionI)
{
    directionVersors_[directionI] *= 1/mag(directionVersors_[directionI]);
}
    
Info<< "Done" << endl;

Info<< "Checking quadrature consistency... again " << endl;

totalDirectionsWeight = 0.0;
maxVersorLenght = 0.0;
minVersorLenght = 2.0;

forAll(directionVersors_, directionI)
{
    totalDirectionsWeight += directionWeights_[directionI];
    maxVersorLenght = max(mag(directionVersors_[directionI]),maxVersorLenght);
    minVersorLenght = min(mag(directionVersors_[directionI]),minVersorLenght);
}
    
Info<< "Done" << endl;


Info<< "Check report:" << endl;
Info<< "totalDirectionsWeight: " << totalDirectionsWeight << endl;

Info<< "maxVersorLenght: " << maxVersorLenght <<  endl;

Info<< "minVersorLenght: " << minVersorLenght << nl << endl;


