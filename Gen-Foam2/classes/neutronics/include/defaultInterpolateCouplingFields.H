neutroToFluid.mapTgtToSrc(*TFuelOrig_, plusEqOp<scalar>(), TFuel_);
neutroToFluid.mapTgtToSrc(*TCladOrig_, plusEqOp<scalar>(), TClad_);
neutroToFluid.mapTgtToSrc(*TCoolOrig_, plusEqOp<scalar>(), TCool_);
neutroToFluid.mapTgtToSrc(*rhoCoolOrig_, plusEqOp<scalar>(), rhoCool_);
if (liquidFuel_)
{
    neutroToFluid.mapTgtToSrc(*UOrig_, plusEqOp<vector>(), UPtr_());
    neutroToFluid.mapTgtToSrc
    (
        *alphaOrig_, 
        plusEqOp<scalar>(), 
        alphaPtr_()
    );
    neutroToFluid.mapTgtToSrc
    (
        *alphatOrig_, 
        plusEqOp<scalar>(), 
        alphatPtr_()
    );
    neutroToFluid.mapTgtToSrc(*muOrig_, plusEqOp<scalar>(), muPtr_());
    phiPtr_() = fvc::flux(UPtr_());
    volScalarField diffCoeffOrig
    (
        (
            *alphatOrig_ 
        +   *muOrig_/xs_.ScNo()
        )/(*rhoCoolOrig_)
    ); 
    neutroToFluid.mapTgtToSrc
    (
        diffCoeffOrig, 
        plusEqOp<scalar>(), 
        diffCoeffPrecPtr_()
    );

    UPtr_().correctBoundaryConditions();
    alphaPtr_().correctBoundaryConditions();
    alphatPtr_().correctBoundaryConditions();
    diffCoeffPrecPtr_().correctBoundaryConditions();
}

TFuel_.correctBoundaryConditions();
TClad_.correctBoundaryConditions();
TCool_.correctBoundaryConditions();
rhoCool_.correctBoundaryConditions();
