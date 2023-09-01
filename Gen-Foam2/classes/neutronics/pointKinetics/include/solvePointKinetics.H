//- Update reactivity coefficients
#include "computeFeedbackFieldValues.H"
#include "correctReactivity.H"

//- Update decayPower_, if available
if (decayPowerPtr_.valid())
{
    decayPower_ = decayPowerPtr_->value
    (
        mesh_.time().timeOutputValue() - decayPowerStartTime_
    );
}

//- Update fissionPower_ and precursorPowers_ via actual pointKinetics model

//- Update old time step values if time step has changed
if (timeIndex_ != mesh_.time().timeIndex())
{
    timeIndex_ = mesh_.time().timeIndex();
    powerOld_ = power_;
    fissionPowerOld_ = fissionPower_;
    for (int i=0; i<delayedGroups_; i++)
    {
        precursorPowersOld_[i] = precursorPowers_[i];
    }
}
label n(delayedGroups_+1);
scalar dt(mesh_.time().deltaT().value());

//- Coefficient matrix, source, solution field
SquareMatrix<scalar> A(n, 0.0);
List<scalar> B(n, 0.0);
List<scalar> x(n, 0.0);

//- First line of matrix has the coefficients of the discretized power
//  differential equation, the first element of the solution vector is the
//  power
A[0][0] = 1.0/dt-(totalReactivity_-beta_)/promptGenerationTime_;
for (int i = 1; i < n; i++)
{
    A[0][i] = -lambdas_[i-1];
    B[0] = fissionPowerOld_/dt;
}

//- Rest of the matrix for the precursor coefficients
for (int i = 1; i < n; i++)
{
    A[i][0] = -betas_[i-1]/promptGenerationTime_;
    A[i][i] = 1.0/dt + lambdas_[i-1];
    B[i] = precursorPowersOld_[i-1]/dt;
}

//- Solve
solve(x, A, B);
fissionPower_ = x[0];
for (int i=1; i<n; i++)
{
    precursorPowers_[i-1] = x[i];
}

//- Update total power
power_ = fissionPower_ + decayPower_;

//- Re-scale fluxes, precursors by change in fissionPower
scalar fluxScaleFactor = fissionPower_/fissionPowerOld_;

for (int i = 0; i < energyGroups_; i++)
{
    fluxes_[i] *= fluxScaleFactor;
    fluxes_[i].correctBoundaryConditions();
}
oneGroupFlux_ *= fluxScaleFactor;
oneGroupFlux_.correctBoundaryConditions();

//- Re-scale powerDensity by change in total power
scalar powerDensityScaleFactor = power_/powerOld_;
powerDensity_ *= powerDensityScaleFactor;
powerDensity_.correctBoundaryConditions();
secondaryPowerDenisty_ *= powerDensityScaleFactor;
secondaryPowerDenisty_.correctBoundaryConditions();

//- Only possible to re-scale precursors if they were initially present,
//  so I'm using a forAllIter rather than a for, as delayedGroups_ is not
//  necessarily representative of the length of the precursors_ PtrList
int i = 0;
forAllIter
(
    PtrList<volScalarField>,
    precursors_,
    iter
)
{
    if (precursorPowersOld_[i] != 0)
    {
        volScalarField& precursori(iter());
        precursori *= precursorPowers_[i]/precursorPowersOld_[i];
    }
    i++;
}

//- Write precursorPowers to reactorState
if (TFuel_.mesh().time().write())
{
    reactorState_.set("pTarget", power_);
    reactorState_.set("precursorPowers", precursorPowers_);
}

Info << endl << "pointKinetics: " << endl;
#include "pointKineticsInfo.H"
