/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Description
    Transient segregated finite-volume solver of linear-elastic,
    small-strain deformation of a solid body, with optional thermal
    diffusion and thermal stresses.

    Simple linear elasticity structural analysis code.
    Solves for the displacement vector field D, also generating the
    stress tensor field sigma.

    Derived from soliddisp_lacementFoam, modified for mu_lti-zone
    applications through use of a dictionary specifying the properties in each
    mesh region

Reference publication (for the modification of soliddisp_lacementFoam)
    Carlo Fiorina, Ivor Clifford, Manuele Aufiero, Konstantin Mikityuk,
    "GeN-Foam: a novel OpenFOAM® based mu_lti-physics solver for 2D/3D transient
    analysis of nuclear reactors", Nuclear Engineering and Design, 2015

Author (for the modification of soliddisp_lacementFoam)
    Carlo Fiorina <carlo.fiorina@outlook.com; carlo.fiorina@epfl.ch;>

\*---------------------------------------------------------------------------*/

int iCorr = 0;
scalar initialResidual = 0;

do
{
    Info << "thermoMechanical iteration no:  " << iCorr << endl;

    {
        fvVectorMatrix DEqn
        (
            fvm::d2dt2(disp_)
         ==
            fvm::laplacian(2*mu_ + lambda_, disp_, "laplacian(DD,D)")
          + divSigmaExp_

        );

        {
            DEqn += fvc::grad(threeKalpha_*(TStruct_-TStructRef_));
        }
        DEqn.relax();
        initialResidual = DEqn.solve().max().initialResidual();

        if (!compactNormalStress_)
        {
            divSigmaExp_ = fvc::div(DEqn.flux());
        }
    }

    {
        volTensorField gradD(fvc::grad(disp_));
        sigmaD_ = mu_*twoSymm(gradD) + (lambda_*I)*tr(gradD);

        if (compactNormalStress_)
        {
            divSigmaExp_ = fvc::div
            (
                sigmaD_ - (2*mu_ + lambda_)*gradD,
                "div(sigmaD)"
            );
        }
        else
        {
            divSigmaExp_ += fvc::div(sigmaD_);
        }
    }

} while (initialResidual > convergenceTolerance_ && ++iCorr < nCorr_);

residual = max(residual,initialResidual);

#include "calculateStress.H"

//- Calculating axial fuel expansion
{

    surfaceScalarField faceArea(fuelOrientation_ & mesh_.Sf());

    fvScalarMatrix fuelDispEqn
    (
        fvm::div(faceArea,fuelDisp_, "div(fuelDisp)")
     ==
        alphaFuel_*(TFuel_-TFuelRef_)
    );

    vector CROrientation = -fuelOrientation_;
    surfaceScalarField faceAreaCR(CROrientation & mesh_.Sf());

    fvScalarMatrix CRDispEqn
    (
        fvm::div(faceAreaCR,CRDisp_, "div(CRDisp)")
     ==
        -alphaCR_*(TStruct_-TCRRef_)
    );

    fuelDispEqn.solve();
    CRDispEqn.solve();

    fuelDispVector_ = (fuelDisp_+CRDisp_)*fuelOrientation_;

    meshDisp_ = 
        (disp_-((disp_ & fuelOrientation_) * fuelOrientation_)) 
    +   fuelDispVector_ ; //- Remove any axial variation

    meshDisp_.correctBoundaryConditions();

}
// ************************************************************************* //
