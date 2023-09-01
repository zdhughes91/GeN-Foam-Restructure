/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

Global
    myMULES

Description
    myMULES: Multidimensional universal limiter for explicit solution.

    Solve a convective-only transport equation using an explicit universal
    multi-dimensional limiter.

    Parameters are the variable to solve, the normal convective flux and the
    actual explicit flux of the variable which is also used to return limited
    flux used in the bounded-solution. 
    This version is virtually indentical to the MULES present in OpenFOAM6.
    The main difference consists in the option to provide a space-dependent
    limiter field, something which is absent in the v1812 OpenFOAM version of
    MULES. 
    Modified by Stefan Radman (stefan.radman@epfl.ch/stefanradman92@gmail.com)

SourceFiles
    myMULES.C
    myMULESTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef myMULES_H
#define myMULES_H

#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "primitiveFieldsFwd.H"
#include "geometricOneField.H"
#include "zero.H"
#include "zeroField.H"
#include "UPtrList.H"
#include "HashSet.H"
#include "UniformField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace myMULES
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class RdeltaTType, class RhoType, class SpType, class SuType>
void explicitSolve
(
    const RdeltaTType& rDeltaT,
    const RhoType& rho,
    volScalarField& psi,
    const surfaceScalarField& phiPsi,
    const SpType& Sp,
    const SuType& Su
);

template<class RhoType>
void explicitSolve
(
    const RhoType& rho,
    volScalarField& psi,
    const surfaceScalarField& phiPsi
);

template<class RhoType, class SpType, class SuType>
void explicitSolve
(
    const RhoType& rho,
    volScalarField& psi,
    const surfaceScalarField& phiPsi,
    const SpType& Sp,
    const SuType& Su
);

template<class RhoType, class PsiMaxType, class PsiMinType>
void explicitSolve
(
    const RhoType& rho,
    volScalarField& psi,
    const surfaceScalarField& phiBD,
    surfaceScalarField& phiPsi,
    const PsiMaxType& psiMax,
    const PsiMinType& psiMin
);

template
<
    class RhoType,
    class SpType,
    class SuType,
    class PsiMaxType,
    class PsiMinType
>
void explicitSolve
(
    const RhoType& rho,
    volScalarField& psi,
    const surfaceScalarField& phiBD,
    surfaceScalarField& phiPsi,
    const SpType& Sp,
    const SuType& Su,
    const PsiMaxType& psiMax,
    const PsiMinType& psiMin
);

template
<
    class RdeltaTType,
    class RhoType,
    class SpType,
    class SuType,
    class PsiMaxType,
    class PsiMinType
>
void explicitSolve
(
    const RdeltaTType& rDeltaT,
    const RhoType& rho,
    volScalarField& psi,
    const surfaceScalarField& phiBD,
    surfaceScalarField& phiPsi,
    const SpType& Sp,
    const SuType& Su,
    const PsiMaxType& psiMax,
    const PsiMinType& psiMin
);

template
<
    class RdeltaTType,
    class RhoType,
    class SpType,
    class SuType,
    class PsiMaxType,
    class PsiMinType
>
void limiter
(
    scalarField& allLambda,
    const RdeltaTType& rDeltaT,
    const RhoType& rho,
    const volScalarField& psi,
    const surfaceScalarField& phiBD,
    const surfaceScalarField& phiCorr,
    const SpType& Sp,
    const SuType& Su,
    const PsiMaxType& psiMax,
    const PsiMinType& psiMin
);

template
<
    class RdeltaTType,
    class RhoType,
    class SpType,
    class SuType,
    class PsiMaxType,
    class PsiMinType
>
void limit
(
    const RdeltaTType& rDeltaT,
    const RhoType& rho,
    const volScalarField& psi,
    const surfaceScalarField& phi,
    surfaceScalarField& phiPsi,
    const SpType& Sp,
    const SuType& Su,
    const PsiMaxType& psiMax,
    const PsiMinType& psiMin,
    const bool returnCorr
);

template
<
    class RhoType,
    class SpType,
    class SuType,
    class PsiMaxType,
    class PsiMinType
>
void limit
(
    const RhoType& rho,
    const volScalarField& psi,
    const surfaceScalarField& phi,
    surfaceScalarField& phiPsi,
    const SpType& Sp,
    const SuType& Su,
    const PsiMaxType& psiMax,
    const PsiMinType& psiMin,
    const bool returnCorr
);

void limitSum(UPtrList<scalarField>& phiPsiCorrs);

template<class SurfaceScalarFieldList>
void limitSum(SurfaceScalarFieldList& phiPsiCorrs);

void limitSum
(
    const UPtrList<const scalarField>& alphas,
    UPtrList<scalarField>& phiPsiCorrs,
    const labelHashSet& fixed
);

template<class SurfaceScalarFieldList>
void limitSum
(
    const SurfaceScalarFieldList& alphas,
    SurfaceScalarFieldList& phiPsiCorrs,
    const labelHashSet& fixed
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace myMULES
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "myMULESTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
