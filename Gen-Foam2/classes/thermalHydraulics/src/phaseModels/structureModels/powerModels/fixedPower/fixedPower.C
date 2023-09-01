/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
|                                                                             |
|    Built on OpenFOAM v2212                                                  |
|    Copyright 2011-2016 OpenFOAM Foundation, 2017-2022 OpenCFD Ltd.         |
-------------------------------------------------------------------------------
License
    This file is part of GeN-Foam.

    GeN-Foam is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    GeN-Foam is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    This offering is not approved or endorsed by the OpenFOAM Foundation nor
    OpenCFD Limited, producer and distributor of the OpenFOAM(R)software via
    www.openfoam.com, and owner of the OPENFOAM(R) and OpenCFD(R) trademarks.

    This particular snippet of code is developed according to the developer's
    knowledge and experience in OpenFOAM. The users should be aware that
    there is a chance of bugs in the code, though we've thoroughly test it.
    The source code may not be in the OpenFOAM coding style, and it might not
    be making use of inheritance of classes to full extent.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "fixedPower.H"
#include "addToRunTimeSelectionTable.H"
#include "structure.H"
#include "fvMatrix.H"
#include "fvmDdt.H"
#include "fvmSup.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace powerModels
{
    defineTypeNameAndDebug(fixedPower, 0);
    addToRunTimeSelectionTable
    (
        powerModel, 
        fixedPower, 
        powerModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::powerModels::fixedPower::fixedPower
(
    structure& structureRef,
    const dictionary& dicts
)
:
    powerModel
    (
        structureRef,
        dicts
    ),
    powerDensity_
    (
        IOobject
        (
            IOobject::groupName("powerDensity", typeName),
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimPower/dimVol, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    timeProfile_(this->toc().size()),
    t0_(this->toc().size()),
    timeDependent_(this->toc().size()),
    T_
    (
        IOobject
        (
            IOobject::groupName("T", typeName),
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimTemperature, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    alphaRhoCp_
    (
        IOobject
        (
            IOobject::groupName("alphaRhoCp", typeName),
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("", dimEnergy/dimVol/dimTemperature, 0.0),
        zeroGradientFvPatchScalarField::typeName
    )
{
    this->setInterfacialArea();
    
    structure_.setRegionField(*this, powerDensity_, "powerDensity");
    structure_.setRegionField(*this, T_, "T");

    forAll(this->toc(), regioni)
    {
        word region(this->toc()[regioni]);
        const dictionary& dict(this->subDict(region));

        //- Setup cellToRegion_ mapping
        const labelList& regionCells
        (
            structure_.cellLists()[region]
        );
        scalar rhoCp(0);
        if (dict.found("rho") and dict.found("Cp"))
        {
            rhoCp =
                dict.get<scalar>("rho")*
                dict.get<scalar>("Cp");
        }
        else if (dict.found("rhoCp"))
        {
            rhoCp = dict.get<scalar>("rhoCp");
        }
        else
        {
            FatalErrorInFunction
                << "fixedPower region: " << region << " -> "
                << "specify either rhoCp or both rho and Cp"
                << exit(FatalError);
        }
        forAll(regionCells, i)
        {
            label celli(regionCells[i]);
            alphaRhoCp_[celli] = rhoCp;
        }

        // Preparing data to update power
        word timeProfileDictName("powerTimeProfile");
        timeDependent_[regioni] = false;
        if (dict.found(timeProfileDictName))
        {
            
            const dictionary& timeProfileDict(dict.subDict(timeProfileDictName));
            word type
            (
                timeProfileDict.get<word>("type")
            );

            timeProfile_.set        
            (
                regioni,
                Function1<scalar>::New
                (
                    type,
                    timeProfileDict,
                    type
                )
            );
            timeDependent_[regioni] = true;
            t0_[regioni] = timeProfileDict.lookupOrDefault("startTime", 0.0);

        }

    }

    //- The alphaRhoCp is read as a rhoCp, alpha is multiplied at this step
    alphaRhoCp_ = 
        max
        (
            alpha_*alphaRhoCp_, 
            dimensionedScalar
            ("", dimEnergy/dimVol/dimTemperature, 1e-69)
        );
    alphaRhoCp_.correctBoundaryConditions();


}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::powerModels::fixedPower::~fixedPower()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::powerModels::fixedPower::correctT(volScalarField& T) const
{
    forAll(cellList_, i)
    {
        label celli(cellList_[i]);
        T[celli] = T_[celli];
    }
    //- Correct boundary conditions called on T in structure.C after it has
    //  been set, cell-by-cell, by all powerModels
}


void Foam::powerModels::fixedPower::powerOff()
{
    powerDensity_ *= 0.0;
}

void Foam::powerModels::fixedPower::powerUpdate()
{
    forAll(this->toc(), regioni)
    {
        word region(this->toc()[regioni]);
        
        //- Setup cellToRegion_ mapping
        const labelList& regionCells
        (
            structure_.cellLists()[region]
        );
        if(timeDependent_[regioni])
        {
            scalar t(mesh_.time().timeOutputValue()-t0_[regioni]);
            scalar timeDependentPowerDensity(timeProfile_[regioni].value(t));
            forAll(regionCells, i)
            {
                label celli(regionCells[i]);
                powerDensity_[celli] =  timeDependentPowerDensity;
            }          
        }
    }
}

void Foam::powerModels::fixedPower::correct
(
    const volScalarField& HTSum,  // == SUM_j [htc_j*T_j*frac_j]
    const volScalarField& HSum    // == SUM_j [htc_j*frac_j]
)
{
    this->powerUpdate();
    scalar dt(mesh_.time().deltaT().value());
    volScalarField& T0(T_.oldTime()); 
    
    forAll(cellList_, i)
    {
        label celli(cellList_[i]);
        scalar alphaRhoCpByDt(alphaRhoCp_[celli]/dt);
        const scalar& iA(iA_[celli]);
        T_[celli] = 
            (
                iA*HTSum[celli] 
            +   alpha_[celli]*powerDensity_[celli]
            +   alphaRhoCpByDt*T0[celli]
            )/
            (alphaRhoCpByDt + iA*HSum[celli]);
    }
    T_.correctBoundaryConditions();

    /*
    What I did above is equivalent to doing this, but possibly faster, given
    that I might not have to solve the equation for the entire mesh, but only
    over the cells over which the fixedPower models exists. The time
    derivative was hardcoded to be forward Euler

    fvScalarMatrix EEqn
    (
            fvm::ddt(alphaRhoCp_, T_)
        ==
            iA_*HTSum
        -   fvm::Sp(iA_*HSum, T_)
        +   alpha_*powerDensity_
    );
    EEqn.solve();
    T_.correctBoundaryConditions();
    */
}


// ************************************************************************* //
