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

#include "porousKEpsilon.H"
#include "fvOptions.H"
#include "bound.H"

#include "fluid.H"
#include "structure.H"
#include "FSPair.H"
#include "myOps.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void porousKEpsilon<BasicTurbulenceModel>::correctNut()
{
    this->nut_ = Cmu_*sqr(k_)/epsilon_;
    this->nut_.correctBoundaryConditions();

    //- Correct alphat before nut is stabilized
    this->Prt_ = dimensioned<scalar>::lookupOrDefault
    (
        "Prt",
        this->coeffDict(),
        1.0
    );
    this->alphat_ = this->rho_*this->nut_/this->Prt_;
    this->alphat_.correctBoundaryConditions();

    if (nutStabilization_)
    {
        this->nut_ += 
            pos(structure_)*FSPair_.fluidRef().magU()*DhStructPtr_()/
            laminarReStructPtr_();
        this->nut_.correctBoundaryConditions();
    }

    fv::options::New(this->mesh_).correct(this->nut_);

    // BasicTurbulenceModel::correctNut();  //- Eh, I don't want other things
                                            //  messing with the stabilized nut
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> porousKEpsilon<BasicTurbulenceModel>::kSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            k_,
            dimVolume*this->rho_.dimensions()*k_.dimensions()
            /dimTime
        )
    );
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> porousKEpsilon<BasicTurbulenceModel>::epsilonSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            epsilon_,
            dimVolume*this->rho_.dimensions()*epsilon_.dimensions()
            /dimTime
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
porousKEpsilon<BasicTurbulenceModel>::porousKEpsilon
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    eddyViscosity<RASModel<BasicTurbulenceModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            this->coeffDict_,
            0.09
        )
    ),
    C1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1",
            this->coeffDict_,
            1.44
        )
    ),
    C2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2",
            this->coeffDict_,
            1.92
        )
    ),
    C3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C3",
            this->coeffDict_,
            0
        )
    ),
    sigmak_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmak",
            this->coeffDict_,
            1.0
        )
    ),
    sigmaEps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaEps",
            this->coeffDict_,
            1.3
        )
    ),

    k_
    (
        IOobject
        (
            IOobject::groupName("k", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    epsilon_
    (
        IOobject
        (
            IOobject::groupName("epsilon", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    fluidName_
    (
        alpha.name() == "alpha" ? 
        "" : myOps::split<word>(alpha.name(), '.')[1]
    ),
    structure_
    (
        this->mesh_.objectRegistry::template 
            lookupObject<structure>("alpha.structure")
    ),
    FSPair_
    (
        this->mesh_.objectRegistry::template 
            foundObject<FSPair>(fluidName_+".structure") ?
        this->mesh_.objectRegistry::template 
            lookupObject<FSPair>(fluidName_+".structure") :
        this->mesh_.objectRegistry::template 
            lookupObject<FSPair>("fluid.structure")
    ),
    porousKEpsilonDict_
    (
        this->subDict("porousKEpsilonProperties")
    ),
    convergenceLength_
    (
        IOobject
        (
            IOobject::groupName("convergenceLength", fluidName_),
            this->mesh_.time().timeName(),
            this->mesh_
        ),
        this->mesh_,
        dimensionedScalar("", dimLength, 1e69), //- Size matters lol
        zeroGradientFvPatchScalarField::typeName
    ),
    turbulenceIntensityCoeff_
    (
        IOobject
        (
            IOobject::groupName("turbulenceIntensityCoeff", fluidName_),
            this->mesh_.time().timeName(),
            this->mesh_
        ),
        this->mesh_,
        dimensionedScalar("", dimless, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    turbulenceIntensityExp_
    (
        IOobject
        (
            IOobject::groupName("turbulenceIntensityExp", fluidName_),
            this->mesh_.time().timeName(),
            this->mesh_
        ),
        this->mesh_,
        dimensionedScalar("", dimless, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    turbulenceLengthScaleCoeff_
    (
        IOobject
        (
            IOobject::groupName("turbulenceLengthScaleCoeff", fluidName_),
            this->mesh_.time().timeName(),
            this->mesh_
        ),
        this->mesh_,
        dimensionedScalar("", dimless, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    equilibriumEpsilon_
    (
        IOobject
        (
            IOobject::groupName("equilibriumEpsilon", fluidName_),
            this->mesh_.time().timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("", dimArea/dimTime/dimTime/dimTime, 0)
    ),
    equilibriumK_
    (
        IOobject
        (
            IOobject::groupName("equilibriumK", fluidName_),
            this->mesh_.time().timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("", dimArea/dimTime/dimTime, 0)
    ),
    nutStabilization_(false),
    Cmu3by4_(pow(Cmu_, 0.75))
{
    bound(k_, this->kMin_);
    bound(epsilon_, this->epsilonMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }

    wordList regionsInDict(0);
    
    forAll(porousKEpsilonDict_.toc(), i)
    {
        word key(porousKEpsilonDict_.toc()[i]);
        const dictionary& regionDict(porousKEpsilonDict_.subDict(key));
        
        wordList regions(myOps::split<word>(key, ':'));

        forAll(regions, j)
        {
            word region(regions[j]);
            
            //- If the region volumeFraction is 0, there is no structure so 
            //  regular kEpsilon applies and there is no need to set 
            //  porousKEpsilon properties
            if (max(structure_.alphaFields()[region]).value() == 0.0) continue;
            
            Info<< "Setting porousKEpsilon parameters for region: " << region 
            << endl;

            regionsInDict.append(region);

            const labelList& cellList(structure_.cellLists()[region]);

            scalar convergenceLength
            (
                regionDict.get<scalar>("convergenceLength")
            );
            scalar turbulenceIntensityCoeff
            (
                regionDict.get<scalar>("turbulenceIntensityCoeff")
            );
            scalar turbulenceIntensityExp
            (
                regionDict.get<scalar>("turbulenceIntensityExp")
            );
            scalar turbulenceLengthScaleCoeff
            (
                regionDict.get<scalar>("turbulenceLengthScaleCoeff")
            );

            forAll(cellList, k)
            {
                label celli(cellList[k]);
                convergenceLength_[celli] = convergenceLength;
                turbulenceIntensityCoeff_[celli] = turbulenceIntensityCoeff;
                turbulenceIntensityExp_[celli] = turbulenceIntensityExp;
                turbulenceLengthScaleCoeff_[celli] = 
                    turbulenceLengthScaleCoeff;
            }
            convergenceLength_.correctBoundaryConditions();
            turbulenceIntensityCoeff_.correctBoundaryConditions();
            turbulenceIntensityExp_.correctBoundaryConditions();
            turbulenceLengthScaleCoeff_.correctBoundaryConditions();  

            //- nut stabilization fields, create only if keyword found
            if (regionDict.found("DhStruct"))
            {
                scalar defaultLaminarReStruct(500);
                nutStabilization_ = true;
                if (!DhStructPtr_.valid())
                {
                    DhStructPtr_.reset
                    (
                        new volScalarField
                        (
                            IOobject
                            (
                                "DhStructure.porousKEpsilon",
                                this->mesh_.time().timeName(),
                                this->mesh_
                            ),
                            this->mesh_,
                            dimensionedScalar("", dimLength, 0.0),
                            zeroGradientFvPatchScalarField::typeName
                        )
                    );
                }
                if (!laminarReStructPtr_.valid())
                {
                    laminarReStructPtr_.reset
                    (
                        new volScalarField
                        (
                            IOobject
                            (
                                "laminarReStructure.porousKEpsilon",
                                this->mesh_.time().timeName(),
                                this->mesh_
                            ),
                            this->mesh_,
                            dimensionedScalar
                            (
                                "", 
                                dimless, 
                                defaultLaminarReStruct
                            ),
                            zeroGradientFvPatchScalarField::typeName
                        )
                    );
                }
                scalar DhStruct(regionDict.get<scalar>("DhStruct"));
                scalar laminarReStruct
                (
                    regionDict.lookupOrDefault<scalar>
                    (
                        "laminarReStruct", 
                        defaultLaminarReStruct
                    )
                );
                forAll(cellList, k)
                {
                    label celli(cellList[k]);
                    DhStructPtr_()[celli] = DhStruct;
                    laminarReStructPtr_()[celli] = laminarReStruct;
                }
                DhStructPtr_().correctBoundaryConditions();
                laminarReStructPtr_().correctBoundaryConditions();
            }  
        }
    }

    //- Check that porousKEpsilon properties were set for all porous regions
    forAll(structure_.regions(), i)
    {
        word region(structure_.regions()[i]);

        //- If the region volumeFraction is 0, there is no structure so regular
        //  kEpsilon applies and there is no need to check for porousKEpsilon
        //  properties
        if (max(structure_.alphaFields()[region]).value() == 0.0) continue;
        
        bool found(false);
        forAll(regionsInDict, j)
        {
            if (regionsInDict[j] == region)
            {
                found = true;
                break;
            }
        }

        if (!found)
        {
            FatalErrorInFunction
                << "porousKEpsilon region: " << region << " -> "
                << "properties were not set / subDictionary not found!"
                << exit(FatalError);
        }
    }

    convergenceLength_ = 
        max
        (
            convergenceLength_, 
            dimensionedScalar("", dimLength, SMALL)
        );
    convergenceLength_.correctBoundaryConditions();

    Info << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool porousKEpsilon<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
        Cmu_.readIfPresent(this->coeffDict());
        C1_.readIfPresent(this->coeffDict());
        C2_.readIfPresent(this->coeffDict());
        C3_.readIfPresent(this->coeffDict());
        sigmak_.readIfPresent(this->coeffDict());
        sigmaEps_.readIfPresent(this->coeffDict());

        Cmu3by4_ = pow(Cmu_, 0.75);

        return true;
    }

    return false;
}


template<class BasicTurbulenceModel>
void porousKEpsilon<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();

    volScalarField::Internal divU
    (
        fvc::div(fvc::absolute(this->phi(), U))().v()
    );

    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField::Internal G
    (
        this->GName(),
        nut.v()*(dev(twoSymm(tgradU().v())) && tgradU().v())
    );
    tgradU.clear();

    // Update epsilon and G at the wall
    epsilon_.boundaryFieldRef().updateCoeffs();

    //- Structure marker fields
    tmp<volScalarField> noStructure(neg0(structure_));
    tmp<volScalarField> yesStructure(pos(structure_));

    //- Compute equilibrium k, epsilon values
    equilibriumK_ =
    (
        1.5*sqr
        (
            mag(U)*
            turbulenceIntensityCoeff_*pow
            (
                FSPair_.Re(), turbulenceIntensityExp_
            )
        )
    );
    //equilibriumEpsilon_.correctBoundaryConditions();
    equilibriumEpsilon_ =
    ( 
        Cmu3by4_*pow(equilibriumK_, 1.5)
        /
        (
            max
            (
                turbulenceLengthScaleCoeff_*structure_.Dh(),
                dimensionedScalar("", dimLength, SMALL)
            )
        )
    );
    //equilibriumK_.correctBoundaryConditions();

    /*
    Info<< "eqEpsilon (avg min max) = "
        << equilibriumEpsilon_.weightedAverage(this->mesh_.V()).value()
        << " " << min(equilibriumEpsilon_).value()
        << " " << max(equilibriumEpsilon_).value()
        << " " << equilibriumEpsilon_.dimensions() << endl;
    Info<< "eqK (avg min max) = "
        << equilibriumK_.weightedAverage(this->mesh_.V()).value()
        << " " << min(equilibriumK_).value()
        << " " << max(equilibriumK_).value()
        << " " << equilibriumK_.dimensions() << endl;
    */

    //- magU/convergenceLength gives the inverse of the convergence time scale,
    //  i.e. what was previously know/provided as convergenceRate
    tmp<volScalarField> alphaRhoConv
    (
        alpha*rho*(mag(U)/convergenceLength_)*yesStructure()
    );

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(alpha, rho, epsilon_)
      + fvm::div(alphaRhoPhi, epsilon_)
      - fvm::laplacian(alpha*rho*DepsilonEff(), epsilon_)
     ==
        C1_*alpha()*rho()*G*epsilon_()/k_()*noStructure()
      - fvm::SuSp
        (
            ((2.0/3.0)*C1_ - C3_)*alpha()*rho()*divU*noStructure(), 
            epsilon_
        )
      - fvm::Sp(C2_*alpha()*rho()*epsilon_()/k_()*noStructure(), epsilon_)
      + epsilonSource()
      - fvm::Sp(alphaRhoConv(), epsilon_)
      + alphaRhoConv()*equilibriumEpsilon_
      + fvOptions(alpha, rho, epsilon_)
    );

    epsEqn.ref().relax();
    fvOptions.constrain(epsEqn.ref());
    epsEqn.ref().boundaryManipulate(epsilon_.boundaryFieldRef());
    solve(epsEqn);
    fvOptions.correct(epsilon_);
    bound(epsilon_, this->epsilonMin_);

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(), k_)
     ==
        alpha()*rho()*G*noStructure()
      - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU*noStructure(), k_)
      - fvm::Sp(alpha()*rho()*epsilon_()/k_()*noStructure(), k_)
      + kSource()
      - fvm::Sp(alphaRhoConv(), k_)
      + alphaRhoConv()*equilibriumK_
      + fvOptions(alpha, rho, k_)
    );
    
    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);

    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
