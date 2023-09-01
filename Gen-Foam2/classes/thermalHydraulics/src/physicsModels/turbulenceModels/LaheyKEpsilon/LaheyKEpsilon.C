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

#include "LaheyKEpsilon.H"
#include "fvOptions.H"
#include "zeroGradientFvPatchFields.H"
#include "fluid.H"
#include "FFPair.H"
#include "myOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
LaheyKEpsilon<BasicTurbulenceModel>::LaheyKEpsilon
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
    kEpsilon<BasicTurbulenceModel>
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName,
        type
    ),
    gasPtr_(nullptr),
    liquidPtr_(nullptr),
    pairPtr_(nullptr),
    alphaInversion_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaInversion",
            this->coeffDict_,
            0.3
        )
    ),
    Cp_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cp",
            this->coeffDict_,
            0.25
        )
    ),
    C3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C3",
            this->coeffDict_,
            this->C2_.value()
        )
    ),
    Cmub_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmub",
            this->coeffDict_,
            0.6
        )
    )
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }

    setFluidNames();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void LaheyKEpsilon<BasicTurbulenceModel>::setFluidNames() const
{
    const fvMesh& mesh(this->mesh_);
    HashTable<const fluid*> fluids(mesh.lookupClass<const fluid>());
    const fluid& fluid1 = *(fluids[fluids.toc()[0]]);
    const fluid& fluid2 = *(fluids[fluids.toc()[1]]);
    if 
    (
        !(fluid1.isLiquid() and fluid2.isGas()) and
        !(fluid2.isLiquid() and fluid1.isGas())
    )
    {
        FatalErrorInFunction
            << "The LaheyKEpsilon turbulence model only works for "
            << "liquid-gas systems. Set the  the stateOfMatter entry in "
            << "phaseProperties." << fluid1.name() << "Properties and/or "
            << "phaseProperties." << fluid2.name() << "Properties) to "
            << "distinguish between gas and liquid"
            << exit(FatalError);
    }
    liquidName_ = (fluid1.isLiquid()) ? fluid1.name() : fluid2.name();
    gasName_ = (fluid1.isGas()) ? fluid1.name() : fluid2.name();
}


template<class BasicTurbulenceModel>
const Foam::fluid& LaheyKEpsilon<BasicTurbulenceModel>::gas() const
{
    if (!gasPtr_)
    {
        const fvMesh& mesh(this->mesh_);
        gasPtr_ = &(mesh.lookupObject<fluid>("alpha."+gasName_));
    }
    return *gasPtr_;
}


template<class BasicTurbulenceModel>
const Foam::fluid& LaheyKEpsilon<BasicTurbulenceModel>::liquid() const
{
    if (!liquidPtr_)
    {
        const fvMesh& mesh(this->mesh_);
        liquidPtr_ = 
            &(mesh.lookupObject<fluid>("alpha."+liquidName_));
    }
    return *liquidPtr_;
}


template<class BasicTurbulenceModel>
const Foam::FFPair& LaheyKEpsilon<BasicTurbulenceModel>::pair() const
{
    if (!pairPtr_)
    {
        const fvMesh& mesh(this->mesh_);
        word keyLG(IOobject::groupName(liquidName_, gasName_));
        word keyGL(IOobject::groupName(gasName_, liquidName_));
        pairPtr_ = 
            &(
                (mesh.foundObject<FFPair>(keyLG)) ?
                mesh.lookupObject<FFPair>(keyLG) :
                mesh.lookupObject<FFPair>(keyGL)
            );
    }
    return *pairPtr_;
}


template<class BasicTurbulenceModel>
const Foam::tmp<Foam::volScalarField> LaheyKEpsilon<BasicTurbulenceModel>::Cd()
const
{
    tmp<volScalarField> tCd
    (
        new volScalarField
        (
            IOobject
            (
                "Cd",
                this->mesh_.time().timeName(),
                this->mesh_
            ),
            this->mesh_,
            dimensionedScalar("", dimless, 1e-3),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    volScalarField& Cd = tCd.ref();

    const FFPair& p(pair());
    const volScalarField& l(liquid());
    const volScalarField& g(gas());

    //- Compute Cd from Kd, cell-by-cell as it's faster (I don't really care 
    //  about BCs)
    forAll(Cd, i)
    {
        const scalar& li(l[i]);
        const scalar& gi(g[i]);
        Cd[i] = 
            (2.0)*p.Kd()[i]*p.DhDispersed()[i]/p.rhoContinuous()[i]/
            max
            (
                (li*gi)/(li+gi)*p.magUr()[i], 1e-3
            );
    }
    Cd.correctBoundaryConditions();

    return tCd;
}


template<class BasicTurbulenceModel>
bool LaheyKEpsilon<BasicTurbulenceModel>::read()
{
    if (kEpsilon<BasicTurbulenceModel>::read())
    {
        alphaInversion_.readIfPresent(this->coeffDict());
        Cp_.readIfPresent(this->coeffDict());
        C3_.readIfPresent(this->coeffDict());
        Cmub_.readIfPresent(this->coeffDict());

        return true;
    }

    return false;
}


template<class BasicTurbulenceModel>
void LaheyKEpsilon<BasicTurbulenceModel>::correctNut()
{
    this->nut_ =
        this->Cmu_*sqr(this->k_)/this->epsilon_
      + Cmub_*gas().Dh()*gas()
       *(mag(pair().magUr()));

    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


template<class BasicTurbulenceModel>
tmp<volScalarField> LaheyKEpsilon<BasicTurbulenceModel>::bubbleG() const
{
    tmp<volScalarField> bubbleG
    (
        Cp_*
        (
            1.0
        +   pow(Cd(), 4.0/3.0)
        )*
        gas()*pow3(pair().magUr())/
        max(gas().Dh(), dimensionedScalar("", dimLength, 1e-3))
    );

    return bubbleG;
}


template<class BasicTurbulenceModel>
tmp<volScalarField>
LaheyKEpsilon<BasicTurbulenceModel>::phaseTransferCoeff() const
{
    const volVectorField& U = this->U_;
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;

    const phaseCompressibleTurbulenceModel& gasTurbulence = gas().turbulence();

    return
    (
        max(alphaInversion_ - alpha, scalar(0))
       *rho
       *min(gasTurbulence.epsilon()/gasTurbulence.k(), 1.0/U.time().deltaT())
    );
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> LaheyKEpsilon<BasicTurbulenceModel>::kSource() const
{
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;

    volScalarField pTC(phaseTransferCoeff());

    return
        alpha*rho*bubbleG()
      + pTC*gas().turbulence().k()
      - fvm::Sp(pTC, this->k_);
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> LaheyKEpsilon<BasicTurbulenceModel>::epsilonSource() const
{
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;

    volScalarField pTC(phaseTransferCoeff());

    return
        alpha*rho*this->C3_*this->epsilon_*bubbleG()/this->k_
      + pTC*gas().turbulence().epsilon()
      - fvm::Sp(pTC, this->epsilon_);
}


template<class BasicTurbulenceModel>
void LaheyKEpsilon<BasicTurbulenceModel>::correct()
{
    kEpsilon<BasicTurbulenceModel>::correct();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
