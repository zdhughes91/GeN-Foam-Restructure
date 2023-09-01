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

#include "mixtureKEpsilon.H"
#include "fvOptions.H"
#include "bound.H"
#include "fixedValueFvPatchFields.H"
#include "inletOutletFvPatchFields.H"
#include "fvmSup.H"

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
mixtureKEpsilon<BasicTurbulenceModel>::mixtureKEpsilon
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
    isGas_(false),
    liquidPtr_(nullptr),
    gasPtr_(nullptr),
    pairPtr_(nullptr),
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
            C2_.value()
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
    )
{
    bound(k_, this->kMin_);
    bound(epsilon_, this->epsilonMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }

    setFluidNames();
    word thisName(myOps::split<word>(alpha.name(), '.')[1]);
    isGas_ = (thisName == gasName_);
}


template<class BasicTurbulenceModel>
wordList mixtureKEpsilon<BasicTurbulenceModel>::epsilonBoundaryTypes
(
    const volScalarField& epsilon
) const
{
    const volScalarField::Boundary& ebf = epsilon.boundaryField();

    wordList ebt = ebf.types();

    forAll(ebf, patchi)
    {
        if (isA<fixedValueFvPatchScalarField>(ebf[patchi]))
        {
            ebt[patchi] = fixedValueFvPatchScalarField::typeName;
        }
    }

    return ebt;
}


template<class BasicTurbulenceModel>
void mixtureKEpsilon<BasicTurbulenceModel>::correctInletOutlet
(
    volScalarField& vsf,
    const volScalarField& refVsf
) const
{
    volScalarField::Boundary& bf = vsf.boundaryFieldRef();
    const volScalarField::Boundary& refBf =
        refVsf.boundaryField();

    forAll(bf, patchi)
    {
        if
        (
            isA<inletOutletFvPatchScalarField>(bf[patchi])
         && isA<inletOutletFvPatchScalarField>(refBf[patchi])
        )
        {
            refCast<inletOutletFvPatchScalarField>
            (bf[patchi]).refValue() =
            refCast<const inletOutletFvPatchScalarField>
            (refBf[patchi]).refValue();
        }
    }
}


template<class BasicTurbulenceModel>
void mixtureKEpsilon<BasicTurbulenceModel>::initMixtureFields()
{
    if (rhom_.valid()) return;

    // Local references to gas-phase properties
    const mixtureKEpsilon<BasicTurbulenceModel>& turbd = 
        refCast<const mixtureKEpsilon<BasicTurbulenceModel>>
        (
            gas().turbulence()
        );
    const volScalarField& kg = turbd.k();
    const volScalarField& epsilong = turbd.epsilon();

    // Local references to liquid-phase properties
    const mixtureKEpsilon<BasicTurbulenceModel>& turbc = 
        refCast<const mixtureKEpsilon<BasicTurbulenceModel>>
        (
            liquid().turbulence()
        );
    const volScalarField& kl = turbc.k();
    const volScalarField& epsilonl = turbc.epsilon();

    word startTimeName
    (
        this->runTime_.timeName(this->runTime_.startTime().value())
    );

    Ct2_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "Ct2",
                startTimeName,
                this->mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            Ct2()
        )
    );

    rhom_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "rhom",
                startTimeName,
                this->mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            rhom()
        )
    );

    km_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "km",
                startTimeName,
                this->mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mix(kl, kg),
            kl.boundaryField().types()
        )
    );
    correctInletOutlet(km_(), kl);

    epsilonm_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "epsilonm",
                startTimeName,
                this->mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mix(epsilonl, epsilong),
            epsilonBoundaryTypes(epsilonl)
        )
    );
    correctInletOutlet(epsilonm_(), epsilonl);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void mixtureKEpsilon<BasicTurbulenceModel>::setFluidNames() const
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
            << "The mixtureKEpsilon turbulence model only works for "
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
const Foam::fluid& mixtureKEpsilon<BasicTurbulenceModel>::gas() const
{
    if (!gasPtr_)
    {
        const fvMesh& mesh(this->mesh_);
        gasPtr_ = &(mesh.lookupObject<fluid>("alpha."+gasName_));
    }
    return *gasPtr_;
}


template<class BasicTurbulenceModel>
const Foam::fluid& mixtureKEpsilon<BasicTurbulenceModel>::liquid() const
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
const Foam::FFPair& mixtureKEpsilon<BasicTurbulenceModel>::pair() const
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
const Foam::tmp<Foam::volScalarField> mixtureKEpsilon<BasicTurbulenceModel>::Cd()
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
bool mixtureKEpsilon<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
        Cmu_.readIfPresent(this->coeffDict());
        C1_.readIfPresent(this->coeffDict());
        C2_.readIfPresent(this->coeffDict());
        C3_.readIfPresent(this->coeffDict());
        Cp_.readIfPresent(this->coeffDict());
        sigmak_.readIfPresent(this->coeffDict());
        sigmaEps_.readIfPresent(this->coeffDict());

        return true;
    }

    return false;
}


template<class BasicTurbulenceModel>
void mixtureKEpsilon<BasicTurbulenceModel>::correctNut()
{
    this->nut_ = Cmu_*sqr(k_)/epsilon_;
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


template<class BasicTurbulenceModel>
tmp<volScalarField> mixtureKEpsilon<BasicTurbulenceModel>::Ct2() const
{
    volScalarField beta
    (
        (6*this->Cmu_/(4*sqrt(3.0/2.0)))*
        pair().Kd()/liquid().rho()*
        (liquid().turbulence().k()/liquid().turbulence().epsilon())
    );
    volScalarField Ct0
    (
        (3 + beta)/(1 + beta + 2*gas().rho()/liquid().rho())
    );
    volScalarField fAlphad((180 + (-4.71e3 + 4.26e4*gas())*gas())*gas());

    return sqr(1 + (Ct0 - 1)*exp(-fAlphad));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> mixtureKEpsilon<BasicTurbulenceModel>::rholEff() const
{
    return liquid().rho();
}


template<class BasicTurbulenceModel>
tmp<volScalarField> mixtureKEpsilon<BasicTurbulenceModel>::rhogEff() const
{
    return
        gas().rho()
      + pair().Vm()*liquid().rho();
}


template<class BasicTurbulenceModel>
tmp<volScalarField> mixtureKEpsilon<BasicTurbulenceModel>::rhom() const
{
    return liquid()*rholEff() + gas()*rhogEff();
}


template<class BasicTurbulenceModel>
tmp<volScalarField> mixtureKEpsilon<BasicTurbulenceModel>::mix
(
    const volScalarField& fc,
    const volScalarField& fd
) const
{
    return (liquid()*rholEff()*fc + gas()*rhogEff()*fd)/rhom_();
}


template<class BasicTurbulenceModel>
tmp<volScalarField> mixtureKEpsilon<BasicTurbulenceModel>::mixU
(
    const volScalarField& fc,
    const volScalarField& fd
) const
{
    return
        (liquid()*rholEff()*fc + gas()*rhogEff()*Ct2_()*fd)/
        (liquid()*rholEff() + gas()*rhogEff()*Ct2_());
}


template<class BasicTurbulenceModel>
tmp<surfaceScalarField> mixtureKEpsilon<BasicTurbulenceModel>::mixFlux
(
    const surfaceScalarField& fc,
    const surfaceScalarField& fd
) const
{
    surfaceScalarField alphalf(fvc::interpolate(liquid()));
    surfaceScalarField alphagf(fvc::interpolate(gas()));

    surfaceScalarField rholEfff(fvc::interpolate(rholEff()));
    surfaceScalarField rhogEfff(fvc::interpolate(rhogEff()));

    return
       (alphalf*rholEfff*fc + alphagf*rhogEfff*fvc::interpolate(Ct2_())*fd)/
       (alphalf*rholEfff + alphagf*rhogEfff*fvc::interpolate(Ct2_()));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> mixtureKEpsilon<BasicTurbulenceModel>::bubbleG() const
{
    tmp<volScalarField> bubbleG
    (
        Cp_*
        //- Differs from the Lahey model as it has this extra term (which also
        //  makes them dimensionally different, but hey, I didn't come up
        //  with this, just porting stuff over to FFSEulerFoam)
        //   |
        //  \|/
        //   *
        liquid()*liquid().rho()*
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
tmp<fvScalarMatrix> mixtureKEpsilon<BasicTurbulenceModel>::kSource() const
{
    return fvm::Su(bubbleG()/rhom_(), km_());
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> mixtureKEpsilon<BasicTurbulenceModel>::epsilonSource() const
{
    return fvm::Su(C3_*epsilonm_()*bubbleG()/(rhom_()*km_()), epsilonm_());
}


template<class BasicTurbulenceModel>
void mixtureKEpsilon<BasicTurbulenceModel>::correct()
{
    // Only solve the mixture turbulence for the gas-phase
    if (!isGas_)
    {
        //- The consistency check is now done for both phases at every call
        //  of initMixtureFields, no need to do it here too
        return;
    }

    if (!this->turbulence_)
    {
        return;
    }

    // Initialise the mixture fields if they have not yet been constructed
    initMixtureFields();

    //- Note: I could modify all this to comply with the new structure
    //  (i.e. using gas() and liquid() for more symmetry and code clarity)
    //  but there is no practical need to and I'd only risk breaking things,
    //  so I won't

    // Local references to gas-phase properties
    tmp<surfaceScalarField> phig = this->phi();
    const volVectorField& Ug = this->U_;
    const volScalarField& alphag = this->alpha_;
    volScalarField& kg = this->k_;
    volScalarField& epsilong = this->epsilon_;
    volScalarField& nutg = this->nut_;

    // Local references to liquid-phase properties
    //- This is a werid way of getting a ref to liquid().turbulence(). Why
    //  so? First, I need it non-const accessible, so I need const_case.
    //  Secondly, I cannot just cast liquid().turbulence() as the casting
    //  fails, so I look for the turbulence model in the mesh registry as
    //  well
    const fvMesh& mesh(this->mesh_);
    mixtureKEpsilon<BasicTurbulenceModel>& liquidTurbulence =
        const_cast<mixtureKEpsilon<BasicTurbulenceModel>&>
        (
            mesh.lookupObject<mixtureKEpsilon<BasicTurbulenceModel>>
            (
                IOobject::groupName
                (
                    turbulenceModel::propertiesName,
                    liquidName_
                )
            )
        );
    tmp<surfaceScalarField> phil = liquidTurbulence.phi();
    const volVectorField& Ul = liquidTurbulence.U_;
    const volScalarField& alphal = liquidTurbulence.alpha_;
    volScalarField& kl = liquidTurbulence.k_;
    volScalarField& epsilonl = liquidTurbulence.epsilon_;
    volScalarField& nutl = liquidTurbulence.nut_;

    // Local references to mixture properties
    volScalarField& rhom = rhom_();
    volScalarField& km = km_();
    volScalarField& epsilonm = epsilonm_();

    fv::options& fvOptions(fv::options::New(this->mesh_));

    eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();

    // Update the effective mixture density
    rhom = this->rhom();

    // Mixture flux
    surfaceScalarField phim("phim", mixFlux(phil, phig));

    // Mixture velocity divergence
    volScalarField divUm
    (
        mixU
        (
            fvc::div(fvc::absolute(phil, Ul)),
            fvc::div(fvc::absolute(phig, Ug))
        )
    );

    tmp<volScalarField> Gc;
    {
        tmp<volTensorField> tgradUl = fvc::grad(Ul);
        Gc = tmp<volScalarField>
        (
            new volScalarField
            (
                this->GName(),
                nutl*(tgradUl() && dev(twoSymm(tgradUl())))
            )
        );
        tgradUl.clear();

        // Update k, epsilon and G at the wall
        kl.boundaryFieldRef().updateCoeffs();
        epsilonl.boundaryFieldRef().updateCoeffs();

        Gc.ref().checkOut();
    }

    tmp<volScalarField> Gd;
    {
        tmp<volTensorField> tgradUg = fvc::grad(Ug);
        Gd = tmp<volScalarField>
        (
            new volScalarField
            (
                this->GName(),
                nutg*(tgradUg() && dev(twoSymm(tgradUg())))
            )
        );
        tgradUg.clear();

        // Update k, epsilon and G at the wall
        kg.boundaryFieldRef().updateCoeffs();
        epsilong.boundaryFieldRef().updateCoeffs();

        Gd.ref().checkOut();
    }

    // Mixture turbulence generation
    volScalarField Gm(mix(Gc, Gd));

    // Mixture turbulence viscosity
    volScalarField nutm(mixU(nutl, nutg));

    // Update the mixture k and epsilon boundary conditions
    km == mix(kl, kg);
    bound(km, this->kMin_);
    epsilonm == mix(epsilonl, epsilong);
    bound(epsilonm, this->epsilonMin_);

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(epsilonm)
      + fvm::div(phim, epsilonm)
      - fvm::Sp(fvc::div(phim), epsilonm)
      - fvm::laplacian(DepsilonEff(nutm), epsilonm)
     ==
        C1_*Gm*epsilonm/km
      - fvm::SuSp(((2.0/3.0)*C1_)*divUm, epsilonm)
      - fvm::Sp(C2_*epsilonm/km, epsilonm)
      + epsilonSource()
      + fvOptions(epsilonm)
    );

    epsEqn.ref().relax();
    fvOptions.constrain(epsEqn.ref());
    epsEqn.ref().boundaryManipulate(epsilonm.boundaryFieldRef());
    solve(epsEqn);
    fvOptions.correct(epsilonm);
    bound(epsilonm, this->epsilonMin_);


    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kmEqn
    (
        fvm::ddt(km)
      + fvm::div(phim, km)
      - fvm::Sp(fvc::div(phim), km)
      - fvm::laplacian(DkEff(nutm), km)
     ==
        Gm
      - fvm::SuSp((2.0/3.0)*divUm, km)
      - fvm::Sp(epsilonm/km, km)
      + kSource()
      + fvOptions(km)
    );

    kmEqn.ref().relax();
    fvOptions.constrain(kmEqn.ref());
    solve(kmEqn);
    fvOptions.correct(km);
    bound(km, this->kMin_);
    km.correctBoundaryConditions();

    volScalarField Cc2(rhom/(alphal*rholEff() + alphag*rhogEff()*Ct2_()));
    kl = Cc2*km;
    kl.correctBoundaryConditions();
    epsilonl = Cc2*epsilonm;
    epsilonl.correctBoundaryConditions();
    liquidTurbulence.correctNut();

    Ct2_() = Ct2();
    kg = Ct2_()*kl;
    kg.correctBoundaryConditions();
    epsilong = Ct2_()*epsilonl;
    epsilong.correctBoundaryConditions();
    nutg = Ct2_()*(liquidTurbulence.nu()/this->nu())*nutl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
