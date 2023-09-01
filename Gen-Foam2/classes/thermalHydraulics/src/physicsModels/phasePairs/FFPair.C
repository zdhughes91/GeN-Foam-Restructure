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

#include "FFPair.H"
#include "fvCFD.H"
#include "myOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(FFPair, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FFPair::FFPair
(
    fluid& fluid1,
    fluid& fluid2,
    const dictionary& dict
)
:
    regIOobject
    (
        IOobject
        (
            IOobject::groupName(fluid1.name(), fluid2.name()),
            fluid1.mesh().time().timeName(),
            fluid1.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true
        )
    ),
    mesh_(fluid1.mesh()),
    dict_(dict),
    pimple_(mesh_.lookupObject<customPimpleControl>("solutionControl")),
    fluid1_(fluid1),
    fluid2_(fluid2),
    alphaDispersed_
    (
        IOobject
        (
            "alphaDispersed",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    alphaContinuous_
    (
        IOobject
        (
            "alphaContinuous",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    DhDispersed_
    (
        IOobject
        (
            "DhDispersed",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("", dimLength, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    DhContinuous_
    (
        IOobject
        (
            "DhContinuous",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("", dimLength, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    rhoContinuous_
    (
        IOobject
        (
            "rhoContinuous",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("", dimDensity, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    nuContinuous_
    (
        IOobject
        (
            "nuContinuous",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("", dimArea/dimTime, SMALL),
        zeroGradientFvPatchScalarField::typeName
    ),
    PrContinuous_
    (
        IOobject
        (
            "PrContinuous",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    Re_
    (
        IOobject
        (
            IOobject::groupName("Re", this->name()),
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimless, 1),
        zeroGradientFvPatchScalarField::typeName
    ),
    minRe_
    (
        dimensionedScalar::lookupOrDefault
        (
            "residualFluidFluidRe",
            dict,
            dimless,
            10
        )
    ),
    magUr_
    (
        IOobject
        (
            "magUr",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("", dimVelocity, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    Kd_
    (
        IOobject
        (
            IOobject::groupName("Kd", fluid1_.name()+"."+fluid2_.name()),
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimDensity/dimTime, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    iT_
    (
        IOobject
        (
            IOobject::groupName("T", "interface"),
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimTemperature, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    iA_
    (
        IOobject
        (
            IOobject::groupName("areaDensity", "interface"),
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimArea/dimVolume, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    firstTimeStepAndIter_(true)
{
    Info << endl;
    
    //- Init htcs
    htcs_.set
    (
        fluid1_.name(),
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("htc", fluid1_.name()+".interface"),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar("", dimPower/dimArea/dimTemperature, 0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    htcs_.set
    (
        fluid2_.name(),
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("htc", fluid2_.name()+".interface"),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar("", dimPower/dimArea/dimTemperature, 0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    //- It is dimPower instead of dimPower/dimVolume as the fvScalarMatrix
    //  expects objects that are either already integrated over the volume
    //  (e.g. what the fvm::Su, SuSp, etc. return) or objects that are not
    //  already volume integrated, yet whose dimensions are equal to the
    //  dimension of the fvScalarMatrix (dimPower in this case) divided by
    //  dimVolume (e.g. explicit terms)
    heSources_.set
    (
        fluid1_.thermo().he().name(), 
        new fvScalarMatrix(fluid1_.thermo().he(), dimPower)
    );
    heSources_.set
    (
        fluid2_.thermo().he().name(), 
        new fvScalarMatrix(fluid2_.thermo().he(), dimPower)
    );

    const dictionary& physicsModelsDict(dict_.subDict("physicsModels"));

    //- Drag
    if (physicsModelsDict.found("dragModels"))
    {
        const dictionary& dragModelsDict
        (
            physicsModelsDict.subDict("dragModels")
        );
        if (foundUnorderedPairSubDict(dragModelsDict))
        {
            KdPtr_.reset
            (
                new FFDragFactor
                (
                    *this,
                    getUnorderedPairSubDict(dragModelsDict)
                )
            );
        }
    }

    //- Virtual mass
    if (physicsModelsDict.found("virtualMassCoefficientModel"))
    {
        const dictionary& virtualMassModelDict
        (
            physicsModelsDict.subDict("virtualMassCoefficientModel")
        );
        virtualMassPtr_.reset
        (
            new virtualMass
            (
                *this,
                virtualMassModelDict
            )
        );
    }

    //- Heat transfer
    if (physicsModelsDict.found("heatTransferModels"))
    {
        const dictionary& heatTransferModelsDict
        (
            physicsModelsDict.subDict("heatTransferModels")
        );
        if (foundUnorderedPairSubDict(heatTransferModelsDict))
        {
            const dictionary& pairDict
            (
                getUnorderedPairSubDict(heatTransferModelsDict)
            );
            if (pairDict.found(fluid1_.name()))
                htc1Ptr_.reset
                (
                    FFHeatTransferCoefficientModel::New
                    (
                        *this,
                        pairDict.subDict(fluid1_.name()),
                        mesh_
                    )
                );
            if (pairDict.found(fluid1_.name()))
                htc2Ptr_.reset
                (
                    FFHeatTransferCoefficientModel::New
                    (
                        *this,
                        pairDict.subDict(fluid2_.name()),
                        mesh_
                    )
                );
        }
    }

    //- Geometry models
    const dictionary& pairGeometryModels
    (
        physicsModelsDict.subDict("pairGeometryModels")
    );
    const dictionary& dispersionModelDict
    (
        getUnorderedPairSubDict(pairGeometryModels).subDict("dispersionModel")
    );
    dispersion1Ptr_.reset
    (
        dispersionModel::New
        (
            *this,
            dispersionModelDict,
            mesh_
        )
    );
    const dictionary& interfacialAreaModelDict
    (
        getUnorderedPairSubDict(pairGeometryModels).subDict
        (
            "interfacialAreaDensityModel"
        )
    );
    iAPtr_.reset
    (
        interfacialAreaModel::New
        (
            *this,
            interfacialAreaModelDict,
            mesh_
        )
    );

    //- Mass transfer
    if (physicsModelsDict.found("phaseChangeModel"))
    {
        const dictionary& phaseChangeModelDict
        (
            physicsModelsDict.subDict("phaseChangeModel")
        );
        phaseChangePtr_.reset
        (
            phaseChangeModel::New
            (
                *this,
                phaseChangeModelDict
            )
        );
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::FFPair::addVirtualMassForce(fvVectorMatrix& UEqn) const
{
    if (virtualMassPtr_.valid())
    {
        UEqn += virtualMassPtr_->force(UEqn.psi().name());
    }
}

void Foam::FFPair::correct
(
    const bool& correctFluidDynamics, 
    const bool& correctEnergy
)
{
    //- Init refs and dispersion marker fields
    const volVectorField& U1(fluid1_.U());
    const volVectorField& U2(fluid2_.U());
    const volScalarField& Dh1(fluid1_.Dh());
    const volScalarField& Dh2(fluid2_.Dh());

    //- Update dispersion/continuity
    scalarField oldDispersion(fluid1_.dispersion());
    dispersion1Ptr_->correctField(fluid1_.dispersion());
    scalar f = myOps::relaxationFactor(mesh_, "dispersion");
    if (f != 1.0)
        fluid1_.dispersion() = (1-f)*oldDispersion + f*fluid1_.dispersion();
    fluid2_.dispersion() = 1.0-fluid1_.dispersion();
    const scalarField& continuity1(fluid2_.dispersion());
    const scalarField& continuity2(fluid1_.dispersion());

    //- Update dipsersed/continuous phase fractions
    alphaDispersed_.primitiveFieldRef() = 
        fluid1_.primitiveField()*fluid1_.dispersion() 
    +   fluid2_.primitiveField()*fluid2_.dispersion();
    alphaDispersed_.correctBoundaryConditions();
    alphaContinuous_.primitiveFieldRef() = 
        fluid1_.primitiveField()*continuity1
    +   fluid2_.primitiveField()*continuity2;
    alphaContinuous_.correctBoundaryConditions();

    //- Update fluid diameter (needs to be done AFTER dispersion update)
    fluid1_.correctDiameter();
    fluid2_.correctDiameter();

    DhDispersed_.primitiveFieldRef() = 
        Dh1.primitiveField()*fluid1_.dispersion() 
    +   Dh2.primitiveField()*fluid2_.dispersion();
    DhDispersed_.correctBoundaryConditions();
    DhContinuous_.primitiveFieldRef() = 
        Dh1.primitiveField()*continuity1
    +   Dh2.primitiveField()*continuity2;
    DhContinuous_.correctBoundaryConditions();

    //- Correct interfacial area density
    iAPtr_->correctField(iA_);

    rhoContinuous_.primitiveFieldRef() =
        fluid1_.rho().primitiveField()*continuity1
    +   fluid2_.rho().primitiveField()*continuity2;
    rhoContinuous_.correctBoundaryConditions();

    nuContinuous_.primitiveFieldRef() = 
        fluid1_.thermo().nu()()*continuity1
    +   fluid2_.thermo().nu()()*continuity2;
    nuContinuous_.correctBoundaryConditions();
    
    PrContinuous_.primitiveFieldRef() = 
        fluid1_.Pr()*continuity1
    +   fluid2_.Pr()*continuity2;
    PrContinuous_.correctBoundaryConditions();
    
    //-
    magUr_ = mag(U1-U2);
    Re_ = max(magUr_*DhDispersed_/nuContinuous_, minRe_);

    //- Correct two-phase quantities. I do this on a by-cell basis to
    //  spare myself the need to do annoying pow operations in cells that
    //  are basically single-phase
    volScalarField& X1(fluid1_.flowQuality());
    volScalarField& X2(fluid2_.flowQuality());
    volScalarField& XLM1(fluid1_.XLM());
    volScalarField& XLM2(fluid2_.XLM());
    const volScalarField& rho1(fluid1_.rho());
    const volScalarField& rho2(fluid2_.rho());
    const volScalarField& mu1(fluid1_.thermo().mu());
    const volScalarField& mu2(fluid2_.thermo().mu());
    const volScalarField& magU1(fluid1_.magU());
    const volScalarField& magU2(fluid2_.magU());
    forAll(mesh_.cells(), i)
    {
        scalar& X1i(X1[i]);
        scalar& X2i(X2[i]);
        scalar& XLM1i(XLM1[i]);
        scalar& XLM2i(XLM2[i]);

        scalar mDot1i(fluid1_[i]*rho1[i]*magU1[i]);
        scalar mDot2i(fluid2_[i]*rho2[i]*magU2[i]);
        scalar mDotSum(max(mDot1i+mDot2i, 1e-69));
        X1i = mDot1i/mDotSum;
        X2i = mDot2i/mDotSum;
        if (X1i < 1e-6)
        {
            XLM1i = fluid1_.minXLM();
        }
        else if (X2i < 1e-6)
        {
            XLM1i = fluid1_.maxXLM();
        }
        else
        {
            XLM1i = 
                min
                (
                    max
                    (
                        pow(mu1[i]/mu2[i], 0.1)*
                        pow(X1i/X2i, 0.9)*
                        sqrt(rho2[i]/rho1[i]),
                        fluid1_.minXLM()
                    ),
                    fluid1_.maxXLM()
                );
        }
        if (X2i < 1e-6)
        {
            XLM2i = fluid2_.minXLM();
        }
        else if (X1i < 1e-6)
        {
            XLM2i = fluid2_.maxXLM();
        }
        else
        {
            XLM2i = 
                min
                (
                    max
                    (
                        pow(mu2[i]/mu1[i], 0.1)*
                        pow(X2i/X1i, 0.9)*
                        sqrt(rho1[i]/rho2[i]),
                        fluid2_.minXLM()
                    ),
                    fluid2_.maxXLM()
                );
        }
        //- Normalize to account for the fact that the clipping for XLM1 and
        //  XLM2 might differ, so that XML1*XLM2 = 1.0 could be violated. So,
        //  restore it via a normalization
        scalar c(sqrt(1.0/(XLM1i*XLM2i)));
        XLM1i *= c;
        XLM2i *= c;
    }
    X1.correctBoundaryConditions();
    X2.correctBoundaryConditions();
    XLM1.correctBoundaryConditions();
    XLM2.correctBoundaryConditions();

    if (correctFluidDynamics)
    {
        //- Correct drag coefficient
        if (KdPtr_.valid())
        {
            myOps::storePrevIterIfRelax(Kd_);
            KdPtr_->correctField(Kd_);
        }
        //- Limit drag coefficient
        scalar KdFF0
        (
            pimple_.dict().lookupOrDefault<scalar>("minKdFF", 1)
        );
        scalar KdFF1
        (
            pimple_.dict().lookupOrDefault<scalar>("maxKdFF", KdFF0)
        );
        scalar a1
        (
            pimple_.dict().lookupOrDefault<scalar>("maxAlphaKdFF", 0)
        );
        if (KdFF1 != KdFF0 and a1 != 0)
        {
            forAll(mesh_.C(), i)
            {
                scalar& Kdi(Kd_[i]);
                scalar Kdlim(KdFF0);
                scalar a(min(fluid1_.normalized()[i], fluid2_.normalized()[i]));
                if (a < a1)
                {
                    scalar c((a1-a)/a1);
                    Kdlim = c*KdFF1 + (1.0-c)*KdFF0;
                }
                Kdi = max(Kdi, Kdlim);   
            }
        }
        else if (KdFF0 != 0)
        {
            forAll(mesh_.C(), i)
            {
                scalar& Kdi(Kd_[i]);
                Kdi = max(Kdi, KdFF0);
            }
        }
        Kd_.correctBoundaryConditions();
        Kd_.relax();

        //- Correct virtual mass forces
        if (virtualMassPtr_.valid())
        {
            virtualMassPtr_->correct();
        }
    }

    if (correctEnergy)
    {
        //- Correct heat transfer coefficient
        volScalarField& htc1(*htcs_[fluid1_.name()]);
        volScalarField& htc2(*htcs_[fluid2_.name()]);
        if (htc1Ptr_.valid())
        {
            myOps::storePrevIterIfRelax(htc1);
            htc1Ptr_->correctField(htc1);
            htc1.relax();
        }
        if (htc2Ptr_.valid())
        {
            myOps::storePrevIterIfRelax(htc2);
            htc2Ptr_->correctField(htc2);
            htc2.relax();
        }
        
        //- Correct interfacial temperature and/or mass transfer
        fluid1_.correctThermoResidualMarkers();
        fluid2_.correctThermoResidualMarkers();
        if (phaseChangePtr_.valid())
        {
            //- Limits interfacial area, updates interfacial temperature, 
            //  latent heat and mass transfer
            phaseChangePtr_->correct();
        }
        else
        {
            const volScalarField& pos1(fluid1_.aboveThermoResidualAlpha());
            const volScalarField& pos2(fluid2_.aboveThermoResidualAlpha());
            const volScalarField& neg1(fluid1_.belowThermoResidualAlpha());
            const volScalarField& neg2(fluid2_.belowThermoResidualAlpha());
            myOps::storePrevIterIfRelax(iT_);
            iT_ = 
                pos1*pos2*
                (
                    htc1*fluid1_.thermo().T() + htc2*fluid2_.thermo().T()
                )/
                (
                    max
                    (   
                        htc1+htc2, 
                        dimensionedScalar
                        (
                            "", 
                            dimPower/dimArea/dimTemperature, 
                            1e-6
                        )
                    )
                )
            +   pos1*neg2*fluid1_.thermo().T()
            +   neg1*pos2*fluid2_.thermo().T();
            if (!firstTimeStepAndIter_)
                iT_.relax();
            else
                firstTimeStepAndIter_ = false;
        }

        //- Add interfacial and mass transfer enthalpy contributions to
        //  heSources
        const volScalarField& he1(fluid1_.thermo().he());
        const volScalarField& he2(fluid2_.thermo().he());
        const volScalarField& Cp1(fluid1_.Cp());
        const volScalarField& Cp2(fluid2_.Cp());
        *(heSources_[he1.name()]) =  
        -   (
                htc1*iA_*(he1/Cp1 + iT_ - fluid1_.thermo().T())
            -   fvm::Sp(htc1*iA_/Cp1, he1)
            )();
        *(heSources_[he2.name()]) =  
        -   (
                htc2*iA_*(he2/Cp2 + iT_ - fluid2_.thermo().T())
            -   fvm::Sp(htc2*iA_/Cp2, he2)
            )();

        if (phaseChangePtr_.valid())
        {
            *heSources_[he1.name()] += phaseChangePtr_->heSource(he1.name()); 
            *heSources_[he2.name()] += phaseChangePtr_->heSource(he2.name()); 
            
            //- Re-set htc to their values before their were modified by
            //  phaseChangePtr_->correctHeSources() (which caches them in
            //  prevIter)
            htc1 = htc1.prevIter();
            htc2 = htc2.prevIter();
        }
    }
}

// ************************************************************************* //
