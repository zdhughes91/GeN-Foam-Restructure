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

#include "FSPair.H"
#include "myOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(FSPair, 0);

    typedef HashTable
    <
        autoPtr<FSDragFactor>,
        word,
        word::hash
    > KdTable;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FSPair::FSPair
(
    fluid& fluidRef,
    structure& structureRef,
    const dictionary& dict
)
:
    regIOobject
    (
        IOobject
        (
            (
                (fluidRef.mesh().thisDb().lookupClass<fluid>().size() > 1) ? 
                IOobject::groupName(fluidRef.name(), "structure") : 
                "fluid.structure"
            ),
            fluidRef.mesh().time().timeName(),
            fluidRef.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true
        )
    ),
    mesh_(fluidRef.mesh()),
    dict_(dict),
    pimple_(mesh_.solutionDict().subDict("PIMPLE")),
    fluid_(fluidRef),
    structure_(structureRef),
    minRe_
    (
        dimensionedScalar::lookupOrDefault
        (
            "residualFluidStructureRe",
            dict,
            dimless,
            10
        )
    ),
    Re_
    (
        IOobject
        (
            (
                (mesh_.thisDb().lookupClass<fluid>().size() > 1) ?
                IOobject::groupName("Re", this->name()) :
                "Re"
            ),
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        minRe_,
        zeroGradientFvPatchScalarField::typeName
    ),
    lRe_
    (
        IOobject
        (
            (
                (mesh_.thisDb().lookupClass<fluid>().size() > 1) ?
                IOobject::groupName("lRe", this->name()) :
                "lRe"
            ),
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("", dimless, vector::one),
        zeroGradientFvPatchVectorField::typeName
    ),
    Kd_
    (
        IOobject
        (
            (
                (mesh_.thisDb().lookupClass<fluid>().size() > 1) ?
                IOobject::groupName("Kd", this->name()) :
                "Kd"
            ),
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedTensor("", dimDensity/dimTime, tensor::zero),
        zeroGradientFvPatchTensorField::typeName
    ),
    htc_
    (
        IOobject
        (
            (
                (mesh_.thisDb().lookupClass<fluid>().size() > 1) ?
                IOobject::groupName("htc", this->name()) :
                "htc"
            ),
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimPower/dimArea/dimTemperature, 0),
        zeroGradientFvPatchScalarField::typeName
    )
{
    Info << endl;
    
    bool onePhase(mesh_.thisDb().lookupClass<fluid>().size() == 1);
    
    const dictionary& physicsModelsDict(dict_.subDict("physicsModels"));
    
    //- Construct drag and heat transfer models
    const dictionary* dragModelsDictPtr(nullptr);
    const dictionary* heatTransferModelsDictPtr(nullptr);
    if (onePhase)
    {
        dragModelsDictPtr = 
            &(physicsModelsDict.subDict("dragModels"));
        heatTransferModelsDictPtr = 
            &(physicsModelsDict.subDict("heatTransferModels"));
    }
    else
    {
        if 
        (
            foundUnorderedPairSubDict
            (
                physicsModelsDict.subDict("dragModels")
            )
        )
        {
            dragModelsDictPtr = 
                &(
                    getUnorderedPairSubDict
                    (
                        physicsModelsDict.subDict("dragModels")
                    )
                );
        }
        if 
        (
            foundUnorderedPairSubDict
            (
                physicsModelsDict.subDict("heatTransferModels")
            )
        )
        {
            heatTransferModelsDictPtr = 
                &(
                    getUnorderedPairSubDict
                    (
                        physicsModelsDict.subDict("heatTransferModels")
                    )
                );
        }
    }
    if (dragModelsDictPtr != nullptr)
    {
        const dictionary& dragModelsDict(*dragModelsDictPtr);
        wordList dragModelKeys(dragModelsDict.toc());
        forAll(dragModelKeys, i)
        {
            word key(dragModelKeys[i]);
            KdPtrs_.set
            (
                key,
                autoPtr<FSDragFactor>
                (
                    new FSDragFactor
                    (
                        *this,
                        dragModelsDict.subDict(key)
                    )
                )
            );
        }
    }
    if (heatTransferModelsDictPtr != nullptr)
    {
        const dictionary& heatTransferModelsDict(*heatTransferModelsDictPtr);
        wordList heatTransferModelKeys(heatTransferModelsDict.toc());
        forAll(heatTransferModelKeys, i)
        {
            word key(heatTransferModelKeys[i]);
            htcPtrs_.set
            (
                key,
                FSHeatTransferCoefficientModel::New
                (
                    *this,
                    heatTransferModelsDict.subDict(key),
                    mesh_
                )
            );
        }
    }
    
    //- Contact fraction models
    if (!onePhase)
    {
        const dictionary& pairGeometryModels
        (
            physicsModelsDict.subDict("pairGeometryModels")
        );
        if 
        (
            foundUnorderedPairSubDict(pairGeometryModels)
        )
        {
            contactPartition_.reset
            (
                contactPartitionModel::New
                (
                    *this,
                    getUnorderedPairSubDict(pairGeometryModels).subDict
                    (
                        "contactPartitionModel"
                    ),
                    mesh_
                )
            );
        }
        else
        {
            dictionary dict;
            dict.set("type", "complementary");
            contactPartition_.reset
            (
                contactPartitionModel::New
                (
                    *this,
                    dict,
                    mesh_
                )
            );
        }

        //- 
        fPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("f", this->name()),
                    mesh_.time().timeName(),
                    mesh_
                ),
                mesh_,
                dimensionedScalar("", dimless, 0.0),
                zeroGradientFvPatchScalarField::typeName
            )
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::FSPair::correct
(
    const bool& correctFluidDynamics, 
    const bool& correctEnergy
)
{
    //- Init refs
    const volVectorField& U(fluid_.U());

    //- Previously fluid_.Dh(). The difference is that fluid_.Dh() will be
    //  different than the structure Dh in regions where fluid_ is dispersed.
    //  Nonetheless, the Reynolds for the fluid-structure pair when the fluid
    //  is dispersed is fairly meaningless, so I changed it to utilize the
    //  structure_.Dh() which does not change in time.
    const volScalarField& Dh(structure_.Dh());
    const volVectorField& lDh(structure_.lDh());
    
    //- Correct Reynolds
    volScalarField nu(fluid_.thermo().nu());
    Re_ = fluid_.normalized()*mag(U)*Dh/nu;
    Re_ = max(Re_, minRe_);

    //- Calculate anisotropic reynolds
    /*  This was the original formulation. Below is the much more efficient one
        which is cell-by-cell, as I do not really care about BC values
    //- Rotate U in the local frame, make an outer product with lDh so that
    //  we get a matrix whose diagonal elements are the components of lRe.
    //  Then directly set these diagonal elemets (I guess this is the most
    //  efficient approch)
    volTensorField tlRet = 
        fluid_.normalized()*(structure_.Rg2l()&U)*lDh/fluid_.thermo().nu()();
    lRe_.replace(0, max(tlRet.component(0), minRe_));
    lRe_.replace(1, max(tlRet.component(4), minRe_));
    lRe_.replace(2, max(tlRet.component(8), minRe_));
    */
    scalar minRe(minRe_.value());
    if (structure_.hasLocalReferenceFrame())
    {
        //- This is to rotate the fluid velocity from the global frame to the
        //  local one
        const volTensorField R(structure_.Rg2l());
        forAll(mesh_.cells(), i)
        {
            scalar alphaByNu(fluid_.normalized()[i]/nu[i]);
            vector Ui(R[i]&U[i]);
            const vector& lDhi(lDh[i]);
            lRe_[i][0] = max(alphaByNu*Ui[0]*lDhi[0], minRe);
            lRe_[i][1] = max(alphaByNu*Ui[1]*lDhi[1], minRe);
            lRe_[i][2] = max(alphaByNu*Ui[2]*lDhi[2], minRe);
        }
    }
    else
    {
        forAll(mesh_.cells(), i)
        {
            scalar alphaByNu(fluid_.normalized()[i]/nu[i]);
            const vector& Ui(U[i]);
            const vector& lDhi(lDh[i]);
            lRe_[i][0] = max(alphaByNu*Ui[0]*lDhi[0], minRe);
            lRe_[i][1] = max(alphaByNu*Ui[1]*lDhi[1], minRe);
            lRe_[i][2] = max(alphaByNu*Ui[2]*lDhi[2], minRe);
        }
    }
    Re_.correctBoundaryConditions();
    lRe_.correctBoundaryConditions();

    if (correctFluidDynamics)
    {
        myOps::storePrevIterIfRelax(Kd_);
        //- Correct drag factor
        forAllIter
        (
            KdTable,
            KdPtrs_,
            iter
        )
        {
            FSDragFactor& dragFactor(iter()());
            dragFactor.correctField(Kd_);
        }

        //- Limit drag factor
        scalar KdFS
        (
            pimple_.lookupOrDefault<scalar>("minKdFS", 0)
        );
        if (KdFS != 0)
        {
            forAll(mesh_.C(), i)
            {
                tensor& Kdi(Kd_[i]);
                Kdi[0] = max(Kdi[0], KdFS);
                Kdi[4] = max(Kdi[0], KdFS);
                Kdi[8] = max(Kdi[0], KdFS);
            }
        }
        Kd_.correctBoundaryConditions();
        Kd_.relax();
        
        //- Rotate drag
        structure_.localToGlobalRotateField(Kd_);
    }

    //- Correct contact partition fraction, if valid
    if (fPtr_.valid())
        contactPartition_->correctField(fPtr_());
    
    if (correctEnergy)
    {
        myOps::storePrevIterIfRelax(htc_);
        //- Correct heat transfer coefficient
        forAllIter
        (
            htcTable,
            htcPtrs_,
            iter
        )
        {
            FSHeatTransferCoefficientModel& heatTransferCoefficient(iter()());
            heatTransferCoefficient.correctField(htc_);
        }
        htc_.correctBoundaryConditions();

        //- Correct htc to account for contact partition, if valid
        if (fPtr_.valid())
            htc_ *= fPtr_();
        htc_.relax();
    }
}


// ************************************************************************* //
