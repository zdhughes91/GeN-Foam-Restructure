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

#include "fluid.H"
#include "fvCFD.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvcFlux.H"
#include "fvcSnGrad.H"
#include "fvcCurl.H"
#include "surfaceInterpolate.H"
#include "zeroGradientFvPatchFields.H"
#include "fixedValueFvsPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "slipFvPatchFields.H"
#include "partialSlipFvPatchFields.H"
#include "myOps.H"

#include "fluidDiameterModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::fluid::stateOfMatter
>
Foam::fluid::stateOfMatterNames_
(
    {
        { 
            stateOfMatter::undetermined, 
            "undetermined" 
        },
        { 
            stateOfMatter::liquid, 
            "liquid" 
        },
        { 
            stateOfMatter::gas, 
            "gas" 
        }
    }
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluid::fluid
(
    const dictionary& dict,
    const fvMesh& mesh,
    const word& phaseName,
    volScalarField& powerDensityNeutronicsToLiquid,
    bool readIfPresentAndWrite
)
:
    phaseBase
    ( 
        dict,
        mesh,
        phaseName,
        calculatedFvPatchScalarField::typeName,
        readIfPresentAndWrite,
        readIfPresentAndWrite
    ),
    stateOfMatter_
    (
        stateOfMatterNames_.get
        (
            dict_.lookupOrDefault<word>
            (
                "stateOfMatter", 
                "undetermined"
            )
        )
    ),
    U_
    (
        IOobject
        (
            IOobject::groupName("U", this->name()),
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("", dimVelocity, vector(0,0,0)),
        slipFvPatchVectorField::typeName
    ),
    magU_
    (
        IOobject
        (
            IOobject::groupName("magU", this->name()),
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimVelocity, scalar(0))
    ),
    alphaPhi_
    (
        IOobject
        (
            IOobject::groupName("alphaPhi", this->name()),
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            (
                (
                    mesh.time().controlDict().lookupOrDefault<bool>
                    (
                        "writeRestartFields", 
                        true
                    )
                ) ?
                IOobject::AUTO_WRITE :
                IOobject::NO_WRITE
            )
        ),
        mesh,
        dimensionedScalar("", dimVol/dimTime, 0)
    ),
    alphaRhoPhi_
    (
        IOobject
        (
            IOobject::groupName("alphaRhoPhi", this->name()),
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            (
                (
                    mesh.time().controlDict().lookupOrDefault<bool>
                    (
                        "writeRestartFields", 
                        true
                    )
                ) ?
                IOobject::AUTO_WRITE :
                IOobject::NO_WRITE
            )
        ),
        mesh,
        dimensionedScalar("", dimMass/dimTime, 0)
    ),
    alphaRhoMagU_
    (
        IOobject
        (
            IOobject::groupName("alphaRhoMagU", this->name()),
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            (
                (
                    mesh.time().controlDict().lookupOrDefault<bool>
                    (
                        "writeRestartFields", 
                        true
                    )
                ) ?
                IOobject::AUTO_WRITE :
                IOobject::NO_WRITE
            )
        ),
        mesh,
        dimensionedScalar("", dimMass/dimArea/dimTime, 0)
    ),
    dgdt_
    (
        IOobject
        (
            IOobject::groupName("dgdt", this->name()),
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            (
                (
                    mesh.time().controlDict().lookupOrDefault<bool>
                    (
                        "writeRestartFields", 
                        true
                    )
                ) ?
                IOobject::AUTO_WRITE :
                IOobject::NO_WRITE
            )
        ),
        mesh,
        dimensionedScalar("", dimless/dimTime, 0)
    ),
    minXLM_
    (
        dict.getOrDefault<scalar>
        (
            "minXLM",
            0.0001
        )
    ),
    maxXLM_
    (
        dict.getOrDefault<scalar>
        (
            "maxXLM",
            1.0/minXLM_
        )
    ),
    normalized_
    (
        IOobject
        (
            IOobject::groupName("normalized.alpha", this->name()),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimless, 0)
    ),
    kappa_
    (
        IOobject
        (
            IOobject::groupName("kappa", this->name()),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimPower/dimLength/dimTemperature, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    Cp_
    (
        IOobject
        (
            IOobject::groupName("Cp", this->name()),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimEnergy/dimMass/dimTemperature, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    mu_
    (
        IOobject
        (
            IOobject::groupName("mu", this->name()),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimPressure*dimTime, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    Pr_
    (
        IOobject
        (
            IOobject::groupName("Pr", this->name()),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimless, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    powerDensityNeutronicsToLiquid_(powerDensityNeutronicsToLiquid),
    contErr_
    (
        IOobject
        (
            IOobject::groupName("contErr", this->name()),
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            (
                (
                    mesh.time().controlDict().lookupOrDefault<bool>
                    (
                        "writeContinuityErrors", 
                        false
                    )
                ) ?
                IOobject::AUTO_WRITE :
                IOobject::NO_WRITE
            )
        ),
        mesh,
        dimensionedScalar("", dimDensity/dimTime, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    cumulContErr_(0),
    Dh_
    (
        IOobject
        (
            IOobject::groupName("Dh", this->name()),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimLength, SMALL),
        zeroGradientFvPatchScalarField::typeName
    ),
    //dispersion_(mesh.C().size(), scalar(0.0)),
    thermoResidualAlpha_
    (
        dimensionedScalar::lookupOrDefault
        (
            "thermoResidualAlpha",
            dict_,
            dimless, 
            0.0
        )
    ),
    Boussinesq_(false),
    rho0Ptr_(nullptr)
{
    if (phaseName != "")
    {
        Info << endl << "Constructing fluid: " << this->name() << endl;
    }
    else
    {
        Info << endl << "Constructing fluid" << endl;
    }

    //- Construct thermodynamics package
    thermo_.reset(rhoThermo::New(mesh, this->name()));

    mesh.setFluxRequired(this->name());

    //- Set initial cellZone powerDensity, 
    //  if no field already available in time folder
    if (dict_.found("initialPowerDensity"))
    {
        if(!powerDensityNeutronicsToLiquid.typeHeaderOk<volScalarField>(true))
        {
            //volScalarField& powerDensity(powerDensityPtr_());
            volScalarField& powerDensity(powerDensityNeutronicsToLiquid);
            const dictionary& powerDensity0s
            (
                dict_.subDict("initialPowerDensity")
            );

            forAllConstIter
            (
                dictionary,
                powerDensity0s,
                powerDensity0Iter
            )
            {
                DynamicList<label> cells(0);
                word zoneName(powerDensity0Iter->keyword());
                scalar powerDensity0(powerDensity0s.get<scalar>(zoneName));
                forAllConstIter
                (
                    DynamicList<label>,
                    mesh.cellZones()[zoneName],
                    cIter
                )
                {
                    cells.append(*cIter);
                }
                forAll(cells, i)
                {
                    powerDensity[cells[i]] = powerDensity0;
                }
            }

            powerDensity.correctBoundaryConditions();
        }
    }

    //- Set initial cellZone phase fractions, if present
    if (dict_.found("initialAlpha"))
    {
        const dictionary& alpha0s
        (
            dict_.subDict("initialAlpha")
        );

        forAllConstIter
        (
            dictionary,
            alpha0s,
            alpha0Iter
        )
        {
            DynamicList<label> cells(0);
            word zoneName(alpha0Iter->keyword());
            scalar alpha0(alpha0s.get<scalar>(zoneName));
            forAllConstIter
            (
                DynamicList<label>,
                mesh.cellZones()[zoneName],
                cIter
            )
            {
                cells.append(*cIter);
            }
            forAll(cells, i)
            {
                (*this)[cells[i]] = alpha0;
            }
        }

        this->correctBoundaryConditions();
    }

    thermo_->validate(phaseName, "h", "e");

    //- Set Boussinesq_ flag by reading thermophysicalProperties data and init
    //  rho0Ptr if using the Boussinesq approx
    List<char> delimiters(3);
    delimiters[0] = '<';
    delimiters[1] = '>';
    delimiters[2] = ',';
    for 
    (
        word const &w : 
        myOps::split<word>
        (
            thermo_->thermoName(), delimiters
        )
    )
    {
        if (w == "Boussinesq")
        {
            Boussinesq_ = true;
            break;
        }
    }
    if (Boussinesq_)
    {
        rho0Ptr_ = 
            new volScalarField
            (
                IOobject
                (
                    "",
                    mesh_.time().timeName(),
                    mesh_
                ),
                mesh_,
                dimensionedScalar
                (
                    "rho0", 
                    dimDensity, 
                    thermo_->subDict("mixture").subDict("equationOfState")
                ),
                zeroGradientFvPatchScalarField::typeName
            );
    }

    //- The rest of the constructor is only for correctly setting the
    //  boundary conditions of phi
    const word phiName = IOobject::groupName("phi", this->name());

    IOobject phiHeader
    (
        phiName,
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ
    );

    if (phiHeader.typeHeaderOk<surfaceScalarField>(true))
    {
        Info<< "Reading face flux field " << phiName << endl;

        phiPtr_.reset
        (
            new surfaceScalarField
            (
                IOobject
                (
                    phiName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            )
        );
    }
    else
    {
        Info<< "Calculating face flux field " << phiName << endl;

        wordList phiTypes
        (
            U_.boundaryField().size(),
            calculatedFvPatchScalarField::typeName
        );

        forAll(U_.boundaryField(), i)
        {
            if
            (
                isA<fixedValueFvPatchVectorField>(U_.boundaryField()[i])
            ||  isA<slipFvPatchVectorField>(U_.boundaryField()[i])
            ||  isA<partialSlipFvPatchVectorField>(U_.boundaryField()[i])
            )
            {
                phiTypes[i] = fixedValueFvsPatchScalarField::typeName;
            }
        }

        phiPtr_.reset
        (
            new surfaceScalarField
            (
                IOobject
                (
                    phiName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    (
                        (
                            mesh.time().controlDict().lookupOrDefault<bool>
                            (
                                "writeRestartFields", 
                                true
                            )
                        ) ?
                        IOobject::AUTO_WRITE :
                        IOobject::NO_WRITE
                    )
                ),
                fvc::flux(U_),
                phiTypes
            )
        );
    }
    
    //- What about alphaPhi, alphaRhoPhi? Well, these depend on the 
    //  phase fraction (unlike phi) but at this step there have been no
    //  phase fractions normalizations (which can only be done by the main
    //  solver). Thus, rather than tentatively set the fields twice (here,
    //  maybe wrong, and then in the main to correct for potential phase
    //  fraction normalization issues), this needs to be handled by the
    //  main solver via the initAlphaPhis function. Needless to say, if
    //  alphaPhi and alphaRhoPhi are found on disk, those are read and that's
    //  the end of it
    
    //- Init placeholder fields
    kappa_ = thermo_->kappa();
    Cp_ = thermo_->Cp();
    mu_ = thermo_->mu();
    magU_ = mag(U_);
    Pr_ = thermo_->Cp()*thermo_->mu()/thermo_->kappa();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluid::~fluid()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fluid::initTwoPhaseFields()
{
    phaseChangeSignPtr_.reset
    (
        new scalarField(mesh_.C().size(), int(0))
    );
    
    flowQualityPtr_.reset
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("flowQuality", this->name()),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                (
                    (
                        mesh_.time().controlDict().lookupOrDefault<bool>
                        (
                            "writeAllFields", 
                            false
                        )
                    ) ?
                    IOobject::AUTO_WRITE :
                    IOobject::NO_WRITE
                )
            ),
            mesh_,
            dimensionedScalar("", dimless, 0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    
    XLMPtr_.reset
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("XLM", this->name()),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                (
                    (
                        mesh_.time().controlDict().lookupOrDefault<bool>
                        (
                            "writeAllFields", 
                            false
                        )
                    ) ?
                    IOobject::AUTO_WRITE :
                    IOobject::NO_WRITE
                )
            ),
            mesh_,
            dimensionedScalar("", dimless, 10000),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    dispersionPtr_.reset
    (
        new scalarField(mesh_.C().size(), scalar(0.0))
    );
}

void Foam::fluid::initAlphaPhis()
{
    IOobject alphaPhiHeader
    (
        alphaPhi_.name(),
        mesh_.time().timeName(),
        mesh_,
        IOobject::NO_READ
    );
    IOobject alphaRhoPhiHeader
    (
        alphaRhoPhi_.name(),
        mesh_.time().timeName(),
        mesh_,
        IOobject::NO_READ
    );
    IOobject alphaRhoMagUHeader
    (
        alphaRhoMagU_.name(),
        mesh_.time().timeName(),
        mesh_,
        IOobject::NO_READ
    );

    if (!alphaPhiHeader.typeHeaderOk<surfaceScalarField>(true))
    {
        alphaPhi_ = fvc::interpolate(*this)*phiPtr_();
    }
    if (!alphaRhoPhiHeader.typeHeaderOk<surfaceScalarField>(true))
    {
        alphaRhoPhi_ = fvc::interpolate(this->rho())*alphaPhi_;
    }
    if (!alphaRhoMagUHeader.typeHeaderOk<volScalarField>(true))
    {
        alphaRhoMagU_ = (*this)*this->rho()*magU_;
    }
}

void Foam::fluid::constructDiameterModel()
{
    diameterPtr_.reset
    (
        fluidDiameterModel::New
        (
            *this,
            dict_.subDict("dispersedDiameterModel"),
            mesh_
        )
    );
}

void Foam::fluid::constructTurbulenceModel()
{
    turbulence_ =
        phaseCompressibleTurbulenceModel::New
        (
            *this,
            this->rho(),
            U_,
            alphaRhoPhi_,
            phi(),
            thermo_
        ); 
}

void Foam::fluid::correctAlphaRhoMagU()
{
    alphaRhoMagU_ = (*this)*this->rho()*magU_;
}

void Foam::fluid::correctDiameter()
{
    if (diameterPtr_.valid())
    {
        diameterPtr_->correctField(Dh_);
    }
    
}

void Foam::fluid::correctThermoResidualMarkers()
{
    if (aboveThermoResidualAlphaPtr_.valid())
    {
        aboveThermoResidualAlphaPtr_() = 
            pos(normalized_-thermoResidualAlpha_);
    }
    else
    {
        aboveThermoResidualAlphaPtr_.reset
        (
            new volScalarField(pos(normalized_-thermoResidualAlpha_))
        );
    }
    if (belowThermoResidualAlphaPtr_.valid())
    {
        belowThermoResidualAlphaPtr_() = 
            neg0(normalized_-thermoResidualAlpha_);
    }
    else
    {
        belowThermoResidualAlphaPtr_.reset
        (
            new volScalarField(neg0(normalized_-thermoResidualAlpha_))
        );
    }
}

volScalarField& Foam::fluid::rho(bool isVariableIfBoussinesq)
{
    if (Boussinesq_ and !isVariableIfBoussinesq)
        return *rho0Ptr_;
    else
        return thermo_->rho();
}

const volScalarField& Foam::fluid::rho(bool isVariableIfBoussinesq) const
{
    if (Boussinesq_ and !isVariableIfBoussinesq)
        return *rho0Ptr_;
    else
        return thermo_->rho();
}


// ************************************************************************* //
