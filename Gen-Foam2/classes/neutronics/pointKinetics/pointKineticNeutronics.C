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

#include "pointKineticNeutronics.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "coordinateSystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pointKineticNeutronics, 0);

    addToRunTimeSelectionTable
    (
        neutronics,
        pointKineticNeutronics,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointKineticNeutronics::pointKineticNeutronics
(
    fvMesh& mesh
)
:
    neutronics(mesh),
    nuclearData_
    (
        IOobject
        (
            "nuclearData",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    power_(reactorState_.get<scalar>("pTarget")),
    fissionPower_(power_),
    decayPower_(0.0),
    decayPowerPtr_(nullptr),
    decayPowerStartTime_(0.0),
    externalReactivityTimeProfile_
    (
        nuclearData_,
        "externalReactivityTimeProfile"
    ),
    precEquilibriumReactivity_(0.0),
    liquidFuelBeta_(0.0),
    totalReactivity_(0.0),
    promptGenerationTime_(nuclearData_.get<scalar>("promptGenerationTime")),
    beta_(0.0),
    betas_
    (
        nuclearData_.get<scalarList>("Beta") // ...
    ),
    lambdas_
    (
        nuclearData_.get<scalarList>("lambda") // ... ...
    ),
    precursorPowers_
    (
        reactorState_.lookupOrDefault<scalarList>
        (
            "precursorPowers",
            scalarList(betas_.size(), 0.0)
        )
    ),
    delayedGroups_(betas_.size()),
    timeIndex_(mesh.time().timeIndex()),
    powerOld_(power_),
    fissionPowerOld_(fissionPower_),
    precursorPowersOld_(precursorPowers_),
    defaultPrec_
    (
        IOobject
        (
            "defaultPrec",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimless/dimVol, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    coeffTFuel_
    (
        nuclearData_.get<scalar>("feedbackCoeffTFuel")
    ),
    coeffTClad_
    (
        nuclearData_.get<scalar>("feedbackCoeffTClad")
    ),
    coeffTCool_
    (
        nuclearData_.get<scalar>("feedbackCoeffTCool")
    ),
    coeffRhoCool_
    (
        nuclearData_.get<scalar>("feedbackCoeffRhoCool")
    ),
    coeffTStruct_
    (
        nuclearData_.get<scalar>("feedbackCoeffTStruct")
    ),
    coeffDrivelineExp_
    (
        nuclearData_.get<scalar>("absoluteDrivelineExpansionCoeff")
    ),
    boronReactivityTimeProfile_
    (
        nuclearData_,
        "boronReactivityTimeProfile"
    ),
    TFuel_
    (
        IOobject
        (
            "pointKinetics.TFuel",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimTemperature, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    TClad_
    (
        IOobject
        (
            "TClad",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimTemperature, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    TCool_
    (
        IOobject
        (
            "TCool",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimTemperature, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    rhoCool_
    (
        IOobject
        (
            "rhoCool",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimDensity, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    TStruct_
    (
        IOobject
        (
            "TStruct",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimTemperature, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    UPtr_(nullptr),
    alphaPtr_(nullptr),
    alphatPtr_(nullptr),
    muPtr_(nullptr),
    phiPtr_(nullptr),
    diffCoeffPrecPtr_(nullptr),
    Dalbedo_
    (
        IOobject
        (
        "Dalbedo",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimLength, 1.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    fluxStarAlbedo_
        (
        IOobject
        (
            "fluxStarAlbedo",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimless/dimArea/dimTime, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    oneGroupFlux_
    (
        IOobject
        (
            "oneGroupFlux",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimless/dimArea/dimTime, 1),
        zeroGradientFvPatchScalarField::typeName
    ),
    initOneGroupFlux_
    (
        IOobject
        (
            "initOneGroupFlux",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        oneGroupFlux_
    ),
    initOneGroupFluxN_
    (
        IOobject
        (
            "initOneGroupFluxN",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        initOneGroupFlux_/fvc::domainIntegrate(initOneGroupFlux_)
    ),
    domainIntegratedInitOneGroupFluxN_(fvc::domainIntegrate(sqr(initOneGroupFluxN_)).value()),
    energyGroups_(0),
    fluxes_(0),
    precursors_(0),
    precPK_(delayedGroups_),
    precPKStar_(delayedGroups_),
    fastNeutrons_
    (
        nuclearData_.lookupOrDefault<bool>("fastNeutrons", false)
    ),
    coeffFastDoppler_
    (
        (fastNeutrons_) ?
        nuclearData_.get<scalar>("feedbackCoeffFastDoppler") :
        0.0
    ),
    fuelFeedbackCellField_
    (
        IOobject
        (
            "pointKinetics.fuelFeedbackCellField",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimless, 0)
    ),
    coolFeedbackCellField_
    (
        IOobject
        (
            "pointKinetics.coolFeedbackCellField",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimless, 0)
    ),
    structFeedbackCellField_
    (
        IOobject
        (
            "pointKinetics.structFeedbackCellField",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimless, 0)
    ),
    drivelineFeedbackCellField_
    (
        IOobject
        (
            "pointKinetics.drivelineFeedbackCellField",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimless, 0)
    ),
    controlRodReactivityMap_
    (
        nuclearData_.lookupOrDefault<List<Pair<scalar>>>
        (
            "controlRodReactivityMap",
            List<Pair<scalar>>()
        )
    ),
    ScNo_
    (
        nuclearData_.lookupOrDefault<scalar>("ScNo",1.0)
    ),
    initPrecursorsLiquidFuel_
    (
        nuclearData_.lookupOrDefault<bool>("initPrecursorsLiquidFuel",false)
    ),
    GEMReactivityMap_
    (
        nuclearData_.lookupOrDefault<List<Pair<scalar>>>
        (
            "GEMReactivityMap",
            List<Pair<scalar>>()
        )
    ),
    intPhiRef_(0),
    phiFaces_(0),
    phiMagSf_(0)
{
    //- Cannot work in eigenvalue mode for obvious reasons, it makes no sense
    if (eigenvalueNeutronics_)
    {
        FatalErrorInFunction
            << "pointKinetics model incompatible with eigenvalueNeutronics"
            << exit(FatalError);
    }

    //- Check if decay power provided
    word decayPowerDictName("decayPowerTimeProfile");
    if (reactorState_.found(decayPowerDictName))
    {
        const dictionary& decayPowerDict
        (
            reactorState_.subDict(decayPowerDictName)
        );
        word type
        (
            decayPowerDict.get<word>("type")
        );
        decayPowerPtr_.reset
        (
            Function1<scalar>::New
            (
                type,
                decayPowerDict,
                type
            )
        );
        decayPowerStartTime_ =
            decayPowerDict.lookupOrDefault<scalar>("startTime", 0.0);
        const scalar& t(mesh_.time().timeOutputValue());
        decayPower_ = decayPowerPtr_->value(t-decayPowerStartTime_);
        fissionPower_ = power_ - decayPower_;
        fissionPowerOld_ = fissionPower_;

    }

    if (liquidFuel_)
    {
        UPtr_.reset
        (
            new volVectorField
            (
                IOobject
                (
                    "U",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedVector("", dimVelocity, vector::zero),
                zeroGradientFvPatchVectorField::typeName
            )
        );
        alphaPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "alpha",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("", dimless, 1.0),
                zeroGradientFvPatchScalarField::typeName
            )
        );
        alphatPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "alphat",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("", dimMass/dimLength/dimTime, 0.0),
                zeroGradientFvPatchScalarField::typeName
            )
        );
        muPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "mu",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("", dimMass/dimLength/dimTime, 0.0),
                zeroGradientFvPatchScalarField::typeName
            )
        );
        phiPtr_.reset
        (
            new surfaceScalarField
            (
                IOobject
                (
                    "phi",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fvc::flux(UPtr_())
            )
        );
        diffCoeffPrecPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "diffCoeffPrec",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("", dimArea/dimTime, 0.0),
                zeroGradientFvPatchScalarField::typeName
            )
        );
    }

    //- Check that provided power is > 0
    if (power_ <= 0)
    {
        FatalErrorInFunction
            << "Set pTarget > 0 !"
            << exit(FatalError);
    }

    //- Check that provided promptGenerationTime is > 0
    if (promptGenerationTime_ <= 0)
    {
        FatalErrorInFunction
            << "Set promptGenerationTime > 0 !"
            << exit(FatalError);
    }

    //- Create fluxes_ PtrList. If pre-existing flux fields exist in the start
    //  time folder, read those and determine the number of energy groups. If
    //  they do not exist, read the number of energy groups from nuclearData
    //  and create the corresponding number of flux fields, scaled so to add
    //  up to oneGroupFlux (if it does not exist, defaults to 1/m2/s). If
    //  the energyGroups keyword is not found, deafults to only one energy
    //  group. Again, most of this has nothing to do with the pointKinetics
    //  model iteself (there is no energy dependence), but it is used to
    //  scale all the fluxes according to pointKinetic results
    while (true)
    {
        word fluxName = "flux"+Foam::name(energyGroups_);

        IOobject fluxHeader
        (
            fluxName,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ
        );

        if (fluxHeader.typeHeaderOk<volScalarField>(true))
        {
            fluxes_.append
            (
                new volScalarField
                (
                    IOobject
                    (
                        fluxName,
                        mesh.time().timeName(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh
                )
            );
            energyGroups_++;
        }
        else
        {
            break;
        }
    }
    if (fluxes_.size() == 0)
    {
        //- One group flux is init from this dict ONLY IF no existing flux
        //  files are already present (otherwise it is just reconstructed
        //  from the sum of those)
        if (nuclearData_.found("initialOneGroupFluxByZone"))
        {
            const dictionary& initialOneGroupFluxes
            (
                nuclearData_.subDict("initialOneGroupFluxByZone")
            );
            oneGroupFlux_ *= 0.0;
            forAllConstIter
            (
                dictionary,
                initialOneGroupFluxes,
                iter
            )
            {
                DynamicList<label> cells(0);
                word zoneName(iter->keyword());
                scalar initialOneGroupFlux
                (
                    initialOneGroupFluxes.get<scalar>(zoneName)
                );
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
                    oneGroupFlux_[cells[i]] = initialOneGroupFlux;
                }
            }
            oneGroupFlux_.correctBoundaryConditions();
            setInitOneGroupFlux();
        }

        energyGroups_ =
            nuclearData_.lookupOrDefault<scalar>("energyGroups", 1);

        for (int i = 0; i < energyGroups_; i++)
        {
            fluxes_.append
            (
                new volScalarField
                (
                    IOobject
                    (
                        "flux"+Foam::name(i),
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        (
                            energyGroups_ == 1 ?
                            IOobject::NO_WRITE :
                            IOobject::AUTO_WRITE
                        )
                    ),
                    (1.0/energyGroups_)*oneGroupFlux_
                )
            );
            fluxes_[i].correctBoundaryConditions();
        }
    }
    else
    {
        oneGroupFlux_ *= 0.0;
        forAllConstIter
        (
            PtrList<volScalarField>,
            fluxes_,
            iter
        )
        {
            const volScalarField& fluxi(iter());
            oneGroupFlux_ += fluxi;
        }
        oneGroupFlux_.correctBoundaryConditions();
    }
    //- Update the initOneGroupFlux related quantities after possible
    //  changes in the oneGroupFlux
    setInitOneGroupFlux();

    //- Read real precursors if present, they only get re-scaled by
    //  pointKinetic results
    label i = 0;
    while (true)
    {
        word precName = "prec"+Foam::name(i);

        IOobject precHeader
        (
            precName,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ
        );

        if (precHeader.typeHeaderOk<volScalarField>(true))
        {
            precursors_.append
            (
                new volScalarField
                (
                    IOobject
                    (
                        precName,
                        mesh.time().timeName(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh
                )
            );
            i++;
        }
        else
        {
            break;
        }
    }

    forAll(precPK_,precI)
    {
        if(liquidFuel_)
        {
            precPK_.set
            (
                precI,
                new volScalarField
                (
                    IOobject
                    (
                        "precPK"+Foam::name(precI),
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::READ_IF_PRESENT,
                        IOobject::AUTO_WRITE
                    ),
                    //chenge the precursors units to power for PK calculations
                    defaultPrec_ //* dimensionedScalar("", dimPower, 1.0)
                )
            );
            precPK_[precI].dimensions().reset(defaultPrec_.dimensions()*dimPower);
            precPKStar_.set
            (
                precI,
                new volScalarField
                (
                    IOobject
                    (
                        "precPKStar"+Foam::name(precI),
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::READ_IF_PRESENT,
                        IOobject::AUTO_WRITE
                    ),
                    //chenge the precursors units to power for PK calculations
                    defaultPrec_ //* dimensionedScalar("", dimPower, 1.0)
                )
            );
            precPKStar_[precI].dimensions().reset(defaultPrec_.dimensions()*dimPower);
        }
        else
        {
            precPK_.set(precI, nullptr);
            precPKStar_.set(precI, nullptr);
        }

    }

    //- Set precursorPowers so that, if they are not found in reactorState, the
    //  system starts from a steady-state.
    if (!reactorState_.found("precursorPowers"))
    {
        for (int i = 0; i < delayedGroups_; i++)
        {
            precursorPowers_[i] =
                (betas_[i]*fissionPower_)/(lambdas_[i]*promptGenerationTime_);
        }
    }

    //- Compute total effective delayed neutron fraction
    forAll(betas_, i)
    {
        beta_ += betas_[i];
    }

    //-
    setFeedbackCellField
    (
        fuelFeedbackCellField_,
        "fuelFeedbackZones"
    );
    setFeedbackCellField
    (
        coolFeedbackCellField_,
        "coolFeedbackZones"
    );
    setFeedbackCellField
    (
        structFeedbackCellField_,
        "structFeedbackZones"
    );
    setFeedbackCellField
    (
        drivelineFeedbackCellField_,
        "drivelineFeedbackZones"
    );

    //- Set flag in base class dict so it can be accessed by the GeN-Foam main
    bool GEMReactivity(GEMReactivityMap_.size() > 1);
    this->IOdictionary::set("GEM", GEMReactivity);

    if (GEMReactivity)
    {
        GEMSodiumLevelRef_ = nuclearData_.get<scalar>("GEMSodiumLevelRef");

        //- Check that GEMReactivityMap is indexed by descending sodium level
        //  values
        for (int i = 0; i < GEMReactivityMap_.size()-1; i++)
        {
            Pair<scalar> X0(GEMReactivityMap_[i]);
            Pair<scalar> X1(GEMReactivityMap_[i+1]);
            if (X0.first() <= X1.first())
            {
                FatalErrorInFunction
                << "GEMReactivityMap should be indexed by descending sodium "
                << "level values!" << exit(FatalError);
            }
        }
    }

    //- Some notes on modelling choices, for clarity
    Info<< "The pointKinetics neutronics model currently computes average "
        << "perturbed values for feedback fields (T fuel, cladding, etc.) "
        << "by weighting those over oneGroupFlux. This is not technically "
        << "correct as they should be weighted also via the adjoint flux."
        << "This will be addressed in future updates, but given that the "
        << "adjoint flux is equal to the flux if dealing with only one energy "
        << "group, this is deemed fine for now." << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::pointKineticNeutronics::setInitOneGroupFlux()
{
    initOneGroupFlux_ = oneGroupFlux_;
    initOneGroupFluxN_ =
        initOneGroupFlux_/fvc::domainIntegrate(initOneGroupFlux_);
    domainIntegratedInitOneGroupFluxN_ =
        fvc::domainIntegrate(sqr(initOneGroupFluxN_)).value();
}

void Foam::pointKineticNeutronics::setFeedbackCellField
(
    volScalarField& feedbackCellField,
    const word& keyword
)
{
    wordList feedbackZones
    (
        nuclearData_.lookupOrDefault(keyword, wordList())
    );
    if (feedbackZones.size() == 0)
    {
        forAll(feedbackCellField, i)
        {
            feedbackCellField[i] = 1.0;
        }
    }
    else
    {
        forAll(feedbackZones, i)
        {
            word zoneName(feedbackZones[i]);
            labelList zoneCells(mesh_.cellZones()[zoneName]);
            forAll(zoneCells, j)
            {
                label celli(zoneCells[j]);
                feedbackCellField[celli] = 1.0;
            }
        }
    }
}

Foam::scalar Foam::pointKineticNeutronics::calcDrivelineReactivity
(
    const scalar& drivelineExpansion
)
{
    scalar drivelineReactivity(0);

    label N(controlRodReactivityMap_.size());
    if (N > 1)
    {
        scalar xFirst(controlRodReactivityMap_[0].first());
        scalar xLast(controlRodReactivityMap_[N-1].first());

        //- If the map is indexed for ascending parameter values
        if (xFirst > xLast)
        {
            for (int i = 0; i < N-1; i++)
            {
                Pair<scalar> p0(controlRodReactivityMap_[i]);
                Pair<scalar> p1(controlRodReactivityMap_[i+1]);
                const scalar& x0(p0.first());
                const scalar& x1(p1.first());
                const scalar& y0(p0.second());
                const scalar& y1(p1.second());
                if
                (
                    drivelineExpansion <= x0
                and drivelineExpansion > x1
                )
                {
                    scalar m((y1-y0)/(x1-x0));
                    drivelineReactivity = m*(drivelineExpansion-x0) + y0;
                    break;
                }
            }
        }

        //- If the map is indexed for descending parameter values
        else
        {
            for (int i = 0; i < N-1; i++)
            {
                Pair<scalar> p0(controlRodReactivityMap_[i]);
                Pair<scalar> p1(controlRodReactivityMap_[i+1]);
                const scalar& x0(p0.first());
                const scalar& x1(p1.first());
                const scalar& y0(p0.second());
                const scalar& y1(p1.second());
                if
                (
                    drivelineExpansion > x0
                and drivelineExpansion <= x1
                )
                {
                    scalar m((y1-y0)/(x1-x0));
                    drivelineReactivity = m*(drivelineExpansion-x0) + y0;
                    break;
                }
            }
        }
    }

    return drivelineReactivity;
}

Foam::Pair<Foam::scalar>
Foam::pointKineticNeutronics::calcGEMLevelAndReactivity()
{
    scalar GEMSodiumLevel(0);
    scalar GEMReactivity(0);
    if (this->get<bool>("GEM"))
    {
        scalar intPhiFrac(0);
        forAll(phiFaces_, i)
        {
            const label& facei(phiFaces_[i]);
            intPhiFrac += (*phiOrig_)[facei]*phiMagSf_[i];
        }
        reduce(intPhiFrac, sumOp<scalar>());
        intPhiFrac /= intPhiRef_;

        //- Specific FFTF relationship between flow fraction and GEM sodium
        //  level
        GEMSodiumLevel =
            (265.0-539504/(2440.13+sqr(intPhiFrac*100)))/100
        -   GEMSodiumLevelRef_;

        label N(GEMReactivityMap_.size()-1);
        if (GEMReactivityMap_.size() > 1)
        {
            if (GEMSodiumLevel >= GEMReactivityMap_[0].first())
            {
                GEMReactivity = GEMReactivityMap_[0].second();
            }
            else if (GEMSodiumLevel <= GEMReactivityMap_[N].first())
            {
                GEMReactivity = GEMReactivityMap_[N].second();
            }
            else
            {
                for (int i = 0; i < GEMReactivityMap_.size()-1; i++)
                {
                    Pair<scalar> X0(GEMReactivityMap_[i]);
                    Pair<scalar> X1(GEMReactivityMap_[i+1]);
                    if
                    (
                        GEMSodiumLevel <= X0.first()
                    and GEMSodiumLevel > X1.first()
                    )
                    {
                        scalar m
                        (
                            (X1.second()-X0.second())/(X1.first()-X0.first())
                        );
                        GEMReactivity =
                            m*(GEMSodiumLevel-X0.first()) + X0.second();
                        break;
                    }
                }
            }
        }
    }

    return Pair<scalar>(GEMSodiumLevel, GEMReactivity);
}

void Foam::pointKineticNeutronics::getCouplingFieldRefs
(
    const objectRegistry& src,
    const meshToMesh& neutroToFluid
)
{
    //- Field names must reflect those defined in createCouplingFields.H
    TFuelOrig_ =
        src.findObject<volScalarField>("bafflelessTFuelAv");
    TCladOrig_ =
        src.findObject<volScalarField>("bafflelessTCladAv");
    TCoolOrig_ =
        src.findObject<volScalarField>("bafflelessTCool");
    rhoCoolOrig_ =
        src.findObject<volScalarField>("bafflelessRhoCool");
    TStructOrig_ =
        src.findObject<volScalarField>("bafflelessTStruct");

    //- Project thermalHydraulic powerDensities onto the neutronic ones to
    //  initialize them if the latters do not exist.
    //  Mirroring what happens in thermalHydraulics.interpolateCouplingFields:
    //  - if(liquidFuel), map powerDensityNeutronicsToLiquid_ to powerDensity_ 
    //    and powerDensityNeutronics_ to secondaryPowerDenisty_
    //  - else, map powerDensityNeutronics_ to powerDensity_ and 
    //    powerDensityNeutronicsToLiquid_ to secondaryPowerDenisty_
    if(liquidFuel_)
    {
        IOobject powerDensityHeader
        (
            "powerDensity",
            mesh_.time().timeName(),
            mesh_.time(),
            IOobject::NO_READ
        );

        if (!powerDensityHeader.typeHeaderOk<volScalarField>(true))
        {
            powerDensityToLiquidOrig_ =
                src.findObject<volScalarField>("bafflelessPowerDensityToLiquid");
            neutroToFluid.mapTgtToSrc
                (
                    *powerDensityToLiquidOrig_,
                    plusEqOp<scalar>(),
                    powerDensity_
                );
            powerDensity_.correctBoundaryConditions();
        }

        IOobject secondaryPowerDensityHeader
        (
            "secondaryPowerDensity",
            mesh_.time().timeName(),
            mesh_.time(),
            IOobject::NO_READ
        );

        if (!secondaryPowerDensityHeader.typeHeaderOk<volScalarField>(true))
        {
            powerDensityOrig_ =
                src.findObject<volScalarField>("bafflelessPowerDensity");
            neutroToFluid.mapTgtToSrc
                (
                    *powerDensityOrig_,
                    plusEqOp<scalar>(),
                    secondaryPowerDenisty_
                );
            secondaryPowerDenisty_.correctBoundaryConditions();
        }
    }
    else
    {
        IOobject powerDensityHeader
        (
            "powerDensity",
            mesh_.time().timeName(),
            mesh_.time(),
            IOobject::NO_READ
        );

        if (!powerDensityHeader.typeHeaderOk<volScalarField>(true))
        {
            powerDensityOrig_ =
                src.findObject<volScalarField>("bafflelessPowerDensity");
            neutroToFluid.mapTgtToSrc
                (
                    *powerDensityOrig_,
                    plusEqOp<scalar>(),
                    powerDensity_
                );
            powerDensity_.correctBoundaryConditions();
        }

        IOobject secondaryPowerDensityHeader
        (
            "secondaryPowerDensity",
            mesh_.time().timeName(),
            mesh_.time(),
            IOobject::NO_READ
        );

        if (!secondaryPowerDensityHeader.typeHeaderOk<volScalarField>(true))
        {
            powerDensityToLiquidOrig_ =
                src.findObject<volScalarField>("bafflelessPowerDensityToLiquid");
            neutroToFluid.mapTgtToSrc
                (
                    *powerDensityToLiquidOrig_,
                    plusEqOp<scalar>(),
                    secondaryPowerDenisty_
                );
            secondaryPowerDenisty_.correctBoundaryConditions();
        }
    }

    if (liquidFuel_)
    {
        UOrig_ =
            src.findObject<volVectorField>("bafflelessU");
        alphaOrig_ =
            src.findObject<volScalarField>("bafflelessAlpha");
        alphatOrig_ =
            src.findObject<volScalarField>("bafflelessAlphat");
        muOrig_ =
            src.findObject<volScalarField>("bafflelessMu");
    }
    else
    {
        UOrig_ = nullptr;
        alphaOrig_ = nullptr;
        alphatOrig_ = nullptr;
        muOrig_ = nullptr;
    }

    //- If doing GEM modelling
    if (this->get<bool>("GEM"))
    {
        phiOrig_ =
            src.findObject<surfaceScalarField>("phi");
        const labelList& faces
        (
            neutroToFluid.tgtRegion().faceZones()
            [
                nuclearData_.get<word>("GEMFlowFaceZone")
            ]
        );
        scalarField magSf(mag(neutroToFluid.tgtRegion().faceAreas()));
        forAll(faces, i)
        {
            const label& facei(faces[i]);
            phiFaces_.append(facei);
            phiMagSf_.append(magSf[facei]);
            intPhiRef_ += (*phiOrig_)[facei]*magSf[facei];
        }
        reduce(intPhiRef_, sumOp<scalar>());
    }
    else
    {
        phiOrig_ = nullptr;
    }

    //- The rest of this function is for initializing the reference values of
    //  the feedback parameters. If they are found in the dictionary, use
    //  those, otherwise compute them from the coupling fields (thus
    //  assuming that the simulation starts from a steady state)

    this->interpolateCouplingFields(neutroToFluid);

    #include "computeFeedbackFieldValues.H"

    //- TFuelRef limited as it appears in a fraction denominator if doing
    //  fastNeutrons (for the Doppler coeff)
    TFuelRef_ =
        max
        (
            reactorState_.lookupOrDefault<scalar>("TFuelRef", TFuelValue),
            SMALL
        );

    TCladRef_ =
        reactorState_.lookupOrDefault<scalar>("TCladRef", TCladValue);

    TCoolRef_ =
        reactorState_.lookupOrDefault<scalar>("TCoolRef", TCoolValue);

    rhoCoolRef_ =
        reactorState_.lookupOrDefault<scalar>("rhoCoolRef", rhoCoolValue);

    TStructRef_ =
        reactorState_.lookupOrDefault<scalar>("TStructRef", TStructValue);

    TDrivelineRef_ =
        reactorState_.lookupOrDefault<scalar>("TDrivelineRef", TDrivelineValue);

    //- Since a T*Ref_ are scalars, the value inside the dictionary is not
    //  updated at runTime automatically. The line below resets the SCALAR IN
    //  THE DICTIONARY at runTime so that the dictionary is written with the
    //  updated value.
    reactorState_.set("TFuelRef", TFuelRef_);
    reactorState_.set("TCladRef", TCladRef_);
    reactorState_.set("TCoolRef", TCoolRef_);
    reactorState_.set("rhoCoolRef", rhoCoolRef_);
    reactorState_.set("TStructRef", TStructRef_);
    reactorState_.set("TDrivelineRef", TDrivelineRef_);
    reactorState_.regIOobject::writeObject
    (
        IOstream::ASCII,
        // IOstream::currentVersion,
        // reactorState_.time().writeCompression(),
        true
    );

    #include "correctReactivity.H"

    //Info << endl << "pointKinetics (initial conditions): " << endl;
    //#include "pointKineticsInfo.H"
}

void Foam::pointKineticNeutronics::interpolateCouplingFields
(
    const meshToMesh& neutroToFluid
)
{
    neutroToFluid.mapTgtToSrc(*TFuelOrig_, plusEqOp<scalar>(), TFuel_);
    neutroToFluid.mapTgtToSrc(*TCladOrig_, plusEqOp<scalar>(), TClad_);
    neutroToFluid.mapTgtToSrc(*TCoolOrig_, plusEqOp<scalar>(), TCool_);
    neutroToFluid.mapTgtToSrc(*rhoCoolOrig_, plusEqOp<scalar>(), rhoCool_);
    neutroToFluid.mapTgtToSrc(*TStructOrig_, plusEqOp<scalar>(), TStruct_);

    if (liquidFuel_)
    {
        neutroToFluid.mapTgtToSrc(*UOrig_, plusEqOp<vector>(), UPtr_());
        neutroToFluid.mapTgtToSrc
        (
            *alphaOrig_,
            plusEqOp<scalar>(),
            alphaPtr_()
        );
        neutroToFluid.mapTgtToSrc
        (
            *alphatOrig_,
            plusEqOp<scalar>(),
            alphatPtr_()
        );
        neutroToFluid.mapTgtToSrc(*muOrig_, plusEqOp<scalar>(), muPtr_());
        phiPtr_() = fvc::flux(UPtr_());
        volScalarField diffCoeffOrig
        (
            (
                *alphatOrig_
            +   *muOrig_/ScNo_
            )/(*rhoCoolOrig_)
        );
        neutroToFluid.mapTgtToSrc
        (
            diffCoeffOrig,
            plusEqOp<scalar>(),
            diffCoeffPrecPtr_()
        );

        UPtr_().correctBoundaryConditions();
        alphaPtr_().correctBoundaryConditions();
        alphatPtr_().correctBoundaryConditions();
        diffCoeffPrecPtr_().correctBoundaryConditions();
    }


    TFuel_.correctBoundaryConditions();
    TClad_.correctBoundaryConditions();
    TCool_.correctBoundaryConditions();
    rhoCool_.correctBoundaryConditions();
    TStruct_.correctBoundaryConditions();
}

void Foam::pointKineticNeutronics::correct
(
    scalar& residual,
    label couplingIter
)
{
    if(!liquidFuel_)
    {
        #include "solvePointKinetics.H"
    }else
    {
        #include "solvePointKineticsLiquidFuel.H"
    }

}

// ************************************************************************* //
