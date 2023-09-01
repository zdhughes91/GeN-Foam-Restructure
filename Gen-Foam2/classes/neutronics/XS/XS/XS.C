/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "XS.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(XS, 0);
    defineRunTimeSelectionTable(XS, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::XS::XS
(
    const fvMesh& mesh
)
:
    mesh_(mesh),
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
    nuclearDataRadialExp_
    (
        IOobject
        (
            "nuclearDataRadialExp",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    nuclearDataAxialExp_
    (
        IOobject
        (
            "nuclearDataAxialExp",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    nuclearDataFuelTemp_
    (
        IOobject
        (
            "nuclearDataFuelTemp",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    nuclearDataRhoCool_
    (
        IOobject
        (
            "nuclearDataRhoCool",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    nuclearDataTCool_
    (
        IOobject
        (
            "nuclearDataTCool",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    nuclearDataCladExp_
    (
        IOobject
        (
            "nuclearDataCladExp",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    energyGroups_(nuclearData_.lookupOrDefault("energyGroups",1)),
    precGroups_(nuclearData_.lookupOrDefault("precGroups",1)),
    legendreMoments_(1+nuclearData_.lookupOrDefault("legendreMoments",0)),
    IV_(energyGroups_),
    D_(energyGroups_),
    nuSigmaEff_(energyGroups_),
    sigmaPow_(energyGroups_),
    sigmaDisapp_(energyGroups_),
    sigmaFromTo_(legendreMoments_),
    chiPrompt_(energyGroups_),
    chiDelayed_(energyGroups_),
    Beta_(precGroups_),
    BetaTot_
    (
        IOobject
        (
            "BetaTot",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimensionSet(0,0,0,0,0,0,0), 1.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    lambda_(precGroups_),
    fuelFraction_
    (
        IOobject
        (
            "fuelFraction",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimensionSet(0,0,0,0,0,0,0), 1.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    secondaryPowerVolumeFraction_
    (
        IOobject
        (
            "secondaryPowerVolumeFraction",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimensionSet(0,0,0,0,0,0,0), 1.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    fractionToSecondaryPower_
    (
        IOobject
        (
            "fractionToSecondaryPower",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimensionSet(0,0,0,0,0,0,0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    discFactor_(energyGroups_),
    ScNo_(nuclearData_.lookupOrDefault("ScNo",1.0)),
    entries_(nuclearData_.lookup("zones")),
    zoneNumber_(entries_.size()),    
    fuelFractionList_(zoneNumber_),
    secondaryPowerVolumeFractionList_(zoneNumber_),
    fractionToSecondaryPowerList_(zoneNumber_),
    dfAdjustList_(zoneNumber_),
    discFactorList_(zoneNumber_),
    integralFluxList_(zoneNumber_),
    fastNeutrons_(nuclearData_.lookupOrDefault("fastNeutrons",true)),
    adjustDiscFactors_(nuclearData_.lookupOrDefault("adjustDiscFactors", false)),
    useGivenDiscFactors_(nuclearData_.lookupOrDefault("useGivenDiscFactors", false)),
    groupsWoDF_(nuclearData_.lookupOrDefault<List<int> >("groupsWoDF", List<int>())),
    doNotParametrize_(nuclearData_.lookupOrDefault<List<int> >("doNotParametrize", List<int>())),
    IVList_(zoneNumber_),
    chiPromptList_(zoneNumber_),
    chiDelayedList_(zoneNumber_),
    BetaList_(zoneNumber_),
    BetaTotList_(zoneNumber_),
    lambdaList_(zoneNumber_),
    DList_(zoneNumber_),
    nuSigmaEffList_(zoneNumber_),
    sigmaPowList_(zoneNumber_),
    sigmaDisappList_(zoneNumber_),
    sigmaFromToList_(zoneNumber_), 
    TfuelRef_(nuclearDataFuelTemp_.lookupOrDefault("TfuelRef",900.0)),
    TfuelPerturbed_(nuclearDataFuelTemp_.lookupOrDefault("TfuelPerturbed",1200.0)),
    fuelTempDList_(zoneNumber_),
    fuelTempNuSigmaEffList_(zoneNumber_),
    fuelTempSigmaPowList_(zoneNumber_),
    fuelTempSigmaDisappList_(zoneNumber_),
    fuelTempSigmaFromToList_(zoneNumber_),
    AxExp_(nuclearDataAxialExp_.lookupOrDefault("expansionFromNominal",1.0)),
    axialExpDList_(zoneNumber_),
    axialExpNuSigmaEffList_(zoneNumber_),
    axialExpSigmaPowList_(zoneNumber_),
    axialExpSigmaDisappList_(zoneNumber_),
    axialExpSigmaFromToList_(zoneNumber_),
    RadExp_(nuclearDataRadialExp_.lookupOrDefault("expansionFromNominal",1.0)),
    axialOrientation_(nuclearDataRadialExp_.lookupOrDefault("axialOrientation",vector(0.0, 0.0, 1.0))),
    radialExpDList_(zoneNumber_),
    radialExpNuSigmaEffList_(zoneNumber_),
    radialExpSigmaPowList_(zoneNumber_),
    radialExpSigmaDisappList_(zoneNumber_),
    radialExpSigmaFromToList_(zoneNumber_),
    rhoCoolRef_(nuclearDataRhoCool_.lookupOrDefault("rhoCoolRef",860.0)),
    rhoCoolPerturbed_(nuclearDataRhoCool_.lookupOrDefault("rhoCoolPerturbed",1.0)),
    rhoCoolDList_(zoneNumber_),
    rhoCoolNuSigmaEffList_(zoneNumber_),
    rhoCoolSigmaPowList_(zoneNumber_),
    rhoCoolSigmaDisappList_(zoneNumber_),
    rhoCoolSigmaFromToList_(zoneNumber_),
    TCoolRef_(nuclearDataTCool_.lookupOrDefault("TCoolRef",900.0)),
    TCoolPerturbed_(nuclearDataTCool_.lookupOrDefault("TCoolPerturbed",1200.0)),
    TCoolDList_(zoneNumber_),
    TCoolNuSigmaEffList_(zoneNumber_),
    TCoolSigmaPowList_(zoneNumber_),
    TCoolSigmaDisappList_(zoneNumber_),
    TCoolSigmaFromToList_(zoneNumber_),
    TcladRef_(nuclearDataCladExp_.lookupOrDefault("TcladRef",900.0)),
    TcladPerturbed_(nuclearDataCladExp_.lookupOrDefault("TcladPerturbed",1200.0)),
    cladExpDList_(zoneNumber_),
    cladExpNuSigmaEffList_(zoneNumber_),
    cladExpSigmaPowList_(zoneNumber_),
    cladExpSigmaDisappList_(zoneNumber_),
    cladExpSigmaFromToList_(zoneNumber_),
    CRmove_
    (
        IOobject
        (
            "CRmove",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    CRentries_(CRmove_.lookup("zones")),
    CRNumber_(CRmove_.lookup("zones").size()),
    CRstart_(CRNumber_),
    CRfinish_(CRNumber_),
    CRspeed_(CRNumber_),
    CRFollowerName_(CRNumber_),
    CRinitialPosition_(CRNumber_),
    CRposition_(CRNumber_),
    initialDistanceFromMeshCR_(CRNumber_)        
{
    #include "readNuclearData.H"
    #include "createXSfields.H"    
    init();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::XS::~XS()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::XS::correct
(
    const volScalarField& Tfuel, 
    const volScalarField& Tclad, 
    const volScalarField& rhoCool, 
    const volScalarField& TCool,
    const volVectorField& Disp
)
{
    #include "setNeutronicsVariables.H"
}

void Foam::XS::init()
{
    #include "setNeutronicsConstants.H"
    #include "setPrecConst.H"
}

void Foam::XS::adjustDiscFactors(const PtrList<volScalarField>& fluxStar)
{
    if (adjustDiscFactors_)
    {
        #include "adjustDiscFactors.H"
    }
}




// ************************************************************************* //
