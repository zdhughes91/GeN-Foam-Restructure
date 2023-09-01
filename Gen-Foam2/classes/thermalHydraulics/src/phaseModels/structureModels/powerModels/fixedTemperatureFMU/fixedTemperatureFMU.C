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

#if defined __has_include
#  if __has_include(<commDataLayer.H>) 
#    include <commDataLayer.H>
#    define isCommDataLayerIncluded
#  endif
#endif

#ifdef isCommDataLayerIncluded


#include "fixedTemperatureFMU.H"
#include "structure.H"
#include "addToRunTimeSelectionTable.H"
#include "commDataLayer.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace powerModels
{
    defineTypeNameAndDebug(fixedTemperatureFMU, 0);
    addToRunTimeSelectionTable
    (
        powerModel, 
        fixedTemperatureFMU, 
        powerModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::powerModels::fixedTemperatureFMU::fixedTemperatureFMU
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
    T_
    (
        IOobject
        (
            "T."+typeName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE //AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimTemperature, 0.0),
        zeroGradientFvPatchScalarField::typeName
    )
{
    this->setInterfacialArea();
    structure_.setRegionField(*this, T_, "T");


    forAll(this->toc(), regioni)
    {
        word region(this->toc()[regioni]);
        const dictionary& dict(this->subDict(region));

        // Preparing data to update temperature
        word temperatureKeyFromFMU("temperatureNameFromFMU");
        if (dict.found(temperatureKeyFromFMU))
        {
            const word temperatureNameFromFMU = dict.get<word>(temperatureKeyFromFMU);

            // Communicating with the FMU
            const Time& runTime = this->db().time();
            commDataLayer& data = commDataLayer::New(runTime); 
            // Store in data layer and set its initial value to the T 
            // in the dictionary      
            data.storeObj(
                dict.get<scalar>("T"),
                temperatureNameFromFMU,
                commDataLayer::causality::in
            );
            Info << "Using FMUs for the temperature in " << dict.dictName() << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::powerModels::fixedTemperatureFMU::~fixedTemperatureFMU()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::powerModels::fixedTemperatureFMU::temperatureUpdate() const
{
    forAll(this->toc(), regioni)
    {
        word region(this->toc()[regioni]);
        const dictionary& dict(this->subDict(region));
        
        word temperatureKeyFromFMU("temperatureNameFromFMU");
        const word temperatureNameFromFMU = dict.get<word>(temperatureKeyFromFMU);

        // Communicating with the FMU
        const Time& runTime = this->db().time();
        commDataLayer& data = commDataLayer::New(runTime);
        const scalar temperatureFromFMU =
            data.getObj<scalar>(temperatureNameFromFMU,commDataLayer::causality::in);
        
        //- Setup cellToRegion_ mapping
        const labelList& regionCells
        (
            structure_.cellLists()[region]
        );

        forAll(regionCells, i)
        {
            label celli(regionCells[i]);
            T_[celli] =  temperatureFromFMU;
        }
    }
}

void Foam::powerModels::fixedTemperatureFMU::correctT(volScalarField& T) const
{
    this->temperatureUpdate();
    forAll(cellList_, i)
    {
        label celli(cellList_[i]);
        T[celli] = T_[celli];
    }
}

void Foam::powerModels::fixedTemperatureFMU::powerOff()
{
    //- If you set iA to 0, the energy contribution from this powerModel to the
    //  fluid energy equation will be 0, equivalent to a "power" off scenario
    iA_ *= 0.0;
}

#endif
// ************************************************************************* //
