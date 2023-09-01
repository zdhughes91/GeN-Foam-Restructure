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

#include "fixedTemperature.H"
#include "structure.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace powerModels
{
    defineTypeNameAndDebug(fixedTemperature, 0);
    addToRunTimeSelectionTable
    (
        powerModel, 
        fixedTemperature, 
        powerModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::powerModels::fixedTemperature::fixedTemperature
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
    ),
    timeProfile_(this->toc().size()),
    t0_(this->toc().size()),
    timeDependent_(this->toc().size())
{
    this->setInterfacialArea();
    structure_.setRegionField(*this, T_, "T");

    forAll(this->toc(), regioni)
    {
        word region(this->toc()[regioni]);
        const dictionary& dict(this->subDict(region));

        // Preparing data to update temperature
        timeDependent_[regioni] = false;
        word timeProfileDictName("temperatureTimeProfile");
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

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::powerModels::fixedTemperature::~fixedTemperature()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::powerModels::fixedTemperature::temperatureUpdate() const
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
            scalar timeDependentTemperature(timeProfile_[regioni].value(t));
            forAll(regionCells, i)
            {
                label celli(regionCells[i]);
                T_[celli] =  timeDependentTemperature;
            }          
        }
    }
}

void Foam::powerModels::fixedTemperature::correctT(volScalarField& T) const
{
    this->temperatureUpdate();
    forAll(cellList_, i)
    {
        label celli(cellList_[i]);
        T[celli] = T_[celli];
    }
}

void Foam::powerModels::fixedTemperature::powerOff()
{
    //- If you set iA to 0, the energy contribution from this powerModel to the
    //  fluid energy equation will be 0, equivalent to a "power" off scenario
    iA_ *= 0.0;
}

// ************************************************************************* //
