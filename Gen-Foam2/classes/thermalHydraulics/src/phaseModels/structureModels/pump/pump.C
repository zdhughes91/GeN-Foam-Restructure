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

#include "pump.H"

//- From forward declarations
#include "structure.H"
//#include "commDataLayer.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pump, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pump::pump
(
    const fvMesh& mesh,
    const dictionary& dict,
    const labelList& cellList
)
:
    IOdictionary
    (
        IOobject
        (
            typeName,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dict
    ),
    mesh_(mesh),
    cellList_(cellList),
    pumpValue_(this->get<vector>("momentumSource")),
    timeProfilePtr_(nullptr),
    t0_(0.0),
    timeProfile_(false),
    pumpMultiplierNameFromFMU_("momentumSourceCoupled"),
    fromFMU_(false)
{
    Info << "Creating pump in " << dict.dictName() << endl;

    word timeProfileDictName("momentumSourceTimeProfile");
    word pumpMultiplierKeyFromFMU("pumpMultiplierNameFromFMU");

    if (this->found(timeProfileDictName))
    {
        const dictionary& timeProfileDict(dict.subDict(timeProfileDictName));
        word type
        (
            timeProfileDict.get<word>("type")
        );
        timeProfilePtr_.reset        
        (
            Function1<scalar>::New
            (
                type,
                timeProfileDict,
                type
            )
        );
        timeProfile_ = true;
        t0_ = timeProfileDict.lookupOrDefault("startTime", 0.0);

        Info << "Using a time profile for the pump in " << dict.dictName() << endl;
    }
    #ifdef isCommDataLayerIncluded
    if (this->found(pumpMultiplierKeyFromFMU))
    {
        pumpMultiplierNameFromFMU_ = this->get<word>(pumpMultiplierKeyFromFMU);
        // Communicating with the FMU
        const Time& runTime = this->db().time();
        commDataLayer& data = commDataLayer::New(runTime); 
        // Store in data layer and set its initial value to 1       
        data.storeObj(
            scalar(1.0),
            pumpMultiplierNameFromFMU_,
            commDataLayer::causality::in
            );
        fromFMU_ = true;

        Info << "Using FMUs for the pump in " << dict.dictName() << endl;
    }

    if(timeProfile_ && fromFMU_)
    {
        Info << "WARNING: Both time profile and FMU coupling provided for the pump in " 
        << dict.dictName() << endl
        << "GeN-Foam will use the FMU" << endl;
    }
    #endif

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::pump::correct
( 
    volVectorField& momentumSource
)
{
    vector pumpValue(pumpValue_);

    if (timeProfile_)
    {
        scalar t(mesh_.time().timeOutputValue()-t0_);
        pumpValue = pumpValue_ * timeProfilePtr_->value(t);
    }
    #ifdef isCommDataLayerIncluded
    if (fromFMU_)
    {
        const Time& runTime = this->db().time();
        commDataLayer& data = commDataLayer::New(runTime);
        const scalar pumpMultiplierFromFMU =
            data.getObj<scalar>(pumpMultiplierNameFromFMU_,commDataLayer::causality::in);
        //update the vector field by adjusting the magnitude
        pumpValue = pumpValue_ * pumpMultiplierFromFMU;
      
    }
    #endif
    forAll(cellList_, i)
    {
        const label& celli(cellList_[i]);
        momentumSource[celli] = pumpValue;
    }

    momentumSource.correctBoundaryConditions();
}

// ************************************************************************* //
