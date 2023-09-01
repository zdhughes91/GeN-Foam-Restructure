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

Class
    Foam::latentHeatModel

Description
    Class used to describe the latent heat of fluid that is undergoing phase
    change. It can depend on some form of system temperature and system
    pressure.

SourceFiles
    latentHeatModel.C
    newSaturationModel.C

\*---------------------------------------------------------------------------*/

#ifndef latentHeatModel_H
#define latentHeatModel_H

#include "fluid.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class FFPair;
class phaseChangeModel;

/*---------------------------------------------------------------------------*\
                         Class latentHeatModel Declaration
\*---------------------------------------------------------------------------*/

class latentHeatModel
:
    public IOdictionary
{
protected:

    //- Ref to mesh
    const fvMesh& mesh_;

    //- Ref to liquid
    const fluid& liquid_;

    //- Ref to vapour
    const fluid& vapour_;

    //- Ref to system pressure
    const volScalarField& p_;

    //- Ref to mass transfer term, only used by the adjust_ feature
    const volScalarField& dmdt_;

    //- Ref to interfacial temperature, only used by the adjust_ feature
    const volScalarField& iT_;

    //  Adjust the latent heat of vaporization in a manner consistent
    //  with how massTransfer-related enthalpy is added to the energy
    //  equations. In particular:
    //  * when evaporating, mass is removed from the liquid at its current
    //    enthalpy and is added to the vapour at vapour saturation
    //    entalphy;
    //  * when condensing, mass is removed from the vapour at its current
    //    enthalpy and is added to the liquid at liquid saturation
    //    enthalpy.
    //    For example, in the case of condensing super-heated steam, its
    //    enthalpy will be larger than saturation enthalpy. However, if
    //    you remove mass from the steam at an enthalpy lower than its
    //    current enthalpy, the remaining steam will inevitably be
    //    hotter.
    //    Thus:
    //  * when evaporating, the latent heat will be decreased by the difference
    //    between the current liquid enthalpy and the saturation liquid
    //    enthalpy;
    //  * when condensing, the latent heat will be increased by the difference
    //    between the current vapour enthalpy and the saturation vapour
    //    enthalpy.
    const bool adjust_;

    //- L > 0 if fluid1_ is liquid and fluid2_ is vapour, < 0 if vice-versa.
    //  This comes from the fact that the rest of the code is hard-coded so
    //  that dmdt represents a mass transfer from fluid1 to fluid2 if positive,
    //  the opposite when negative
    const int LSign_;

public:

    TypeName("latentHeatModel");

    //- Declare runtime construction
    declareRunTimeSelectionTable
    (
        autoPtr,
        latentHeatModel,
        latentHeatModels,
        (
            const phaseChangeModel& pcm,
            const dictionary& dict, 
            const objectRegistry& objReg
        ),
        (pcm, dict, objReg)
    );

    //- Constructors

        latentHeatModel
        (
            const phaseChangeModel& pcm,
            const dictionary& dict, 
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~latentHeatModel(){}

    //- Selectors

        static autoPtr<latentHeatModel> New
        (
            const phaseChangeModel& pcm,
            const dictionary& dict, 
            const objectRegistry& objReg
        );

    //- Member Functions

        //- Update latent heat
        virtual scalar value(const label& celli) const;

        //-
        virtual void correctField(volScalarField& L) const;

        //-
        void adjust(volScalarField& L) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
