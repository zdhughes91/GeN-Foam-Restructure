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
    Foam::phaseChangeModels::heatDrivenPhaseChange

Description
    Phase change models that compute mass transfer based on interfacial
    heat fluxes. There are three possible modes that can be selected:

    -   conductionLimited: dmdt_ is computed so to conserve total energy
        transfer across interface if the interfacial heat fluxes do not
        balance out. It q1i and q2i are the respective heat fluxes form the
        bulk of phase 1/2 to the interface, then dmdt_ = (q1i+q2i)/L with L
        being the latent heat. This is the approach found in TRACE 
        (check theory manual, https://www.nrc.gov/docs/ML0710/ML071000097.pdf)
        but it is not suited if there is a massive difference in the volumetric
        heat capacity (J/K/m3) of the phases (such as in the case of sodium, due
        to the density difference of a factor ~ 2000). It can be very unstable 
        and not converge at all in such scenarios;

    -   twoPhaseDriven: evaporation is driven uniquely by liquid superheat and
        condensation is driven uniquely by vapor undercooling. There is no need
        to specify which phase is liquid and which is vapour in the
        phaseProperties dict, as this information is contained in the latent 
        heat L (positive if phase1 is liquid and 2 is vapour, negative 
        otherwise);

    -   onePhaseDriven: both evaporation and condensation are driven uniquely
        by one phase's superheat or cooling. Implemented as a counterpart to
        twoPhaseDriven.

    The latent heat is computed as a difference of the enthalpies of formation
    (specified in the thermoPhysicalProperties of each phase under the Hf
    keyword). Thus, the vapour phase is the phase with the highest enthalpy of
    formation. The latent heat can be adjusted as described in the TRACE theory
    manual to avoid thermal-run-aways in particular scenarios. This feature can
    be enabled via the correctLatentHeat flag in the phaseChangeModel subDict,
    but it is not recommended if dealing with phases with a low volumetric
    heat capacity, as described before (e.g. Sodium).

    Dict example
    \verbatim
    phaseChangeModel
    {
        type            heatDriven;
        mode            onePhaseDriven; //- conductionLimited;//twoPhaseDriven;
        drivingPhase    "phaseName";    //- (required only for onePhaseDriven)
        correctLatentHeat   false;      //- defaults to false

        saturationModel                 //- Unrelated to the heatDriven model
        {
            type       "BrowningPotter"; 
        }
    }
    \endverbatim

SourceFiles
    heatDrivenPhaseChange.C

\*---------------------------------------------------------------------------*/

#ifndef heatDrivenPhaseChange_H
#define heatDrivenPhaseChange_H

#include "phaseChangeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseChangeModels
{

/*---------------------------------------------------------------------------*\
                         Class heatDrivenPhaseChange Declaration
\*---------------------------------------------------------------------------*/

class heatDrivenPhaseChange
:
    public phaseChangeModel
{
public:

    enum mode
    {
        conductionLimited, 
        onePhaseDriven,
        twoPhaseDriven,
        mixedDriven
    };

    static const Enum<mode> modeNames_;

protected:

    //-
    const mode mode_;

    //-
    const word drivingPhaseName_;

public:

    TypeName("heatDriven");

    //- Constructors

        //- Construct from a dictionary
        heatDrivenPhaseChange
        (
            FFPair& pair,
            const dictionary& dict
        );

        //- Destructor
        virtual ~heatDrivenPhaseChange(){}

    //- Member Functions

        //-
        virtual void correctInterfacialDmdt();
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace phaseChangeModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
