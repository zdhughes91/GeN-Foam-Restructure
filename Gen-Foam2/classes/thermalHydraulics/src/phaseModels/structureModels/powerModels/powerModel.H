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
    Foam::powerModel

Description
    Run-time selectable class to handle the thermal description (i.e. energy 
    source/sinks) of a structure.
    It returns a heat source term to the main solver be used in the energy 
    equations.

    Author: 
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

SourceFiles
    powerModel.C

\*---------------------------------------------------------------------------*/

#ifndef powerModel_H
#define powerModel_H

#include "dictionary.H"
#include "volFields.H"
#include "runTimeSelectionTables.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Forward declarations
class structure;

/*---------------------------------------------------------------------------*\
                           Class powerModel Declaration
\*---------------------------------------------------------------------------*/

class powerModel
:
    public IOdictionary
{
protected:

    //- Ref to parent structure
    structure& structure_;

    //- Reference to mesh
    const fvMesh& mesh_;

    //- List of cells in which this powerModel exists
    labelList cellList_;

    //- Marker field for this power model, = 1 in all the cells that are
    //  listed in cellsList, and 0 elsewhere
    scalarField cellField_;

    //- Interfacial area of the powerModel
    volScalarField iA_;

    //- Volume fraction of the powerModel (can in principle be different from
    //  the structure alpha as not the entirety of a structure region might
    //  be producing power)
    volScalarField alpha_;

public:

    //- Runtime type information
    TypeName("powerModel");

    //- Runtime selection table declaration
    declareRunTimeSelectionTable
    (
        autoPtr,
        powerModel,
        powerModels,
        (
            structure& structure,
            const dictionary& dicts
        ),
        (structure, dicts)
    );

    // Constructors

        //- Construct over specific cells
        powerModel
        (
            structure& structure,
            const dictionary& dicts
        );

        //- Destructor
        virtual ~powerModel(){}

    // Selectors

        static autoPtr<powerModel> New
        (
            structure& structure,
            const dictionary& dicts
        );

    // Members

        //- For some models it is mandatory, others (e.g. the pin ones) should
        //  allow for the option for the interfacial area to be computed
        //  starting from other powerModel geometric properties (e.g. again,
        //  for a pin geometry, the interfacial area is directly linked to the
        //  volume-averaged pin volume fraction and the pin outermost radius).
        //  This, this functions inits interfacialArea as a must-read when
        //  called from the derived calsses constructors, while derived classes
        //  that do not use this will construct the interfacial area their own
        //  way
        virtual void setInterfacialArea();

        /*
        volScalarField& initFieldInTable
        (
            word name, 
            dimensionedScalar value,
            bool readIfPresent,
            bool autoWrite
        );
        */

        //- In two-phase solvers :
        //  HT = H1*frac1*T1 + H2*frac2*T2
        //  H =  H1*frac1 + H2*frac2
        //- In a mono-phase solver :
        //  HT = H*T
        //  H = H
        virtual void correct
        ( 
            const volScalarField& HT, 
            const volScalarField& H
        ) = 0;

        //- Set provided field with structure surface temperature
        virtual void correctT(volScalarField& T) const = 0;

        //- Turn off power
        virtual void powerOff(){}

        //- Remove abstractness inherited from regIOobject
        virtual bool writeData(Ostream& os) const 
        {
            return os.good();
        }

        // Access

            const volScalarField& iA() const
            {
                return iA_;
            }

            const scalarField& cellField() const
            {
                return cellField_;
            }

            structure& structureRef() 
            {
                return structure_;
            }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
