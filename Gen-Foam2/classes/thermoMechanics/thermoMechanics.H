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
    Foam::thermoMechanics

Description
    Base class for thermo-mechanical models.

    Reference publications:

    Carlo Fiorina, Ivor Clifford, Manuele Aufiero, Konstantin Mikityuk,
    "GeN-Foam: a novel OpenFOAM® based multi-physics solver for 2D/3D transient
    analysis of nuclear reactors", Nuclear Engineering and Design, submitted

    Carlo Fiorina, Konstantin Mikityuk, " Application of the new GeN-Foam multi-physics
    solver to the European Sodium Fast Reactor and verification against available codes",
    Proceedings of ICAPP 2015, May 03-06, 2015 - Nice (France), Paper 15226

    Author of this file (and associated .C or included .H files): 
    Carlo Fiorina <carlo.fiorina@outlook.com; carlo.fiorina@epfl.ch;>, EPFL (Switzerland)

SourceFiles
    thermoMechanics.C






\*---------------------------------------------------------------------------*/

#ifndef thermoMechanics_H
#define thermoMechanics_H

#include "volFields.H"
#include "rhoThermo.H"
#include "fvMatrices.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "coordinateSystem.H"
#include "fvCFD.H"
#include "turbulentFluidThermoModel.H"
#include "volPointInterpolation.H"
#include "meshToMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class thermoMechanics Declaration
\*---------------------------------------------------------------------------*/

class thermoMechanics
:
    public IOdictionary
{

protected:

    // Protected data

        //- Reference to mesh
        fvMesh& mesh_;

        //- Mesh displacement
        volVectorField meshDisp_;

        //- Initialize residual
        double initialResidual_;        

private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        thermoMechanics(const thermoMechanics&);

        //- Disallow default bitwise assignment
        void operator=(const thermoMechanics&);


public:

    //- Runtime type information
    TypeName("thermoMechanics")

#ifndef SWIG
        declareRunTimeSelectionTable
        (
            autoPtr,
            thermoMechanics,
            dictionary,
            (
                fvMesh& mesh
            ),
            (mesh)
        );
#endif

    // Constructors

        //- Construct from mesh
        thermoMechanics(fvMesh& mesh);


    // Selectors

        //- Return a reference to the selected model
        static autoPtr<thermoMechanics> New
        (
            fvMesh& mesh
        );


    //- Destructor
    	virtual ~thermoMechanics(){}


    // Member Functions   	

        //- Substitutes the old getFields. Set local class references to global
        //  coupling fields, to be called only once after construction. The 
        //  fields are set via lookup in the src objectRegistry (i.e. 
        //  fluidMeshNB for the current state of GeN-Foam)
        virtual void getCouplingFieldRefs
        (
            const objectRegistry& src,
            const meshToMesh& mechToFluid
        ) = 0;

        //- The old getFields, yet operates on the local refs without having
        //  the need for all the fields to be passed. This effectively bypasses
        //  the old limitation of getFields, i.e. if someone wanted to
        //  implement a new thermoMechanical model that requires more fields 
        //  for the coupling, then they'd need to modify the getFields function
        //  in ALL of the other thermoMechanical models too
        virtual void interpolateCouplingFields
        (
            const meshToMesh& mechToFluid
        ) = 0;

        //- Correct/update the properties
        virtual void correct(scalar& residual) = 0;    

        //- Needed to remove abstractedness from IOdictionary
        virtual bool writeData(Ostream& os) const
        {
            return os.good();
        }

        //- Access

            //- Return mesh displacement
            const volVectorField& meshDisp() const
            {
                return meshDisp_;
            } 

            //- Return initialResidual_
            double initialResidual() const
            {
                return initialResidual_;
            };   
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
