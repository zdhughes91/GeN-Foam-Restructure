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
    Foam::structure

Description
    Class to handle the description of a mesh-unresolved structure. It is 
    derived from phaseBase, to represent the volume fraction of the structure.
    The key member for the understanding of this class is cellFields_. 
    The structure can be in fact though of as being composed by different
    regions, with properties defined only on certain cells of the mesh.
    However, as these region are meant NOT to spatially overlap, the 
    overall structure can be described by a single set of fields (Dh_, iApow_, 
    etc.) that are initialized/set/modified on a cell-by-cell basis.
    The cellFields_ member is thus a HashTable of scalarFields indexed by
    region name. These scalarFields are = 1.0 where the given
    region exists, and 0 elsewhere. What the region names are, that
    is set by the actual run-time-selected instance of the structure
    class.
    Currently, the only available structure is byZone, in which
    the regions are effectively cellZones, and the cellFields_ is
    thus a dicts of scalarFields representing where each cellZone exists.

    Author:  
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

SourceFiles
    structure.C

\*---------------------------------------------------------------------------*/

#ifndef structure_H
#define structure_H

#include "fvMatrix.H"
#include "phaseBase.H"
#include "HashPtrTable.H"
#include "runTimeSelectionTables.H"
#include "zeroGradientFvPatchFields.H"
#include "pump.H"
#include "heatExchanger.H"
#include "powerModel.H"
#include "powerOffCriterionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Forward declarations
class fluid;

/*---------------------------------------------------------------------------*\
                           Class structure Declaration
\*---------------------------------------------------------------------------*/

class structure
:
    public phaseBase
{
protected:

    //- Typedefs
    typedef HashTable
    <
        labelList,
        word,
        word::hash
    >   labelListTable;

    typedef HashTable
    <
        scalarField,
        word,
        word::hash
    >   scalarFieldTable;

    typedef HashTable
    <
        volScalarField,
        word,
        word::hash
    >   volScalarFieldTable;

    typedef HashTable
    <
        pump,
        word,
        word::hash
    >   pumpTable;

    typedef HashTable
    <
        heatExchanger,
        word,
        word::hash
    >   heatExchangerTable;

    typedef HashTable
    <
        autoPtr<powerModel>,
        word,
        word::hash
    >   powerModelTable;

    //- List of names of each group of cells into which this structure is split
    wordList regions_;

    //- Table of limited volume fraction fields indexed by region name
    //  Every field is equal to the overall alpha of structure only in those
    //  cells where the corresponding region is present, and is 0 
    //  elsewhere
    volScalarFieldTable alphaFields_;

    //- Table of cell scalarFields indexed by region name
    scalarFieldTable cellFields_;

    //- Table of cell labelLists indexed by region name
    labelListTable cellLists_;

    //- All cells in which the structure exists ()
    labelList cells_;

    //- Ref to powerDensity from neutronics (projected on fluid mesh)
    volScalarField& powerDensityNeutronics_;

    //- Hydraulic diameter
    volScalarField Dh_;

    //- Hydraulic diameter in the X, Y, Z directions in the local structure 
    //  frame of reference. Used exclusively for the computation of the local
    //  vector Reynolds number.
    volVectorField lDh_;

    //- Table of momentum sources (e.g. pumps)
    pumpTable pumps_;

    //- Table of heat exchangers
    heatExchangerTable heatExchangers_;

    //- Table of power models indexed by powerModel name
    powerModelTable powerModels_;

    //- Heat flux from active structures only. Positive if from structure to
    //  fluid
    volScalarField heatFlux_;

    //- Surface temperature of power regions
    volScalarField Tact_;

    //- Interfacial area of the powerModel region
    volScalarField iAact_;

    //- Heat exchanger surface temperature
    autoPtr<volScalarField> THXPtr_;

    //- Heat exchanger interfacial area (no difference between inner and outer,
    //  don't ask too much pls thnx. Maybe in the future)
    autoPtr<volScalarField> iAHXPtr_;

    //- Volume fraction of the passive structure (might be different than
    //  global structure volumefraction)
    volScalarField alphapas_;

    //- Surface temperature of the passive region
    volScalarField Tpas_;

    //- Temperature of the fuel (for coupling purposes)
    volScalarField TFuelAv_;

    //- Temperature of the cladding (for coupling purposes)
    volScalarField TCladAv_;

    //- Interfacial area of the passive region
    volScalarField iApas_;

    //- Volumetric heat capacity of the passive region times volume
    //  fraction of passive region
    volScalarField alphaRhoCppas_;

    //- Max between Tact_ and Tpas_;
    volScalarField Twall_;

    //- Tortuosity of the structure in the global reference frame
    autoPtr<volTensorField> tortuosityPtr_;

    //- Rotation matrix to move from global to local reference frame
    autoPtr<volTensorField> Rg2lPtr_;

    //- Rotation matrix to move from local to global reference frame
    autoPtr<volTensorField> Rl2gPtr_;

    //- Total momentum source
    autoPtr<volVectorField> momentumSourcePtr_;

    //- Power off criterion
    autoPtr<powerOffCriterionModel> powerOffCriterionModelPtr_;

public:

    TypeName("structure");

    //- Constructors

        structure
        (
            const dictionary& dict,
            const fvMesh& mesh,
            volScalarField& powerDensityNeutronics
        );

        //- Destructor
        virtual ~structure(){}

    //- Member Functions

        template<class Type>
        void setRegionField
        (
            const dictionary& dicts,
            GeometricField<Type, fvPatchField, volMesh>& field,
            word key
        ) const;

        //- Adjust pumps and return momentum source
        const volVectorField& momentumSource();

        //- By default heatExchangers operate on cellZones, so the way they are
        //  constructed is handled by a function here rather than in byZone
        void constructHeatExchangers();

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
        );

        tmp<volScalarField> explicitHeatSource
        (
            const fluid& fluid,
            const volScalarField& H
        ) const;

        tmp<fvScalarMatrix> linearizedSemiImplicitHeatSource
        (
            fluid& fluid,
            const volScalarField& H
        ) const;

        //- Function to access any member field of a powerModel, needed e.g.
        //  to get temperatures from the fuelModel

        //- Turn off all powerDensity fields of all the powerModels if the
        //  powerOffCriterionModel criterion is satisfied
        void checkPowerOff();

        //- Rotate a tensor field from the local reference frame to the global
        //  mesh one
        void localToGlobalRotateField(volTensorField& field) const;

        //- Rotate a vector field from the local reference frame to the global
        //  mesh one
        void localToGlobalRotateField(volVectorField& field) const;

        //- Rotate a tensor field from the global reference frame to the local 
        //  one
        void globalToLocalRotateField(volTensorField& field) const;

        //- Rotate a vector field from the global reference frame to the local 
        //  one
        void globalToLocalRotateField(volVectorField& field) const;

        //- Const access

            const wordList& regions() const
            {
                return regions_;
            }

            const volScalarFieldTable& alphaFields() const
            {
                return alphaFields_;
            }

            const scalarFieldTable& cellFields() const
            {
                return cellFields_;
            }

            const labelListTable& cellLists() const
            {
                return cellLists_;
            }

            const labelList& cells() const
            {
                return cells_;
            }

            const autoPtr<volTensorField>& tortuosityPtr() const
            {
                return tortuosityPtr_;
            }

            const volTensorField& tortuosity() const
            {
                return tortuosityPtr_();
            }

            volScalarField& powerDensityNeutronics()
            {
                return powerDensityNeutronics_;
            }

            const volScalarField& Dh() const
            {
                return Dh_;
            }

            const volVectorField& lDh() const
            {
                return lDh_;
            }

            const volScalarField& Tact() const
            {
                return Tact_;
            }

            const volScalarField& iAact() const
            {
                return iAact_;
            }

            const volScalarField& iApas() const
            {
                return iApas_;
            }

            const volScalarField& Tpas() const
            {
                return Tpas_;
            }

            volScalarField& TFuelAv()
            {
                return TFuelAv_;
            }

            volScalarField& TCladAv()
            {
                return TCladAv_;
            }

            const volScalarField& Twall() const
            {
                return Twall_;
            }

            const autoPtr<volTensorField>& Rg2lPtr() const
            {
                return Rg2lPtr_;
            }

            const autoPtr<volTensorField>& Rl2gPtr() const
            {
                return Rl2gPtr_;
            }

            bool hasLocalReferenceFrame() const
            {
            	return (Rg2lPtr_.valid() and Rl2gPtr_.valid());
            }

            const volTensorField& Rg2l() const
            {
                return Rg2lPtr_();
            }

            const volTensorField& Rl2g() const
            {
                return Rl2gPtr_();
            }

            const autoPtr<volVectorField>& momentumSourcePtr() const
            {
                return momentumSourcePtr_;
            }
};

//- Definition of templated static function in header, it's best practice
//  according to C++ 11 guidelines
template<class Type>
void structure::setRegionField 
(
    const dictionary& dicts,
    GeometricField<Type, fvPatchField, volMesh>& field,
    word key
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    IOobject fieldHeader
    (
        field.name(),
        mesh_.time().timeName(),
        mesh_,
        IOobject::NO_READ
    );

    //- field is a reference to a field that already exists. If it was
    //  created with a READ_IF_PRESENT or MUST_READ option, do nothing
    if 
    (
        ((
            field.readOpt() == IOobject::MUST_READ
        or  field.readOpt() == IOobject::READ_IF_PRESENT
        )
        and fieldHeader.typeHeaderOk<fieldType>(true))
    )
    {
        // DO NOTHING
        /*
        //  This part used to spatially limit the field only over the cells 
        //  overwhich it should be defined.
        //  Removed. This creates problems to setting powerDenistyNeutronics
        //  since each powerModel will set it to zero outside its cells 
        scalarField cellField(mesh_.cells().size(), 0);
        forAll(dicts.toc(), i)
        {
            word region(dicts.toc()[i]);
            const labelList& cells(cellLists_[region]);
            forAll(cells, j)
            {
                cellField[cells[j]] = 1.0;
            }
        }
        field.primitiveFieldRef() *= cellField;
        */
    }
    //- Otherwise, read the value the field should have on a 
    //  region-by-region basis
    else
    {
        forAll(dicts.toc(), i)
        {
            word region(dicts.toc()[i]);
            const labelList& cells(cellLists_[region]);
            Type value(dicts.subDict(region).get<Type>(key));
            forAll(cells, j)
            {
                field[cells[j]] = value;
            }
        }
    }
    field.correctBoundaryConditions();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
