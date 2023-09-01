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

#include "heatExchanger.H"

//- From forward declarations
#include "fluid.H"
#include "structure.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(heatExchanger, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatExchanger::heatExchanger
(
    const fvMesh& mesh,
    const dictionary& dict
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
    iA_(this->get<scalar>("volumetricArea")),
    Hw_(max(this->get<scalar>("wallConductance"), 1e-69)),
    primaryCells_(0),
    secondaryCells_(0)
{
    Info << "Creating heatExchanger: " << dict.dictName() << endl;

    //- Read primary and secondary cellZones and their cells
    primaryCells_ = mesh.cellZones()[this->get<word>("primary")];
    secondaryCells_ = mesh.cellZones()[this->get<word>("secondary")];

    //- Calculate displacement vector that, if applied to the primary
    //  cellZone, would translate it to the secondary cellZone. This assumes
    //  that the cellZones have the same shape and volume. If that is not
    //  the case, a check is done afterwards (currently only on the volume, not
    //  on the shape, which I will implement via comparing bounding boxes at
    //  some point)
    vector pCOV(vector::zero);
    vector sCOV(vector::zero);
    const vectorField& C(mesh_.C());
    const scalarField& V(mesh_.V());
    scalar pV(0);
    scalar sV(0);
    forAll(primaryCells_, i)
    {
        const label& celli(primaryCells_[i]);
        pCOV += V[celli]*C[celli];
        pV += V[celli];
    }
    reduce(pCOV, sumOp<vector>());
    reduce(pV, sumOp<scalar>());
    pCOV /= pV;
    forAll(secondaryCells_, i)
    {
        const label& celli(secondaryCells_[i]);
        sCOV += V[celli]*C[celli];
        sV += V[celli];
    }
    reduce(sCOV, sumOp<vector>());
    reduce(sV, sumOp<scalar>());
    sCOV /= sV;
    vector delta(sCOV-pCOV);
    scalar errV((pV-sV)/pV);
    if ((pV-sV)/pV >= 1e-2)
    {
        FatalErrorInFunction
            << "Primary and secondary zone volumes differ by " << errV
            << exit(FatalError);
    }
    Info<< "Heat exchanger " << dict.dictName() << " mapping: applied "
        << "translation of " << delta << " m" << endl;

    //- Read mesh geometric data
    const pointField& pointsRef(mesh_.points());
    const faceList& facesRef(mesh_.faces());
    const cellList& cellsRef(mesh_.cells());

    //- Init new mesh geometric data
    pointField points(pointsRef.size());
    faceList faces(facesRef.size());
    cellList cells(cellsRef.size());

    //- Create translated mesh
    //  First, copy all points and apply translation
    forAll(points, i)
    {
        points[i] = pointsRef[i] + delta;
    }

    //- Set faces, cells as copies of the original indexing (yet the indexing
    //  applies to the translated points, so all the resulting faces and cells
    //  will result translated as well)
    forAll(faces, i)
    {
        faces[i] = facesRef[i];
    }
    forAll(cells, i)
    {
        cells[i] = cellsRef[i];
    }

    //- Assemble new translated mesh
    autoPtr<fvMesh> translatedMeshPtr;
    translatedMeshPtr.reset
    (
        new Foam::fvMesh
        (
            IOobject
            (
                "translatedMesh",
                mesh_.time().constant(),
                mesh_.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            std::move(points),
            std::move(faces),
            std::move(cells)
        )
    );
    fvMesh& translatedMesh = translatedMeshPtr();

    //- Assemble mapping
    mappingPtr_.reset
    (
        new meshToMesh
        (
            translatedMesh,
            mesh,
            Foam::meshToMesh::interpolationMethod
            (
                this->lookupOrDefault
                (
                    "interpolationMethod", 2
                )
            ),
            Foam::meshToMesh::procMapMethod::pmAABB,
            false
        )
    );

    //  I do not care about the translatedMesh anymore after this point, I only
    //  needed the mapping!
    translatedMeshPtr.clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::heatExchanger::iA() const
{
    tmp<volScalarField> tiA
    (
        new volScalarField
        (
            IOobject
            (
                "",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("", dimArea/dimVol, 0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    volScalarField& iA(tiA.ref());
    
    forAll(primaryCells_, i)
    {
        const label& celli(primaryCells_[i]);
        iA[celli] = iA_;
    }
    forAll(secondaryCells_, i)
    {
        const label& celli(secondaryCells_[i]);
        iA[celli] = iA_;
    }
    iA.correctBoundaryConditions();
    
    return tiA;
}

/*

Theoretical explanation & implementation

I assume that the wall separating the two fluids (which might as well be 
two-phase mixtures) has a certain thermal conductance Hw (provided by the
user, which for e.g. a planar sheet of metal, or a thin-walled tube,
reduced to Hw = thermal conductivity/wall thickness) and that can be
characterised by two surface temperatures, one on the primary side (I will
call it Tp) and one on the secondary side (Ts). The THX field that is being
passed is exactly the HX surface temperature, so that, it those mesh cells
that belong to the primary cellZone, it will be Tp, and Ts in the mesh
cells that belong to the secondary. The key idea is to set Tp and Ts
so to ensure heat flow (flux to, as I assume the same HX volumetric area
on both the primary and secondary sides) conservation across the HX wall.
Furthermore, I assume that the HX wall has no thermal inertia per-se 
(even though I suppose it would not be overly complicated to model that as
well, given the current framework). For exaple, let us focus on the primary
side. The heat flux through the wall must be equal to the heat flux
through the fluid. For a two-phase mixture, I assume that the heat is
distributed between the two phases (labelled 1 and 2) according to the
fluidPartitionModel (the frac variables seen elsewhere in the code).
In particular the heat flux from the secondary side of the wall to the
primary will be:

Q_(s->p) = Hw*(Ts-Tp)                                                       (1)

while the heat flow from the primary wall surface to the primary two-phase
mixture will be:

Q_(p->mix) = f1p*H1p*(Tp-T1p)+f2p*H2p*(Tp-T2p)                              (2)

in which the p subscript of all the quantities on the right refers to them
being the values of the fluid mixture on the primary side. This heat flux
repartition is the same one used for powerModels and passiveStructure source
terms elsewhere in the code. By equating (1) and (2) to ensure flux 
conservation, we can obtain an expression that relates the secondary side wall
temperature Ts to the primary wall temperature Ts as well as two-phase mixture
properties on the primary side. In particular, if we consider the expression
for both wall sides, we obtain the following system of two equations for two 
unknowns (Tp, Ts).

Ts = ((Hw+f1p*H1p+f2p*H2p)/Hw)*(Tp) - (f1p*H1p*T1p+f2p*H2p*T2p)/Hw          (3)
Tp = ((Hw+f1s*H1s+f2s*H2s)/Hw)*(Ts) - (f1p*H1s*T1s+f2s*H2s*T2s)/Hw          

The solutions to this system of equations consist of:

Tp = (Bp*As+Bs)/(Ap*As-1)                                                   (4)
Ts = (Bs*Ap+Bp)/(As*Ap-1) 

with (i=p or i=s):

Ai = (Hw+(f1i*H1i+f2i*H2i))/Hw
Bi = (f1i*H1i*T1i+f2i*H2i*T2i)/Hw

Oh, guess what are HT and H being passed to the correct function? In the 
onePhase solver version these are:

HT = H*T
H = H

with H being the fluid-structure heat transfer coeff and T the fluid
temperature. In the two-phase version:

HT = f1*H1*T1+f2*H2*T2
H = f1*H1+f2*H2

So this whole framework applies without issues to both the onePhase and 
twoPhase solvers. For clarity I will call the arguments passed to the correct
as sumH (=H) and sumHT (=HT). One thus always has:

Ai = (Hw+(sumHi))/Hw                                                        (5)
Bi = (sumHTi)/Hw

Now, again, what is the deal with the subscript i? Well, I remind you that it
simply means on which side of the HX wall I am evaluating the variable. Let us 
make an example. Let us assume we want to compute T on the primary side (thus 
Tp), and I doing this cell-by-cell. Let us assume we are in cell celli, which 
belongs to the primary. By accessing the fields sumHT[celli] and sumH[celli] I
am accessing what correponds to sumHTp and sumHp from the perspective of the
subscripting of equations 4 and 5. But so, how do we access the variables
sumHT and sumH in the cell 'opposite' to the primary, on the other side of the
wall, in the secondary cell zone? Well, that is what the mapping is for. In
particular, without going into the details of how the mapping was constructed (
which you can figure out from some comments here and there in the constructor),
the name extensions of the mapped variables in the code below (i.e. pOs or sOp)
stand for primary over secondary, secondary over primary respectively. What 
that means is that, if a cell cellip on the primary side corrsponds to a cell
cellis on the secondary side, accessing HTsOp[cellip] will yield the value
that HT has in cellis on the secondary, namely HT[cellis]. Conversely, 
accessing HTpOs[cellis] will yield HT[cellip]. Now, if it was just about an
indexing change, why used meshToMesh mapping? Well the point is that the the 
two meshes of the heat exchanger could be non-conformal, so that you loose a
clear primary-secondary correspondence between individual cells. Plus, even if
mesh were conformal, I see no point in not using already existing OpenFOAM
mechanics. Either way, this should help you understand to some degree the
implementation of the correct function. Lastly though, THX 
is the surface temperature of the heat exchanger, which is set by the correct.
In particular, THX = Tp in the primary and THX = Ts in the secondary. By 
considering this, as well as the form of the system of equations 4, it should
be easy to understand the double-loop implementation for setting THX.
*/

void Foam::heatExchanger::correct
( 
    const volScalarField& HT, 
    const volScalarField& H,
    volScalarField& THX
)
{
    tmp<scalarField> tHTpOs = mappingPtr_->mapSrcToTgt(HT.internalField());
    tmp<scalarField> tHTsOp = mappingPtr_->mapTgtToSrc(HT.internalField());
    tmp<scalarField> tHpOs = mappingPtr_->mapSrcToTgt(H.internalField());
    tmp<scalarField> tHsOp = mappingPtr_->mapTgtToSrc(H.internalField());
    const scalarField& HTpOs = tHTpOs.ref();
    const scalarField& HTsOp = tHTsOp.ref();
    const scalarField& HpOs = tHpOs.ref();
    const scalarField& HsOp = tHsOp.ref();

    //- The max against 1e-69 is just to avoid the retarded case in which the
    //  heat transfer coefficient on both side is 0 (which results in 
    //  Ap*As = 1.0), a case in which no heat is to be transferred via the HX
    //  but that would result in a division by 0 if not accounted for. Please
    //  note that in such a scenario THX is set to 0 yet it does not matter
    //  as said earlier, this only happens where the heat transfer coefficient
    //  0 or very close to it
    forAll(primaryCells_, i)
    {
        const label& celli(primaryCells_[i]);
        scalar Ap((Hw_+H[celli])/Hw_);
        scalar Bp(HT[celli]/Hw_);
        scalar As((Hw_+HsOp[celli])/Hw_);
        scalar Bs(HTsOp[celli]/Hw_);
        THX[celli] = (Bp*As+Bs)/max(Ap*As-1.0, 1e-69);
    }
    forAll(secondaryCells_, i)
    {
        const label& celli(secondaryCells_[i]);
        scalar Ap((Hw_+HpOs[celli])/Hw_);
        scalar Bp(HTpOs[celli]/Hw_);
        scalar As((Hw_+H[celli])/Hw_);
        scalar Bs(HT[celli]/Hw_);
        THX[celli] = (Bs*Ap+Bp)/max(As*Ap-1.0, 1e-69);
    }

    //- Kinda useless but whatever, at least the fields looks nicer in paraView
    //  if you set the boundaries
    THX.correctBoundaryConditions();
}

// ************************************************************************* //
