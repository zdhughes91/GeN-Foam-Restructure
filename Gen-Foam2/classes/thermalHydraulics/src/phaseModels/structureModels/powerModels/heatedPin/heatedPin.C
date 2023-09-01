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

#include "heatedPin.H"
#include "structure.H"
#include "addToRunTimeSelectionTable.H"
#include "SquareMatrix.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace powerModels
{
    defineTypeNameAndDebug(heatedPin, 0);
    addToRunTimeSelectionTable
    (
        powerModel, 
        heatedPin, 
        powerModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::powerModels::heatedPin::heatedPin
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
    Trad_
    (
        IOobject
        (
            "Trad."+typeName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_.cells().size()
    ),
    powerDensity_
    (
        IOobject
        (
            "powerDensity."+typeName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("powerDensity", dimPower/dimVol, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    Ti_
    (
        IOobject
        (
            "Ti."+typeName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimTemperature, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    To_
    (
        IOobject
        (
            "To."+typeName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimTemperature, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    Tav_
    (
        IOobject
        (
            "Tav."+typeName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimTemperature, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    Tmax_(0),
    Tmin_(1e69),
    meshSize_(0),
    r_(0),
    ri_(0),
    ro_(0),
    dr_(0),
    rhoCp_(0),
    k_(0),
    hollow_(0),
    cellToRegion_(mesh_.cells().size(), 0),
    pi_(constant::mathematical::pi),
    dA_(0)
{   
    structure_.setRegionField(*this, powerDensity_, "powerDensity");

    bool foundBoundaryTemperatures
    (
        Ti_.typeHeaderOk<volScalarField>(true)
    and To_.typeHeaderOk<volScalarField>(true)
    );

    typedef IOFieldField<Field, scalar> scalarFieldField;
    bool foundTrad
    (
        Trad_.typeHeaderOk<scalarFieldField>(true)
    );

    scalarList T0(0);

    forAll(this->toc(), regioni)
    {
        word region(this->toc()[regioni]);
        const dictionary& dict(this->subDict(region));
        
        //- Setup cellToRegion_ mapping
        const labelList& regionCells
        (
            structure_.cellLists()[region]
        );
        forAll(regionCells, i)
        {
            label celli(regionCells[i]);
            cellToRegion_[celli] = regioni;
        }

        //- Read region dict entries
        scalar ri(dict.get<scalar>("innerRadius"));
        scalar ro(dict.get<scalar>("outerRadius"));
        label meshSize(dict.get<label>("meshSize"));
        scalar dr((ro-ri)/(meshSize-1));
        scalar rhoCp(0);
        if (dict.found("rho") and dict.found("Cp"))
        {
            rhoCp =
                dict.get<scalar>("rho")*
                dict.get<scalar>("Cp");
        }
        else if (dict.found("rhoCp"))
        {
            rhoCp = dict.get<scalar>("rhoCp");
        }
        else
        {
            FatalErrorInFunction
                << "heatedPin region: " << region << " -> "
                << "specify either fuelRhoCp or both fuelRho and fuelCp"
                << exit(FatalError);
        }
        
        scalar k(dict.get<scalar>("k"));
        bool hollow((ri >= 1e-5) ? true : false);

        if (!foundBoundaryTemperatures and !foundTrad)
        {
            T0.append(dict.get<scalar>("T"));
        }

        //- Calc mesh array
        scalarList r(0);
        r.append(ri);
        for (int i = 0; i < meshSize-1; i++)
        {
            r.append(r.last() + dr);
        }

        //- Calc ring areas
        scalarList dA(0);
        dA.append(pi_*(sqr(r[0]+dr/2.0)-sqr(r[0])));
        for (int i = 1; i < meshSize-1; i++)
        {
            dA.append(pi_*(sqr(r[i]+dr/2.0)-sqr(r[i]-dr/2.0)));
        }
        dA.append(pi_*(sqr(r[meshSize-1])-sqr(r[meshSize-1]-dr/2.0)));

        //- Fill in lists for this region
        meshSize_.append(meshSize);
        r_.append(r);
        ri_.append(ri);
        ro_.append(ro);
        dr_.append(dr);
        rhoCp_.append(rhoCp);
        k_.append(k);
        hollow_.append(hollow);
        dA_.append(dA);
    }

    //- If Trad not found, init it from either the boundary temperatures
    //  (I mean boundary in a mathematical sense, i.e. inner/outer fuel/clad
    //  temperature) or from dictionary values Tf0, Tc0 read previously
    if (!foundTrad)
    {
        forAll(mesh_.cells(), i)
        {
            Trad_.set(i, new Field<scalar>(0, 0));
        }

        //- If the files are present, reconstruct initial Trad_ profile 
        //  analytically. The analytical form is:
        //  T(r) = -(1/4)*powerDensity_(r)*r^2/k + C*ln(r)/k + D
        //  with C and D coming from imposing fixedValue BC on all sides,
        //  equal to the starting temperatures found in the files.
        if (foundBoundaryTemperatures)
        {
            Info<< "Found heatedPin temperatures "
                << "reconstructing profiles " << endl;
            forAll(this->cellList_, i)
            {
                label celli(this->cellList_[i]);

                label regioni(cellToRegion_[celli]);
                const scalarList& rList(r_[regioni]);
                scalar k(k_[regioni]);
                scalar ri(ri_[regioni]);
                scalar ro(ro_[regioni]);
                bool hollow(hollow_[regioni]);

                Trad_.set(celli, new Field<scalar>(meshSize_[regioni], 0));
                scalar q = powerDensity_[celli];
                scalar ti = Ti_[celli];
                scalar to = To_[celli];
                scalar C;
                scalar D;
                
                if (hollow)
                {
                    C = 
                        (k*(ti-to)-0.25*q*(sqr(ro)-sqr(ri)))/
                        log(ri/ro);
                }
                else
                {
                    C = 0.0;
                }
                D = to+(0.25*q*sqr(ro)-C*log(ro))/k;
                
                forAll(Trad_[celli], j)
                {
                    scalar r(rList[j]);
                    Trad_[celli][j] = -0.25*q*sqr(r)/k + D;
                    Trad_[celli][j] +=
                        (hollow) ? C*log(r)/k : 0.0;
                }
            }
        }
        else //- Otherwise, read from dict
        {
            Info<< "Reading heatedPin initial temperatures from "
                << "dictionary" << endl;

            forAll(this->cellList_, i)
            {
                label celli(this->cellList_[i]);
                label regioni(cellToRegion_[celli]);
                Trad_.set(celli, new Field<scalar>(meshSize_[regioni], 0));
                forAll(Trad_[celli], subCelli)
                {
                    Trad_[celli][subCelli] = T0[regioni];
                }
            }
        }
    }
    else
    {
        Info<< "Setting heatedPin initial temperatures from "
                << Trad_.name() << endl;
    }
    
    //- Set I/O fields and compute initial scalar max, min
    scalar Tavav(0);
    scalar totV(0.0);
    const scalarList& V(mesh_.V());
    forAll(this->cellList_, i)
    {
        label celli(this->cellList_[i]);
        label regioni(cellToRegion_[celli]);
        const label& n(meshSize_[regioni]);
        const scalarField& Trad(Trad_[celli]);
        Ti_[celli] = Trad[0];
        To_[celli] = Trad[n-1];

        const scalarList& rRegion(r_[regioni]);
        scalar& Tavi(Tav_[celli]);
        
        updateLocalAvgGlobalMinMaxT
        (
            0,
            n,
            rRegion,
            dr_[regioni],
            Trad,
            Tavi,
            Tmin_,
            Tmax_
        );

        //- This is for updating the global averages, not the local cell ones!
        const scalar& dV(V[celli]);
        totV += dV;
        Tavav += Tavi*dV;
    }

    //- Sync across processors
    reduce(totV, sumOp<scalar>());
    reduce(Tavav, sumOp<scalar>());
    reduce(Tmax_, maxOp<scalar>());
    reduce(Tmin_, minOp<scalar>());

    Tavav /= totV;

    //- Initialize in dict
    this->IOdictionary::set("Tavav", Tavav);
    this->IOdictionary::set("Tmax", Tmax_);
    this->IOdictionary::set("Tmin", Tmin_);

    Ti_.correctBoundaryConditions();
    To_.correctBoundaryConditions();

    //- Finally, set up interfacial area
    this->setInterfacialArea();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::powerModels::heatedPin::~heatedPin()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::powerModels::heatedPin::setInterfacialArea()
{
    forAll(this->cellList_, i)
    {
        const label& celli(this->cellList_[i]);
        iA_[celli] = 2.0*alpha_[celli]/(ro_[cellToRegion_[celli]]);
    }
    iA_.correctBoundaryConditions();
}

void Foam::powerModels::heatedPin::updateLocalAvgGlobalMinMaxT
(
    const label& starti,
    const label& endi,
    const scalarList& r,
    const scalar& dr,
    const scalarField& Trad,
    scalar& Tavi,
    scalar& Tmin,
    scalar& Tmax
)
{
    scalar intr(0);
    scalar intTr(0);
    for(int j = starti; j < endi; j++)
    {   
        const scalar& T(Trad[j]);
        scalar rdr(r[j]*dr);
        
        //- Cells at the mesh ends are only half as wide (the other half
        //  belongs to the ghost node). Thus, weigh temperatures at the extrema
        //  by a factor 0.5
        if (j == starti or j == endi-1)
        {
            intr += rdr/2.0;
            intTr += T*rdr/2.0;
        }
        else
        {
            intr += rdr;
            intTr += T*rdr;
        }
        if (T > Tmax) Tmax = T;
        if (T < Tmin) Tmin = T;
    }
    Tavi = intTr/intr;
}

void 
Foam::powerModels::heatedPin::updateLocalTemperatureProfile
(
    const label& celli,
    const scalar& HTSumi,
    const scalar& HSumi
)
{
    //-
    scalarField& Trad(Trad_[celli]);

    //- Read region values
    const label& regioni(cellToRegion_[celli]);
    const label& meshSize(meshSize_[regioni]);
    const scalarList& rRegion(r_[regioni]);
    const scalarList& dARegion(dA_[regioni]);
    const scalar& dr(dr_[regioni]);
    const scalar& k(k_[regioni]);
    const scalar& q(powerDensity_[celli]);
    const scalarField& TOld = Trad_.oldTime()[celli];

    //- Recurrent quantities
    scalar drh(dr/2.0);
    scalar dt(mesh_.time().deltaT().value());
    scalar twoPkByDr(2.0*pi_*k/dr);
    scalar X(rhoCp_[regioni]/dt);

    //- Init matrix, source
    SquareMatrix<scalar> M(meshSize, meshSize, Foam::zero());
    List<scalar> S(meshSize, 0.0);
    
    //- Fill in matrix, source coeffs
    {
        //- Set zeroGradient BC at pin inner surface
        {
            const scalar& r(rRegion[0]);
            const scalar& dA(dARegion[0]);
            scalar B(twoPkByDr*(r+drh));
            scalar XdA(X*dA);
            M[0][1] =   -B;
            M[0][0] =   B+XdA;
            S[0] =      q*dA+TOld[0]*XdA;
        }
        
        //- Set bulk
        for (int i = 1; i < meshSize-1; i++)
        {
            const scalar& r(rRegion[i]);
            const scalar& dA(dARegion[i]);
            scalar B(twoPkByDr*(r+drh));
            scalar C(twoPkByDr*(r-drh));
            scalar XdA(X*dA);
            M[i][i+1] =     -B;
            M[i][i-1] =     -C;
            M[i][i] =       B+C+XdA;
            S[i] =          q*dA+TOld[i]*XdA;
        }

        //- Set convective BC at pin outer surface
        {
            label i(meshSize-1);
            const scalar& r(rRegion[i]);
            const scalar& dA(dARegion[i]);
            scalar C(twoPkByDr*(r-drh));
            scalar D(2.0*pi_*r);
            scalar XdA(X*dA);
            M[i][i-1] =     -C;
            M[i][i] =       C+D*HSumi+XdA;
            S[i] =          q*dA+TOld[i]*XdA + D*HTSumi;
        }
    }

    //- Solve linear system
    solve(Trad, M, S);

    //- Set inner/outer pin temperature fields
    Ti_[celli] = Trad[0];
    To_[celli] = Trad[meshSize-1];

    //- Check energy conservation via linear power
    /*
    scalar aLP(q*pi_*(sqr(rRegion[0])-sqr(rRegion[meshSize-1])));
    scalar nLP((HSumi*To_[celli]-HTSumi)*2.0*pi_*rRegion[meshSize-1]);
    Info<< celli << " " << nLP << " " << aLP << " W/m" << endl;
    */

    //- Update local T averages and local min/max
    scalar& Tavi(Tav_[celli]);
    updateLocalAvgGlobalMinMaxT
    (
        0,
        meshSize,
        rRegion,
        dr,
        Trad,
        Tavi,
        Tmin_,
        Tmax_
    );
}


void Foam::powerModels::heatedPin::correct
(
    const volScalarField& HTSum,  // == SUM_j [htc_j*T_j*frac_j]
    const volScalarField& HSum    // == SUM_j [htc_j*frac_j]
)
{
    //- Reset min, max, fuel, clad temperatures
    Tmax_ = 0.0;
    Tmin_ = 1e69;
    
    //- Update temperatures cell-by-cell and compute averages over the entire
    //  spatial extent of the heatedPin model (what I call global 
    //  averages, opposed to local averages, which are the average temperature
    //  values, fuel and clad, of the local cell radial pin temperature 
    //  profile)
    const scalarField& V(mesh_.V());
    scalar totV(0);
    scalar Tavav(0);
    forAll(this->cellList_, i)
    {
        label celli(this->cellList_[i]);
        updateLocalTemperatureProfile(celli, HTSum[celli], HSum[celli]);
        const scalar& dV(V[celli]);
        totV += dV;
        Tavav += Tav_[celli]*dV;
    }
    reduce(totV, sumOp<scalar>());
    reduce(Tavav, sumOp<scalar>());
    Tavav /= totV;

    reduce(Tmax_, maxOp<scalar>());
    reduce(Tmin_, minOp<scalar>());

    Info<< "T.heatedPin (avg min max) = " 
        << Tavav << " " << Tmin_ << " " << Tmax_ << " K" << endl;

    //- Save these to the dictionary
    this->IOdictionary::set("Tavav", Tavav);
    this->IOdictionary::set("Tmax", Tmax_);
    this->IOdictionary::set("Tmin", Tmin_);
}


void Foam::powerModels::heatedPin::correctT(volScalarField& T) const
{
    //- Set T to pin surface temperature, i.e. Tco_
    forAll(cellList_, i)
    {
        label celli(cellList_[i]);
        T[celli] = To_[celli];
    }
}


void Foam::powerModels::heatedPin::powerOff()
{
    powerDensity_ *= 0.0;
}

// ************************************************************************* //