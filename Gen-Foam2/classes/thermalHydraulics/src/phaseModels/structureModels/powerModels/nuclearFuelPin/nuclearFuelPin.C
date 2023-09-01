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

#include "nuclearFuelPin.H"
#include "structure.H"
#include "addToRunTimeSelectionTable.H"
#include "SquareMatrix.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace powerModels
{
    defineTypeNameAndDebug(nuclearFuelPin, 0);
    addToRunTimeSelectionTable
    (
        powerModel, 
        nuclearFuelPin, 
        powerModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::powerModels::nuclearFuelPin::nuclearFuelPin
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
    Tfi_
    (
        IOobject
        (
            "Tfi."+typeName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimTemperature, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    Tfo_
    (
        IOobject
        (
            "Tfo."+typeName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimTemperature, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    Tci_
    (
        IOobject
        (
            "Tci."+typeName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimTemperature, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    Tco_
    (
        IOobject
        (
            "Tco."+typeName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimTemperature, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    Tfav_
    (
        IOobject
        (
            "Tfav."+typeName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimTemperature, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    Tcav_
    (
        IOobject
        (
            "Tcav."+typeName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimTemperature, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    Tfmax_(0),
    Tfmin_(1e69),
    Tcmax_(0),
    Tcmin_(1e69),
    fractionOfPowerFromNeutronics_(0),
    fuelMeshSize_(0),//this->get<label>("fuelMeshSize")),
    cladMeshSize_(0),//this->get<label>("cladMeshSize")),
    meshSize_(0),
    r_(0),
    rfi_(0),//this->get<scalar>("fuelInnerRadius")),
    rfo_(0),//this->get<scalar>("fuelOuterRadius")),
    rci_(0),//this->get<scalar>("cladInnerRadius")),
    rco_(0),//this->get<scalar>("cladOuterRadius")),
    drf_(0),//(rfo_-rfi_)/(fuelSubMeshSize_-1)),
    drc_(0),//(rco_-rci_)/(cladSubMeshSize_-1)),
    drg_(0),//(rci_-rfo_),
    rhoCpf_(0),//this->get<scalar>("fuelDensity")),
    rhoCpc_(0),//this->get<scalar>("cladDensity")),
    kf_(0),//this->get<scalar>("fuelConductivity")),
    kc_(0),//this->get<scalar>("cladConductivity")),
    gapH_(0),//this->get<scalar>("gapConductance")),
    hollowFuel_(0),// (rfi_ >= 1e-5) ? true : false),
    cellToRegion_(mesh_.cells().size(), 0),
    regionIndexToRegionName_(0),
    pi_(constant::mathematical::pi),
    dA_(0),
    gapHPowerDensityTable_(0),
    useGapHPowerDensityTable_(0)
{   
    structure_.setRegionField(*this, structureRef.powerDensityNeutronics(), "powerDensity");

    bool foundBoundaryTemperatures
    (
        Tfi_.typeHeaderOk<volScalarField>(true)
    and Tfo_.typeHeaderOk<volScalarField>(true)
    and Tci_.typeHeaderOk<volScalarField>(true)
    and Tco_.typeHeaderOk<volScalarField>(true)
    );

    typedef IOFieldField<Field, scalar> scalarFieldField;
    bool foundTrad
    (
        Trad_.typeHeaderOk<scalarFieldField>(true)
    );

    scalarList Tf0(0);
    scalarList Tc0(0);

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

        //- Add to regionIndexToRegionName_ mapping
        regionIndexToRegionName_.append(region);

        //- Read region dict entries
        scalar fractionOfPowerFromNeutronics(dict.lookupOrDefault<scalar>("fractionOfPowerFromNeutronics",1.0));        
        scalar rfi(dict.get<scalar>("fuelInnerRadius"));
        scalar rfo(dict.get<scalar>("fuelOuterRadius"));
        scalar rci(dict.get<scalar>("cladInnerRadius"));
        scalar rco(dict.get<scalar>("cladOuterRadius"));
        label fuelMeshSize(dict.get<label>("fuelMeshSize"));
        label cladMeshSize(dict.get<label>("cladMeshSize"));
        label meshSize(fuelMeshSize+cladMeshSize);
        scalar drf((rfo-rfi)/(fuelMeshSize-1));
        scalar drc((rco-rci)/(cladMeshSize-1));
        scalar drg(rci-rfo);
        scalar rhoCpf(0);
        if (dict.found("fuelRho") and dict.found("fuelCp"))
        {
            rhoCpf =
                dict.get<scalar>("fuelRho")*
                dict.get<scalar>("fuelCp");
        }
        else if (dict.found("fuelRhoCp"))
        {
            rhoCpf = dict.get<scalar>("fuelRhoCp");
        }
        else
        {
            FatalErrorInFunction
                << "nuclearFuelPin region: " << region << " -> "
                << "specify either fuelRhoCp or both fuelRho and fuelCp"
                << exit(FatalError);
        }
        scalar rhoCpc(0);
        if (dict.found("cladRho") and dict.found("cladCp"))
        {
            rhoCpc =
                dict.get<scalar>("cladRho")*
                dict.get<scalar>("cladCp");
        }
        else if (dict.found("cladRhoCp"))
        {
            rhoCpc = dict.get<scalar>("cladRhoCp");
        }
        else
        {
            FatalErrorInFunction
                << "nuclearFuelPin region: " << region << " -> "
                << "specify either cladRhoCp or both cladRho and cladCp"
                << exit(FatalError);
        }
        scalar kf(dict.get<scalar>("fuelK"));
        scalar kc(dict.get<scalar>("cladK"));
        
        bool hollowFuel((rfi >= 1e-5) ? true : false);

        if (!foundBoundaryTemperatures and !foundTrad)
        {
            Tf0.append(dict.get<scalar>("fuelT"));
            Tc0.append(dict.get<scalar>("cladT"));
        }

        //- Calc mesh array
        scalarList r(0);
        r.append(rfi);
        for (int i = 0; i < fuelMeshSize-1; i++)
        {
            r.append(r.last() + drf);
        }
        r.append(r.last() + drg);
        for (int i = 0; i < cladMeshSize-1; i++)
        {
            r.append(r.last() + drc);
        }

        //- Calc dA
        scalarList dA(0);
        dA.append(pi_*(sqr(r[0]+drf/2.0)-sqr(r[0])));
        for(int i = 1; i < fuelMeshSize-1; i++)
        {
            dA.append(pi_*(sqr(r[i]+drf/2.0)-sqr(r[i]-drf/2.0)));
        }
        dA.append
        (
            pi_*(sqr(r[fuelMeshSize-1])-sqr(r[fuelMeshSize-1]-drf/2.0))
        );
        dA.append
        (
            pi_*(sqr(r[fuelMeshSize]+drf/2.0)-sqr(r[fuelMeshSize]))
        );
        for(int i = 1; i < cladMeshSize-1; i++)
        {
            dA.append(pi_*(sqr(r[i]+drc/2.0)-sqr(r[i]-drc/2.0)));
        }
        dA.append
        (
            pi_*(sqr(r[meshSize-1])-sqr(r[meshSize-1]-drc/2.0))
        );

        //- Fill in lists for this region
        fractionOfPowerFromNeutronics_.append(fractionOfPowerFromNeutronics), 
        fuelMeshSize_.append(fuelMeshSize);
        cladMeshSize_.append(cladMeshSize);
        meshSize_.append(meshSize);
        r_.append(r);
        rfi_.append(rfi);
        rfo_.append(rfo);
        rci_.append(rci);
        rco_.append(rco);
        drf_.append(drf);
        drc_.append(drc);
        drg_.append(drg);
        rhoCpf_.append(rhoCpf);
        rhoCpc_.append(rhoCpc);
        kf_.append(kf);
        kc_.append(kc);
        
        hollowFuel_.append(hollowFuel);
        dA_.append(dA);

        //- Construct gapHPowerDensityTable if found, otherwise use the
        //  provided constant value
        scalar gapH(0);
        word tableName("gapHPowerDensity");
        bool foundTable(dict.found(tableName));
        bool foundValue(dict.found("gapH"));
        if (foundTable)
        {
            /*
                The table Function1 requires a particular input. In general,
                it should be something like this

                topLevelDictName
                {
                    type                table;
                    topLevelDictName    table
                    (
                        (0 0)
                        (1 1)
                        (...)
                    );
                }

                however, I only want the user to give an input in the form

                powerModel
                {
                    type        nuclearFuelPin;

                    ...

                    gapHPowerDensity table
                    (
                        (0 0)
                        (1 1)
                        (...)
                    );
                }

                in which powerModel is the topLevelDictName.

                The code below does just this, by creating a copy dict of
                powerModel renaming it to gapHPowerDensity table, resetting
                type to table, and passing that to the Function1 table 
                selector
            */
            dictionary tableDict(tableName);
            tableDict.merge(dict);
            tableDict.set("type", "table");
            gapHPowerDensityTable_.insert
            (
                region,
                autoPtr<Function1<scalar>>
                (
                    Function1<scalar>::New
                    (
                        tableName,
                        tableDict,
                        "table"
                    )
                )
            );
        }
        if (foundValue)
        {
            gapH = dict.get<scalar>("gapH");
        }

        //- The gapH list needs to have the same length as the number of 
        //  regions no matter what, or the indexing will stop working as 
        //  intended
        gapH_.append(gapH);

        useGapHPowerDensityTable_.append(foundTable);

        if (foundTable and foundValue)
        {
            FatalErrorInFunction
                << "nuclearFuelPin region: " << region << " -> "
                << "provide either a gapH value or a gapHPowerDensityTable "
                << "but not both!"
                << exit(FatalError);
        }
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
            Info<< "Found nuclearFuelPin temperatures "
                << "reconstructing profiles " << endl;
            forAll(this->cellList_, i)
            {
                label celli(this->cellList_[i]);

                label regioni(cellToRegion_[celli]);
                const scalarList& rList(r_[regioni]);
                scalar kf(kf_[regioni]);
                scalar kc(kc_[regioni]);
                scalar rfi(rfi_[regioni]);
                scalar rfo(rfo_[regioni]);
                scalar rci(rci_[regioni]);
                scalar rco(rco_[regioni]);
                scalar fractionOfPowerFromNeutronics(fractionOfPowerFromNeutronics_[regioni]);
                bool hollowFuel(hollowFuel_[regioni]);

                Trad_.set(celli, new Field<scalar>(meshSize_[regioni], 0));
                scalar q = structure_.powerDensityNeutronics()[celli] * fractionOfPowerFromNeutronics;
                scalar tfi = Tfi_[celli];
                scalar tfo = Tfo_[celli];
                scalar tci = Tci_[celli];
                scalar tco = Tco_[celli];
                scalar Cf;
                scalar Df;
                scalar Cc;
                scalar Dc;
                
                if (hollowFuel)
                {
                    Cf = 
                        (kf*(tfi-tfo)-0.25*q*(sqr(rfo)-sqr(rfi)))/
                        log(rfi/rfo);
                }
                else
                {
                    Cf = 0.0;
                }
                Df = tfo+(0.25*q*sqr(rfo)-Cf*log(rfo))/kf;
                Cc = (tco-tci)*kc/(log(rco/rci));
                Dc = tci - log(rci)*Cc/kc;
                
                forAll(Trad_[celli], j)
                {
                    scalar r(rList[j]);
                    if (j < fuelMeshSize_[regioni])
                    {
                        Trad_[celli][j] = -0.25*q*sqr(r)/kf + Df;
                        Trad_[celli][j] +=
                            (hollowFuel) ? Cf*log(r)/kf : 0.0;
                    }
                    else
                    {   
                        Trad_[celli][j] = Cc*log(r)/kc + Dc;
                    }
                }
           }
        }
        else //- Otherwise, read from dict
        {
            Info<< "Reading nuclearFuelPin initial temperatures from "
                << "dictionary" << endl;

            forAll(this->cellList_, i)
            {
                label celli(this->cellList_[i]);
                label regioni(cellToRegion_[celli]);
                Trad_.set(celli, new Field<scalar>(meshSize_[regioni], 0));
                forAll(Trad_[celli], subCelli)
                {
                    Trad_[celli][subCelli] = 
                        (subCelli < fuelMeshSize_[regioni]) ?
                        Tf0[regioni] : Tc0[regioni];
                }
            }
        }
    }
    else
    {
        Info<< "Setting nuclearFuelPin initial temperatures from "
                << Trad_.name() << endl;
    }
    
    //- Set I/O fields and compute initial scalar max, min
    scalar Tfavav(0);
    scalar Tcavav(1e69);
    scalar totV(0.0);
    const scalarList& V(mesh_.V());
    forAll(this->cellList_, i)
    {
        label celli(this->cellList_[i]);
        label regioni(cellToRegion_[celli]);
        const label& nf(fuelMeshSize_[regioni]);
        const label& n(meshSize_[regioni]);
        const scalarField& Trad(Trad_[celli]);
        Tfi_[celli] = Trad[0];
        Tfo_[celli] = Trad[nf-1];
        Tci_[celli] = Trad[nf];
        Tco_[celli] = Trad[n-1];

        const scalarList& rRegion(r_[regioni]);
        scalar& Tfavi(Tfav_[celli]);
        scalar& Tcavi(Tcav_[celli]);
        
        updateLocalAvgGlobalMinMaxT
        (
            0,
            nf,
            rRegion,
            drf_[regioni],
            Trad,
            Tfavi,
            Tfmin_,
            Tfmax_
        );
        updateLocalAvgGlobalMinMaxT
        (
            nf,
            n,
            rRegion,
            drc_[regioni],
            Trad,
            Tcavi,
            Tcmin_,
            Tcmax_
        );

        // Update average fuel and clad temp used for coupling
        this->structureRef().TFuelAv()[celli] = Tfav_[celli];
        this->structureRef().TCladAv()[celli] = Tcav_[celli];

        //- This is for updating the global averages, not the local cell ones!
        const scalar& dV(V[celli]);
        totV += dV;
        Tfavav += Tfavi*dV;
        Tcavav += Tcavi*dV;
    }

    //- Sync across processors
    reduce(totV, sumOp<scalar>());
    reduce(Tfavav, sumOp<scalar>());
    reduce(Tcavav, sumOp<scalar>());
    reduce(Tfmax_, maxOp<scalar>());
    reduce(Tfmin_, minOp<scalar>());
    reduce(Tcmax_, maxOp<scalar>());
    reduce(Tcmin_, minOp<scalar>());

    Tfavav /= totV;
    Tcavav /= totV;

    //- Initialize in dict
    this->IOdictionary::set("Tfavav", Tfavav);
    this->IOdictionary::set("Tcavav", Tcavav);
    this->IOdictionary::set("Tfmax", Tfmax_);
    this->IOdictionary::set("Tfmin", Tfmin_);
    this->IOdictionary::set("Tcmax", Tcmax_);
    this->IOdictionary::set("Tcmin", Tcmin_);

    Tfi_.correctBoundaryConditions();
    Tfo_.correctBoundaryConditions();
    Tci_.correctBoundaryConditions();
    Tco_.correctBoundaryConditions();

    //- Finally, set up interfacial area
    this->setInterfacialArea();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::powerModels::nuclearFuelPin::~nuclearFuelPin()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::powerModels::nuclearFuelPin::setInterfacialArea()
{
    forAll(this->cellList_, i)
    {
        const label& celli(this->cellList_[i]);
        iA_[celli] = 2.0*alpha_[celli]/(rco_[cellToRegion_[celli]]);
    }
    iA_.correctBoundaryConditions();
}

void Foam::powerModels::nuclearFuelPin::updateLocalAvgGlobalMinMaxT
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
Foam::powerModels::nuclearFuelPin::updateLocalTemperatureProfile
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
    const word& region(regionIndexToRegionName_[regioni]);
    const label& fuelMeshSize(fuelMeshSize_[regioni]);
    const label& meshSize(meshSize_[regioni]);
    const scalarList& rRegion(r_[regioni]);
    const scalarList& dARegion(dA_[regioni]);
    const scalar& drf(drf_[regioni]);
    const scalar& drc(drc_[regioni]);
    const scalar& kf(kf_[regioni]);
    const scalar& kc(kc_[regioni]);
    const scalar& rfo(rfo_[regioni]);
    const scalar& rci(rci_[regioni]);
    const scalar& fractionOfPowerFromNeutronics(fractionOfPowerFromNeutronics_[regioni]);
    
    const scalarField& TOld = Trad_.oldTime()[celli];

    //- Update power density
    const scalar& qRef(structure_.powerDensityNeutronics()[celli]);
    scalar q = qRef * fractionOfPowerFromNeutronics;

    scalar gapH
    (
        (useGapHPowerDensityTable_[regioni]) ?
        gapHPowerDensityTable_[region]->value(q) :
        gapH_[regioni]
    );

    //- Init matrix, source
    SquareMatrix<scalar> M(meshSize, meshSize, Foam::zero());
    List<scalar> S(meshSize, 0.0);

    //- Recurrent quantities
    scalar dt(mesh_.time().deltaT().value());
    scalar Xf(rhoCpf_[regioni]/dt);
    scalar Xc(rhoCpc_[regioni]/dt);
    scalar twoPkByDrf(2.0*pi_*kf/drf);
    scalar twoPkByDrc(2.0*pi_*kc/drc);
    scalar drhf(drf/2.0);
    scalar drhc(drc/2.0);
    
    //- Construct matrix, source
    {
        //- Set zeroGradient BC at fuel inner surface
        {
            const scalar& r(rRegion[0]);
            const scalar& dA(dARegion[0]);
            scalar B(twoPkByDrf*(r+drhf));
            scalar XdA(Xf*dA);
            M[0][1] =   -B;
            M[0][0] =   B+XdA;
            S[0] =      q*dA+TOld[0]*XdA;
        }

        //- Fuel bulk
        for (int i = 1; i < fuelMeshSize-1; i++)
        {
            const scalar& r(rRegion[i]);
            const scalar& dA(dARegion[i]);
            scalar B(twoPkByDrf*(r+drhf));
            scalar C(twoPkByDrf*(r-drhf));
            scalar XdA(Xf*dA);
            M[i][i+1] =     -B;
            M[i][i-1] =     -C;
            M[i][i] =       B+C+XdA;
            S[i] =          q*dA+TOld[i]*XdA;
        }

        //- Fuel outer surface, convective BC with inner cladding surface via
        //  gap conductance
        {
            label i(fuelMeshSize-1);
            const scalar& r(rRegion[i]);
            const scalar& dA(dARegion[i]);
            scalar C(twoPkByDrf*(r-drhf));
            scalar D(2.0*pi_*r*gapH);
            scalar XdA(Xf*dA);
            M[i][i+1] =     -D;
            M[i][i-1] =     -C;
            M[i][i] =       C+D+XdA;
            S[i] =          q*dA+TOld[i]*XdA;
        }

        //- Cladding inner surface, convective BC with outer fuel surface via
        //  gap conductance adjusted by radii ratio (to preserve total heat
        //  flow as geometry is cylindrical)
        {
            label i(fuelMeshSize);
            const scalar& r(rRegion[i]);
            const scalar& dA(dARegion[i]);
            scalar C(twoPkByDrc*(r+drhc));
            scalar D(2.0*pi_*r*(rfo/rci)*gapH);
            scalar XdA(Xc*dA);
            M[i][i+1] =     -C;
            M[i][i-1] =     -D;
            M[i][i] =       C+D+XdA;
            S[i] =          TOld[i]*XdA;
        }

        //- Cladding bulk
        for (int i = fuelMeshSize+1; i < meshSize-1; i++)
        {
            const scalar& r(rRegion[i]);
            const scalar& dA(dARegion[i]);
            scalar B(twoPkByDrc*(r+drhc));
            scalar C(twoPkByDrc*(r-drhc));
            scalar XdA(Xc*dA);
            M[i][i+1] =     -B;
            M[i][i-1] =     -C;
            M[i][i] =       B+C+XdA;
            S[i] =          TOld[i]*XdA;
        }

        //- Cladding outer surface, convective BC with fluid(s) wetting the pin
        {
            label i(meshSize-1);
            const scalar& r(rRegion[i]);
            const scalar& dA(dARegion[i]);
            scalar C(twoPkByDrc*(r-drhc));
            scalar D(2.0*pi_*r);
            scalar XdA(Xc*dA);
            M[i][i-1] =     -C;
            M[i][i] =       C+D*HSumi+XdA;
            S[i] =          TOld[i]*XdA + D*HTSumi;
        }
    }

    //- Solve linear system
    solve(Trad, M, S);

    //- Set fields (inner/outer fuel/clad)
    Tfi_[celli] = Trad[0];
    Tfo_[celli] = Trad[fuelMeshSize-1];
    Tci_[celli] = Trad[fuelMeshSize];
    Tco_[celli] = Trad[meshSize-1];

    //- Update local T averages and local min/max
    scalar& Tfavi(Tfav_[celli]);
    scalar& Tcavi(Tcav_[celli]);
    updateLocalAvgGlobalMinMaxT
    (
        0,
        fuelMeshSize,
        rRegion,
        drf,
        Trad,
        Tfavi,
        Tfmin_,
        Tfmax_
    );
    updateLocalAvgGlobalMinMaxT
    (
        fuelMeshSize,
        meshSize,
        rRegion,
        drc,
        Trad,
        Tcavi,
        Tcmin_,
        Tcmax_
    );

    // Update average fuel and clad temp used for coupling
    this->structureRef().TFuelAv()[celli] = Tfav_[celli];
    this->structureRef().TCladAv()[celli] = Tcav_[celli];

    /*
    //- Check energy conservation via linear power comparison (analytic
    //  vs numerical at cladding surface)
    const scalar& rfi(rfi_[regioni]);
    const scalar& rco(rco_[regioni]);
    scalar analyticalLP
    (
        q*pi_*(sqr(rfo)-sqr(rfi))
    );
    scalar numericalLP
    (
        (HSumi*Tco_[celli]-HTSumi)*2.0*pi_*rco
    );
    Info<< celli << " " << numericalLP << " " << analyticalLP << " W/m" 
        << endl;

    //- Check energy conservation via heat flux comparison. However, since
    //  you can't directly compare heat fluxes (it's total heat that is
    //  conserved, not heat fluxes), the heat flux trhough the cladding is
    //  adjusted to take into consideration the difference in surface areas
    //  between outer cladding and outer fuel
    scalar heatFluxg(gapH*(Tfo_[celli]-Tci_[celli]));
    scalar heatFluxAdjc((HSumi*Tco_[celli]-HTSumi)*(rco/rfo));
    Info<< celli << " " << heatFluxg << " " << heatFluxAdjc << " W/m2" 
        << endl;
    */
}


void Foam::powerModels::nuclearFuelPin::correct
(
    const volScalarField& HTSum,  // == SUM_j [htc_j*T_j*frac_j]
    const volScalarField& HSum    // == SUM_j [htc_j*frac_j]
)
{
    //- Reset min, max, fuel, clad temperatures
    Tfmax_ = 0.0;
    Tfmin_ = 1e69;
    Tcmax_ = 0.0;
    Tcmin_ = 1e69;
    
    //- Update temperatures cell-by-cell and compute averages over the entire
    //  spatial extent of the nuclearFuelPin model (what I call global 
    //  averages, opposed to local averages, which are the average temperature
    //  values, fuel and clad, of the local radial pin temperature profile)
    const scalarField& V(mesh_.V());
    scalar totV(0);
    scalar Tfavav(0);
    scalar Tcavav(0);
    forAll(this->cellList_, i)
    {
        label celli(this->cellList_[i]);
        updateLocalTemperatureProfile(celli, HTSum[celli], HSum[celli]);
        const scalar& dV(V[celli]);
        totV += dV;
        Tfavav += Tfav_[celli]*dV;
        Tcavav += Tcav_[celli]*dV;
    }
    reduce(totV, sumOp<scalar>());
    reduce(Tfavav, sumOp<scalar>());
    reduce(Tcavav, sumOp<scalar>());
    Tfavav /= totV;
    Tcavav /= totV;

    reduce(Tfmax_, maxOp<scalar>());
    reduce(Tfmin_, minOp<scalar>());
    reduce(Tcmax_, maxOp<scalar>());
    reduce(Tcmin_, minOp<scalar>());

    Info<< "T.nuclearFuelPin.fuel (avg min max) = " 
        << Tfavav << " " << Tfmin_ << " " << Tfmax_ << " K" << endl;
    Info<< "T.nuclearFuelPin.clad (avg min max) = " 
        << Tcavav << " " << Tcmin_ << " " << Tcmax_ << " K" << endl;

    //- Save these to the dictionary
    this->IOdictionary::set("Tfavav", Tfavav);
    this->IOdictionary::set("Tcavav", Tcavav);
    this->IOdictionary::set("Tfmax", Tfmax_);
    this->IOdictionary::set("Tfmin", Tfmin_);
    this->IOdictionary::set("Tcmax", Tcmax_);
    this->IOdictionary::set("Tcmin", Tcmin_);
}


void Foam::powerModels::nuclearFuelPin::correctT(volScalarField& T) const
{
    //- Set T to pin surface temperature, i.e. Tco_
    forAll(cellList_, i)
    {
        label celli(cellList_[i]);
        T[celli] = Tco_[celli];
    }
}


// ************************************************************************* //


/*-------------------------- EQUATION DERIVATION ----------------------------*\

//- WARNING: THIS IS THE DERIVATION OF THE PREVIOUS FINITE-DIFFERENCE EQUATIONS
//  WHICH ARE NO LONGER USED. THE CURRENT CODE IS BASED ON FINITE VOLUME SCHEME
//  WHICH IS SIGNIFICANTLY MORE CONSERVATIVE THAN THE FINITE-DIFFERENCE ONE
//  FOR THE SAME MESH SIZE. THE FINITE VOLUME DERIVATION IS NOT REPORTED AS...
//  I SIMPLY DON'T HAVE TIME NOW

I added this section in the hope that it might spare quite some time to anyone
who will unfortunately have to deal with this code after me. It shows how the
1-D matrices in updateLocalTemperatureProfile are built.

Let us start with the heat diffusion equation in cylindrical coordinates with
temporally constant cp, rho and spatially constant k. For simplicity,
let us not consider a pin right now, just a regular uniform hollow cylinder:

rho*cp*ddtT - div(k*grad(T)) = q    =>    rho*cp*ddt - (k/r)*ddr(r*ddr(T)) = q

Here is an example of a mesh with 5 nodes, represented by (o) and two "ghost"
nodes before and after the boundary nodes 0 and 4.

    radius=   rIn  ...                           ...  rOut

    index =   0         1         2         3         4 (= N-1)

    x_ _ _ _ _o_________o_________o_________o_________o_ _ _ _ _x



Let us discretise the eq. term by term with a finite difference scheme. For a 
radial node i, one has:

    1) rho*cp*ddtT -> rho_i*cp_i*( T_i - T_i_old )/dt =

        =       T_i *       ( rho_i*cp_i/dt )
            -   T_(i, old)* (rho_i*cp_i/dt)


    
    2) (k/r)*ddr(r*ddr(T)) = k*( (1/r)*(ddr(T)) + d2dr2(T) ) =>
        
        

        2.1) (k/r)*ddr(T) -> k*( T_(i+1) - T_(i-1) )/(2*dr*r_i) [CDS scheme] =

            =       T_(i+1) *   (  k/(2*dr*r_i) )
                +   T_(i-1) *   ( -k/(2*dr*r_i) )
        
        

        2.2) k*d2dr2(T) = ddr(ddr(T)) -> k*( ddr(T)_(i+1) - ddr(T)_i )/dr
            
            = k_i*( (T_(i+1) - T_(i))/dr - (T_(i) - T_(i-1))/dr )/dr =

            =       T_(i+1) *   (   k/(dr^2) )
                +   T_i *       ( -2k/(dr^2) )
                +   T_(i-1) *   (   k/(dr^2) )



        So (k/r)*ddr(r*ddr(T)) ->

            ->      T_(i+1) *   (   k*( 1/(dr^2)   + 1.0/(2*dr*r_i) ) )
                +   T_i *       ( -2k*( 1/(dr^2) )                    )
                +   T_(i-1) *   (   k*( 1/(dr^2)   - 1.0/(2*dr*r_i) ) )



    Recall that (k/r)*ddr(r*ddr(T)) appears with the - sign in the heat
    equation, so all the coefficients from 2) need to be changed in sign. 
    Combining the coefficients from 1) and 2) and 2.2) we get equation [I]:

        rho*cp*ddtT - div(k*grad(T)) = q ->

        ->      T_(i+1) *   ( -k*( 1/(dr^2)     + 1.0/(2*dr*r_i) ) )
            +   T_i *       ( 2k*( 1/(dr^2) )   + rho_i*cp_i/dt    )
            +   T_(i-1) *   ( -k*( 1/(dr^2)     - 1.0/(2*dr*r_i) ) )
            =
                q_i
            +   T_(i, old)*(rho_i*cp_i/dt)

    If we set:

        A = k/(dr^2)
        B = 1.0/(2*dr*r_i)
        X = rho_i*cp_i/dr

    Then we can re-write the whole thing as:                                 

        rho*cp*ddtT - div(k*grad(T)) = q ->

        ->      T_(i+1) *   ( -A + B )
            +   T_i *       ( 2A + X )
            +   T_(i-1) *   ( -A - B )
            =
                q_i
            +   T_(i, old)*X                                                [I]



Equation I is valid everywhere, even at the boundaries. However, at the
boundaries, T_N and T_(-1) do not exists. This is where ghost nodes come into
play. Let us have the following boundary conditions:

    1) Zero gradient at i = 0 : (i.e. Neumann)

        ddr(T)_(r_0) = 0 -> (T_1 - T(_(-1))/2*dr = 0 - > T_(-1) = T_1

    T(-1) does not exists and is our left ghost node. It does not matter,
    we obtained a relationship between T_1 and T_(-1)! Let us substitute this
    relationship in equation [I] and we get the equation at the left boundary:

        rho*cp*ddtT - div(k*grad(T)) = q -> (zeroGradient at r = rIn) ->

        ->      T_1 *       ( -2A )
            +   T_0 *       ( 2A + X )
            =
                q_i  
            +   T_(0, old)*X                                               [II]



    2) Convective boundary condition at i = N-1 : (i.e. Robin)

        k*ddr(T)_N = hExt*(TExt - T_(N-1))

    with hExt being the heat transfer coefficients of the medium after
    rOut and TExt being the temperature of the medium. As before, we can
    discretize it as:

        k*( T_N - T_(N-2) )/(2*dr) = hExt*(TExt - T_(N-1)) ->

        -> T_N = T_(N-2) + (2*dr/k)*hExt*(TExt -T_(N-1))

    If we substitute the expression for the ghost node T_N in equation I, we 
    obtain the BC. For convenience, let us rename C = (2*dr/k):

        rho*cp*ddtT - div(k*grad(T)) = q -> (convective at r = rOut) ->

                T_(N-1) * ( 2A + hExt*C*(A+B) + X )
            +   T_(N-2) * ( -2A )
            =
                q_(N-1)
            +   T_(N-1, old)*X
            +   TExt*hExt*C*(A+B)                                         [III]



So, there you have it, the boundary coefficients (II, III) and the bulk
coefficients (I).
Ok, now what if you want to do a cylinder? Well, you don't need to change
anything as there are no terms that contain 1/r_i for the zeroGradient
expression at r = r_1. Ok again, what if you want to do a proper nuclear fuel 
pin? Well, you should be able to understand what I did in the code. I 
used a convective boundary condition at the fuel outer surface where TExt is 
the  actual T_(i+1), while the fuel T_(i+1) is treated as a ghost node. Same 
for the inner side of the cladding, yet the coefficients are swapped for 
obvious reasons. For the outer cladding surface, again, it is a convective 
boundary conditions, yet it is complicated by the fact that I might have a mix 
of vapour and liquid, at different temperatures, both contacting the pin. For 
simplicity, I assume a single cladding outer temperature. Then, due to the
additivity of heat transfer phenomena, I simply consider:

    k*( T_N - T_(N-2) )/2*dr = 
        =   hVap*fracVap*(TVap - T_(N-1)) + hLiq*fracLiq*(TLiq - T_(N-1))

with fracVap and fracLiq being the fractions of structure interfacial area
that is in contact with either the vapour or the liquid. More in general, the
treatment of the convective boundary conditions can be generalized for any
number of fluids contacting the pin.
Let us have:

    k*( T_N - T_(N-2) )/2*dr = 
        =   SUM_j [hf_j(T_(ext, j) - T_(N-1))]

Where SUM_j denotes a summation over the j indices. hf_j is the heat
transfer coefficient betwee the j-th fluid and the structure, multiplied by
the fraction of interfacial area of the structure that is in contact with the
j-th fluid. Then, this can be re-written as:

    k*( T_N - T_(N-2) )/2*dr = 
        =   SUM_j[hf_j*T_(ext, j)] - SUM_j[hf_j]*T_(N-1)
        =   HTSum - HSum*T_(N-1)

which is the notation that is used in the code implementation. This grants more
generality to this class, which does not need to now any details on how many 
fluid are there, nor what are the interfacial area fractions and so on, as
these need to be passed by the user at a higher level (i.e., when calling the
correct(HTSum, HSum) function in the main program).

\*---------------------------------------------------------------------------*/