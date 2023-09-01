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

#include "structure.H"
#include "fvmDdt.H"
#include "fvcDiv.H"
#include "fvmSup.H"
#include "fluid.H"
#include "myOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(structure, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::structure::structure
(
    const dictionary& dict,
    const fvMesh& mesh,
    volScalarField& powerDensityNeutronics
)
:
    phaseBase
    ( 
        dict,
        mesh,
        "structure",
        zeroGradientFvPatchScalarField::typeName
    ),
    regions_(0),
    cells_(0),
    powerDensityNeutronics_(powerDensityNeutronics),
    Dh_
    (
        IOobject
        (
            "Dh",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("", dimLength, SMALL),
        zeroGradientFvPatchScalarField::typeName
    ),
    lDh_
    (
        IOobject
        (
            "lDh",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedVector("", dimLength, vector(SMALL, SMALL, SMALL)),
        zeroGradientFvPatchScalarField::typeName
    ),
    heatFlux_
    (
        IOobject
        (
            "heatFlux.structure",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimPower/dimArea, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    Tact_
    (
        IOobject
        (
            "T.activeStructure",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ, //- Constructed from powerModel initial Ts
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimTemperature, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    iAact_
    (
        IOobject
        (
            "volumetricArea.activeStructure",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("", dimArea/dimVol, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    alphapas_
    (
        IOobject
        (
            "alpha.passiveStructure",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        *this
    ),
    Tpas_
    (
        IOobject
        (
            "T.passiveStructure",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("T", dimTemperature, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    TFuelAv_
    (
        IOobject
        (
            "T.fuelAvForNeutronics",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("T", dimTemperature, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    TCladAv_
    (
        IOobject
        (
            "T.cladAvForNeutronics",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("T", dimTemperature, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    iApas_
    (
        IOobject
        (
            "volumetricArea.passiveStructure",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("", dimArea/dimVol, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    alphaRhoCppas_
    (
        IOobject
        (
            "alphaRhoCp.passiveStructure",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar
        ("", dimEnergy/dimVol/dimTemperature, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    Twall_
    (
        IOobject
        (
            "T.wall",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("T", dimTemperature, 0),
        zeroGradientFvPatchScalarField::typeName
    )
{
    bool foundAtLeastOnePassivePropertiesDict(false);

    //- Correct header class names. This is necessary for all volScalarField
    //  of the structure that are read/written from disk, as, for whatever
    //  reason I can't figure out, their class header defaults to "byZone"
    //  instead of "volScalarField". This can mess up field reading.
    //  Now the fun part. Changing the headerClassName via the headerClassName
    //  method only works for member fields BUT NOT for the base volScalarField
    //  (i.e. for the alpha.structure field). I tried all of the three variants
    //  below but literally no variant works... What the actual fuck in the
    //  name of sweet baby Jesus? Jasak help pls
    this->headerClassName() = "volScalarField";
    this->IOobject::headerClassName() = "volScalarField";
    this->volScalarField::headerClassName() = "volScalarField";

    //- At least it worked for this field
    Tpas_.headerClassName() = "volScalarField";

    bool alphaHeaderOk(this->typeHeaderOk<volScalarField>(true));
    bool TpasHeaderOk(Tpas_.typeHeaderOk<volScalarField>(true));

    forAll(dict.toc(), i)
    {
        //- The sturctureProperties dictionary keys consist in cellZone names.
        //  These keys can consist in either the name of a single cellZone, or
        //  for brevitiy, multplie cellZone names in a single string (e.g. if 
        //  multiple cellZones share the same structure properties, e.g. 
        //  powerModels, void fraction, hydraulic diameter, etc.). The format
        //  for the latter is "zone0:zone1:zone2:...:zoneN", i.e. the colon is
        //  the separation character between zone names
        word key(dict.toc()[i]);
        if (key == "type") continue;
        if (key == "powerOffCriterionModel") continue;
        if (key == "heatExchangers") continue;
        const dictionary& zoneDict(dict.subDict(key));
        
        wordList zones(myOps::split<word>(key, ':'));

        forAll(zones, i)
        {
            word zone(zones[i]);
            const labelList& zoneCellList(mesh.cellZones()[zone]);
        
            //- Construct cellLists_, cells_, cellFields_
            regions_.append(zone);
            cellLists_.insert
            (
                zone,
                zoneCellList
            );
            forAll(zoneCellList, i)
            {
                cells_.append(zoneCellList[i]);
            }
            scalarField zoneCellField(mesh.cells().size(), 0.0);
            forAll(zoneCellList, j)
            {
                zoneCellField[zoneCellList[j]] = 1.0;
            }
            cellFields_.insert
            (
                zone,
                zoneCellField
            );
            
            //- Set volumeFraction of the structure, hydraulic diameter.
            //  volumeFraction set here only if alpha.structure not found
            //  on disk
            scalar Dh(zoneDict.get<scalar>("Dh"));
            if (!alphaHeaderOk)
            {
                scalar alpha(zoneDict.get<scalar>("volumeFraction"));
                forAll(zoneCellList, j)
                {   
                    label cellj(zoneCellList[j]);
                    (*this)[cellj] = alpha;
                    Dh_[cellj] = Dh;
                }
            }
            else
            {
                forAll(zoneCellList, j)
                {   
                    label cellj(zoneCellList[j]);
                    Dh_[cellj] = Dh;
                }
            }
            

            //- Set HashTable of volumeFraction volScalarField indexed by 
            //  region (i.e. zone) name
            volScalarField zoneAlphaField(*this);
            zoneAlphaField.primitiveFieldRef() *= zoneCellField;
            zoneAlphaField.correctBoundaryConditions();
            alphaFields_.insert
            (
                zone,
                zoneAlphaField
            );

            //- Read momentum source
            if (zoneDict.found("momentumSource"))
            {
                pumps_.insert
                (
                    zone,
                    pump
                    (
                        mesh_,
                        zoneDict,
                        zoneCellList
                    )
                );

                if (!momentumSourcePtr_.valid())
                {
                    momentumSourcePtr_.reset
                    (
                        new volVectorField
                        (
                            IOobject
                            (
                                "momentumSource",
                                mesh_.time().timeName(),
                                mesh_
                            ),
                            mesh_,
                            dimensionedVector
                            (
                                "momentumSource", 
                                dimDensity*dimVelocity/dimTime, 
                                vector::zero
                            ),
                            zeroGradientFvPatchVectorField::typeName
                        )
                    );
                }
            }

            //- Set passive properties fields if keywords present
            if (zoneDict.isDict("passiveProperties"))
            {
                foundAtLeastOnePassivePropertiesDict = true;

                const dictionary& pasDict
                (
                    zoneDict.subDict("passiveProperties")
                );

                //- Set interfacial area and rhoCp (later, alphaRhoCp) of the
                //  passive subStructure
                scalar iApas(pasDict.get<scalar>("volumetricArea"));
                scalar rhoCppas(0);
                if 
                (
                    pasDict.found("rho") 
                and pasDict.found("Cp") 
                and pasDict.found("rhoCp")
                )
                {
                    FatalErrorInFunction
                        << "Structure region: " << zone << " -> "
                        << "specify either rhoCp or both rho and Cp but not "
                        << "all of them"
                        << exit(FatalError);
                }
                else
                {
                    if (pasDict.found("rho") and pasDict.found("Cp"))
                    {
                        rhoCppas =
                            pasDict.get<scalar>("rho")*
                            pasDict.get<scalar>("Cp");
                    }
                    else if (pasDict.found("rhoCp"))
                    {
                        rhoCppas = pasDict.get<scalar>("rhoCp");
                    }
                    else
                    {
                        FatalErrorInFunction
                            << "Structure region: " << zone << " -> "
                            << "specify either rhoCp or both rho and Cp"
                            << exit(FatalError);
                    }
                }
                
                forAll(zoneCellList, j)
                {
                    label cellj(zoneCellList[j]);
                    iApas_[cellj] = iApas;
                    alphaRhoCppas_[cellj] = rhoCppas;   //- The multiplication
                                                        //  by alphapas_ will
                                                        //  be done at the end
                }

                //- Only adjust passive subStructure volume fraction if keyword
                //  found in the passive properties dictionary
                if (pasDict.found("volumeFraction"))
                {
                    scalar alphapas(pasDict.get<scalar>("volumeFraction"));
                    forAll(zoneCellList, j)
                    {
                        label cellj(zoneCellList[j]);
                        alphapas_[cellj] = alphapas;
                    }
                }
                //- If the alpha.structure file does not exist in the initial 
                //  time step folder, and a volumeFraction keyword is not found
                //  in the passive subStructure dict, init value to region 
                //  alpha value, read before

                //- NOTE: Oddly enough, the class keyword in the file headers
                //  of alphapas_, Tpas_, are set to byZone rather than 
                //  volScalarField. This is due to some weird dark magic of
                //  the runTimeSelection mechanism, which I have no will to
                //  investigate. Thus, to check that fields are present, the
                //  headerType that needs to be looked for is byZone, not
                //  volScalarField. Thanks to Carlo for finding out about this
                //  odd behaviour
                else if (!alphaHeaderOk)
                {
                    scalar alpha(zoneDict.get<scalar>("volumeFraction"));
                    forAll(zoneCellList, j)
                    {
                        label cellj(zoneCellList[j]);
                        alphapas_[cellj] = alpha;
                    }
                }

                //- If the passive subStructure temperature field does not 
                //  exist in the initial time step folder, get it from dict
                if (!TpasHeaderOk)
                {
                    scalar Tpas(pasDict.get<scalar>("T"));
                    forAll(zoneCellList, j)
                    {
                        label cellj(zoneCellList[j]);
                        Tpas_[cellj] = Tpas;
                    }
                }
            }

            //- Now for the rotation matrices to move from the global to the
            //  local reference frame

            //- Lambda function to rotate a vector around an axis by a certain
            //  angle in radians, counter-clockwise
            auto rotateCCWAroundAxisByAngle = []
            (
                vector& v,
                const vector& axis,
                scalar angle
            )
            {
                scalar c(Foam::cos(angle));
                scalar s(Foam::sin(angle));
                scalar x(axis[0]);
                scalar y(axis[1]);
                scalar z(axis[2]);
                tensor R
                (
                    c+sqr(x)*(1-c),     x*y*(1-c)-z*s,      x*z*(1-c)+y*s,

                    y*x*(1-c)+z*s,      c+sqr(y)*(1-c),     y*z*(1-c)-x*s,

                    z*x*(1-c)-y*s,      z*y*(1-c)+x*s,      c+sqr(z)*(1-c)

                );
                v = R & v;
            };

            //- Construct rotation matrices
            bool foundLocalX(zoneDict.found("localX"));
            bool foundLocalZ(zoneDict.found("localZ"));
            if (foundLocalX or foundLocalZ)
            {
                if (!Rl2gPtr_.valid())
                    Rl2gPtr_.reset
                    (
                        new volTensorField
                        (
                            IOobject
                            (
                                "R.localToGlobal",
                                mesh.time().timeName(),
                                mesh
                            ),
                            mesh,
                            dimensionedTensor
                            (
                                "", 
                                dimless, 
                                tensor
                                (
                                    1,0,0,
                                    0,1,0,
                                    0,0,1
                                )
                            ),
                            zeroGradientFvPatchTensorField::typeName
                        )
                    );
                if (!Rg2lPtr_.valid())
                    Rg2lPtr_.reset
                    (
                        new volTensorField
                        (
                            IOobject
                            (
                                "R.globalToLocal",
                                mesh.time().timeName(),
                                mesh
                            ),
                            mesh,
                            dimensionedTensor
                            (
                                "", 
                                dimless, 
                                tensor
                                (
                                    1,0,0,
                                    0,1,0,
                                    0,0,1
                                )
                            ),
                            zeroGradientFvPatchTensorField::typeName
                        )
                    );
                volTensorField& Rl2g(Rl2gPtr_());
                volTensorField& Rg2l(Rg2lPtr_());

                vector localX
                (
                    zoneDict.lookupOrDefault<vector>("localX", vector(1,0,0))
                );
                localX /= mag(localX);
                vector localZ
                (
                    zoneDict.lookupOrDefault<vector>("localZ", vector(0,0,1))
                );
                localZ /= mag(localZ);
                //- Non orthogonality absorption
                if (foundLocalX and !foundLocalZ)
                {
                    //- Absorb non-orthogonalities in Z
                    localZ -= (localX&localZ)*localX/mag(localX);
                    localZ /= mag(localZ);
                }
                else if (!foundLocalX and foundLocalZ)
                {
                    //- Absorb non-orthogonalities in X
                    localX -= (localX&localZ)*localZ/mag(localZ);
                    localX /= mag(localX);
                }
                else
                {
                    //- Split non-orthogonalities equally among X and Z by 
                    //  rotating them in the plane they lie in by an angle
                    //  computed so that, after the rotation, they will be
                    //  orthogonal
                    scalar deltaTheta
                    (
                        (
                            constant::mathematical::pi/2.0
                        -   Foam::acos(localX&localZ)
                        )/2.0
                    );
                    vector axis(localX ^ localZ);
                    axis /= mag(axis);
                    rotateCCWAroundAxisByAngle(localX, axis, -deltaTheta);
                    rotateCCWAroundAxisByAngle(localZ, axis, deltaTheta);
                }
                
                //- Compute third axis
                vector localY(localZ ^ localX);
                localY /= mag(localY);

                //- Construct transformation matrices

                //- The basis change matrix is the transformation to move from
                //  the local reference frame to the global one. It is constructed
                //  by simply arranging the local basis vectors (expressed in 
                //  global reference frame coordinates) in columns. Since these are 
                //  orthonormal, the matrix is orthonormal and its inverse is equal
                //  to its transpose. Thus, the transformation matrix to move from
                //  the global to the local frame is the transpose of the one to
                //  move from the local to the global frame
                tensor Rl2gi
                (
                    localX[0], localY[0], localZ[0],
                    localX[1], localY[1], localZ[1],
                    localX[2], localY[2], localZ[2]
                );
                tensor Rg2li = Rl2gi.T();

                forAll(zoneCellList, j)
                {
                    label cellj(zoneCellList[j]);
                    Rl2g[cellj] = Rl2gi;
                    Rg2l[cellj] = Rg2li;
                }
                Rl2g.correctBoundaryConditions();
                Rg2l.correctBoundaryConditions();
            }

            //- Construct lDh_ (for isotropic structures each component of
            //  lDh_ is equal to Dh cell by cell)
            vector lDhAnisotropy
            (
                zoneDict.lookupOrDefault<vector>
                (
                    "localDhAnisotropy", 
                    vector::one
                )
            );
            forAll(zoneCellList, j)
            {
                label cellj(zoneCellList[j]);
                lDh_[cellj][0] = lDhAnisotropy[0]*Dh_[cellj];
                lDh_[cellj][1] = lDhAnisotropy[1]*Dh_[cellj];
                lDh_[cellj][2] = lDhAnisotropy[2]*Dh_[cellj];
            }
            
            //- Construct global tortuosity tensor by transforming it from the
            //  local frame (as provided in the dictionary) to the global one.
            //  Recall that if R is the transformation matrix to rotate a 
            //  vector from the local to the global frame, a local tensor Q can
            //  be rotated to the global frame via R & Q & R.T(). In this case,
            //  R = Rl2g. Note that while Rl2g.T() = Rg2l, Rl2g.T() was kept 
            //  for clarity
            if (zoneDict.found("localTortuosity"))
            {
                if (!tortuosityPtr_.valid())
                    tortuosityPtr_.reset
                    (
                        new volTensorField
                        (
                            IOobject
                            (
                                "tortuosity",
                                mesh.time().timeName(),
                                mesh
                            ),
                            mesh,
                            dimensionedTensor
                            (
                                "", 
                                dimless, 
                                tensor(1,0,0,0,1,0,0,0,1)
                            ),
                            zeroGradientFvPatchScalarField::typeName
                        )
                    );
                volTensorField& tortuosity = tortuosityPtr_();
                vector lTortuosityVector
                (
                    zoneDict.get<vector>("localTortuosity")
                );
                tensor lTortuosity
                (
                    lTortuosityVector[0], 0, 0,
                    0, lTortuosityVector[1], 0,
                    0, 0, lTortuosityVector[2]
                );
                if (this->hasLocalReferenceFrame())
                {
                    forAll(zoneCellList, j)
                    {
                        label cellj(zoneCellList[j]);
                        tortuosity[cellj] = 
                            Rl2g()[cellj] & lTortuosity & Rg2l()[cellj];
                    }
                }
                else
                {
                    forAll(zoneCellList, j)
                    {
                        label cellj(zoneCellList[j]);
                        tortuosity[cellj] = lTortuosity;
                    }
                }
            }
        }
    }

    //- Don't write passiveStructure temperature field to files if no passive
    //  properties were specified in the phasePropertiesDict
    if (!foundAtLeastOnePassivePropertiesDict)
    {
        Tpas_.writeOpt() = IOobject::NO_WRITE;
    }

    //- Correct BCs of fields set cell-by-cell
    this->correctBoundaryConditions();
    Dh_.correctBoundaryConditions();
    Tpas_.correctBoundaryConditions();
    iApas_.correctBoundaryConditions();
    alphapas_.correctBoundaryConditions();
    lDh_.correctBoundaryConditions();
    if (tortuosityPtr_.valid())
    {
        tortuosityPtr_().correctBoundaryConditions();
    }
    if (momentumSourcePtr_.valid())
    {
        momentumSourcePtr_().correctBoundaryConditions();
    }

    //- Multiply rhoCppas by alphapas to get actual volumetric heat capacity
    //  of the passive subStructure, limit to avoid 0 matrix coefficients when
    //  solving for passive subSubstructure energy equation
    alphaRhoCppas_ = 
        Foam::max
        (
            alphapas_*alphaRhoCppas_, 
            dimensionedScalar("", dimEnergy/dimVolume/dimTemperature, 1e-69)
        );
    alphaRhoCppas_.correctBoundaryConditions();

    //- Construct powerModels
    //  First, I need a list of all the typeNames in the various subDicts
    wordList powerModelTypes(0);
    forAll(dict.toc(), i)
    {
        word key(dict.toc()[i]);
        if (key == "type") continue;
        const dictionary regionDict(dict_.subDict(key));
        if (regionDict.isDict("powerModel"))
        {
            const dictionary& powerModelDict
            (
                regionDict.subDict("powerModel")
            );

            word powerModelType(powerModelDict.get<word>("type"));

            bool found(false);
            forAll(powerModelTypes, i)
            {
                if (powerModelTypes[i] == powerModelType)
                {
                    found = true;
                    break;
                }
            }
            if (!found)
            {
                powerModelTypes.append(powerModelType);
            }
        }
    }

    forAll(powerModelTypes, i)
    {
        //- Then, for each typeName, create a dict of subDicts
        dictionary dicts;
        word powerModelType(powerModelTypes[i]);

        forAll(dict.toc(), i)
        {
            word key(dict.toc()[i]);
            if (key == "type") continue;
            const dictionary& regionDict(dict_.subDict(key));
            if (regionDict.isDict("powerModel"))
            {
                const dictionary& powerModelDict
                (
                    regionDict.subDict("powerModel")
                );

                wordList zones(myOps::split<word>(key, ':'));

                if (powerModelDict.get<word>("type") == powerModelType)
                {
                    forAll(zones, i)
                    {
                        dicts.add(zones[i], powerModelDict);
                    }
                }
            }
        }

        //- Finally, construct one powerModel per type
        powerModels_.insert
        (
            powerModelType,
            powerModel::New
            (
                *this,
                dicts
            )
        );
    }
 
    //- Adjust iAact that was left out as it is a member of powerModel
    //  and set intitial Tact_
    forAllIter
    (
        powerModelTable,
        powerModels_,
        iter
    )
    {
        iAact_ += iter()->iA();
        iter()->correctT(Tact_);
    }
    iAact_.correctBoundaryConditions();
    Tact_.correctBoundaryConditions();

    this->constructHeatExchangers();

    if (this->dict().isDict("powerOffCriterionModel"))
    {
        powerOffCriterionModelPtr_.reset
        (
            powerOffCriterionModel::New
            (
                mesh,
                this->dict().subDict("powerOffCriterionModel")
            )
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::structure::localToGlobalRotateField(volTensorField& field) const
{
    if(this->hasLocalReferenceFrame())
    {
        field = Rl2g() & field & Rg2l();
        field.correctBoundaryConditions();
    }
}

void Foam::structure::localToGlobalRotateField(volVectorField& field) const
{
    if(this->hasLocalReferenceFrame())
    {
        field = Rl2g() & field;
        field.correctBoundaryConditions();
    }
}

void Foam::structure::globalToLocalRotateField(volTensorField& field) const
{
    if(this->hasLocalReferenceFrame())
    {
        field = Rg2l() & field & Rl2g();
        field.correctBoundaryConditions();
    }
}

void Foam::structure::globalToLocalRotateField(volVectorField& field) const
{
    if(this->hasLocalReferenceFrame())
    {
        field = Rg2l() & field;
        field.correctBoundaryConditions();
    }
}

const Foam::volVectorField& Foam::structure::momentumSource()
{
    momentumSourcePtr_() *= 0.0;
    forAllIter
    (
        pumpTable,
        pumps_,
        iter
    )
    {
        iter().correct(momentumSourcePtr_());
    }
    return momentumSourcePtr_();
}

void Foam::structure::constructHeatExchangers()
{
    if (this->dict().found("heatExchangers"))
    {
        const dictionary& HXDicts(this->dict().subDict("heatExchangers"));
        forAll(HXDicts.toc(), j)
        {
            word HXKey(HXDicts.toc()[j]);
            const dictionary& HXDict(HXDicts.subDict(HXKey));
            heatExchangers_.insert
            (
                HXKey,
                heatExchanger
                (
                    mesh_,
                    HXDict
                )
            );
        }

        //- If HXs were specified, created the HX-specific fields, THXPtr_ and
        //  iAHXPtr_
        THXPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "T.HX",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar("", dimTemperature, 0),
                zeroGradientFvPatchScalarField::typeName
            )
        );
        iAHXPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar("", dimArea/dimVol, 0),
                zeroGradientFvPatchScalarField::typeName
            )
        );
        volScalarField& iAHX(iAHXPtr_());
        forAllIter
        (
            heatExchangerTable,
            heatExchangers_,
            iter
        )
        {
            iAHX += iter().iA();
        }
        iAHX.correctBoundaryConditions();
    }
}

void Foam::structure::correct
(
    const volScalarField& HT,
    const volScalarField& H
)
{
    //- Correct powerModels & correct active structure surface temperature
    forAllIter
    (
        powerModelTable,
        powerModels_,
        iter
    )
    {
        iter()->correct(HT, H);
        iter()->correctT(Tact_);
    }
    Tact_.correctBoundaryConditions();

    //- Correct HX surface temperature
    if (THXPtr_.valid())
    {
        forAllIter
        (
            heatExchangerTable,
            heatExchangers_,
            iter
        )
        {
            //- Unlike powerModels, the functionalities of correct and correctT
            //  are merged into a single correct. Why didn't I do so in 
            //  powerModels? A very good question for past Stefan which present
            //  Stefan does not know the answer to. Probably for added
            //  """flexibility""" given that powerModel is runTimeSelectable 
            //  and heatExchanger is not?
            iter().correct(HT, H, THXPtr_());
        }
    }
    
    if (Tpas_.writeOpt() == IOobject::AUTO_WRITE)
    {
        //- Correct inert subStructure. What follows is the equivalent of doing
        //  the following:
        /*
            fvScalarMatrix pasEqn
            (
                    fvm::ddt(alphaRhoCppas_, Tpas_)
                ==
                    iApas_*HT
                -   fvm::Sp(iApas_*H, Tpas_)
            );
            pasEqn.solve();
        */
        //  Except, it is faster like this rather than to solve an equation 
        //  over the entire mesh, as the passive subStructure might not exist
        //  everywhere

        scalar dt(mesh_.time().deltaT().value());
        const volScalarField& Tpas0(Tpas_.oldTime()); 
        forAll(cells_, i)
        {
            label celli(cells_[i]);
            const scalar& iA(iApas_[celli]);
            if (iA == 0) continue;  //- Avoid solving where the passive 
                                    //  structure does not exist
            scalar alphaRhoCpByDt(alphaRhoCppas_[celli]/dt);
            Tpas_[celli] = 
                (
                    iA*HT[celli] 
                +   alphaRhoCpByDt*Tpas0[celli]
                )/
                (alphaRhoCpByDt + iA*H[celli]);
        }
        Tpas_.correctBoundaryConditions();
    }

    //- Update Twall, heatFlux    
    forAll(cells_, i)
    {
        const label& celli(cells_[i]);
        scalar& Twall(Twall_[celli]);

        //- Set wall temperature as  power structure surface temperature
        Twall = Tact_[celli];

        //- Update heat flux (mostly for extra info purposes, maybe only
        //  used by the Shah pool boiling model under some circumstances).
        heatFlux_[celli] = H[celli]*Twall-HT[celli];
    }
    Twall_.correctBoundaryConditions();
    heatFlux_.correctBoundaryConditions();
}

Foam::tmp<Foam::volScalarField> Foam::structure::explicitHeatSource
(
    const fluid& fluid,
    const volScalarField& H
) const
{
    const volScalarField& T(fluid.thermo().T());
    tmp<volScalarField> tQ
    (
        new volScalarField
        (
                iAact_*H*(Tact_-T)
            +   iApas_*H*(Tpas_-T)
        )
    );

    if (THXPtr_.valid())
    {
        //- Read comment on the same piece of code in
        //  linearizedSemiImplicitHeatSource
        if (Foam::max(THXPtr_()).value() > 1e-69)
        {
            volScalarField& Q = tQ.ref();
            const volScalarField& iAHX(iAHXPtr_());
            const volScalarField& THX(THXPtr_());
            Q += iAHX*H*(THX-T);
        }
    }
    return tQ;
}

Foam::tmp<Foam::fvScalarMatrix> 
Foam::structure::linearizedSemiImplicitHeatSource
(
    fluid& fluid,
    const volScalarField& H
) const
{
    volScalarField& he(fluid.thermo().he());
    const volScalarField& T(fluid.thermo().T());
    const volScalarField& Cp(fluid.Cp());
    tmp<fvScalarMatrix> tQ
    (
        new fvScalarMatrix
        (
                //- Source/sink due to active subStructure
                iAact_*H*(Tact_-T+he/Cp)
            -   fvm::Sp(iAact_*H/Cp, he)
                //- Source/sink due to passive subStructure
            +   iApas_*H*(Tpas_-T+he/Cp)
            -   fvm::Sp(iApas_*H/Cp, he)
        )
    );
    
    //- Add HX contribution
    if (THXPtr_.valid())
    {
        //- What is the deal with this funky max? Well, the heatExchanger class
        //  does not set THX at construction (as the class itself only accesses
        //  THX though the correct function), so that, if the 
        //  linearizedSemiImplicitHeatSource function is called BEFORE
        //  structure correct (depending on how the EEqn of the solver is
        //  implemented), then at the VERY FIRST PIMPLE iteration of the VERY
        //  FIRST time-step, THX = 0 and the resulting massive heat transfer 
        //  with the fluid will make the simulation explode. Thus, if THX is 0
        //  everywhere (i.e. if it was not set yet at all), do not add HX
        //  contributions. I could do this by checking if I am in the first 
        //  PIMPLE iteration of the first simulation time step, but whatever.
        //  Ok, ok, what if I wanted to set it at construction? Well, you 
        //  cannot because you inherently need T and H*T to correct the HX
        //  temperatures, and those need to be passed somehow. Wanna use 
        //  registry lookup? Very bad idea, all the heat transfer tables
        //  that contain the Hs are constructed only after the structure (and
        //  it cannot be otherwise). So, just shut up and go with the flow!
        if (Foam::max(THXPtr_()).value() > 1e-69)
        {
            fvScalarMatrix& Q = tQ.ref();
            const volScalarField& iAHX(iAHXPtr_());
            const volScalarField& THX(THXPtr_());
            Q += 
                    iAHX*H*(THX-T+he/Cp)
                -   fvm::Sp(iAHX*H/Cp, he);
        }
    }

    return tQ;
}

void Foam::structure::checkPowerOff()
{
    if (powerOffCriterionModelPtr_.valid())
    {
        if (powerOffCriterionModelPtr_->powerOffCriterion())
        {
            forAllIter
            (
                powerModelTable,
                powerModels_,
                iter
            )
            {
                iter()->powerOff();
            }
        }
    }
}

// ************************************************************************* //
