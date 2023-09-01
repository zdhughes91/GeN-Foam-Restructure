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

#include "byZoneStructure.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcDiv.H"
#include "myOps.H"

//- From forward declarations in structureModel.C
#include "fluid.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace structureModels
{
    defineTypeNameAndDebug(byZone, 0);
    addToRunTimeSelectionTable
    (
        structureModel,
        byZone,
        structureModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::structureModels::byZone::byZone
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    structureModel
    ( 
        dict,
        mesh
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

            bool foundLocalX(zoneDict.found("localX"));
            bool foundLocalZ(zoneDict.found("localZ"));
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
            tensor Rl2g
            (
                localX[0], localY[0], localZ[0],
                localX[1], localY[1], localZ[1],
                localX[2], localY[2], localZ[2]
            );
            tensor Rg2l = Rl2g.T();

            forAll(zoneCellList, j)
            {
                label cellj(zoneCellList[j]);
                Rl2g_[cellj] = Rl2g;
                Rg2l_[cellj] = Rg2l;
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
            vector lTortuosityVector
            (
                zoneDict.lookupOrDefault<vector>
                (
                    "localTortuosity", 
                    vector::one
                )
            );
            tensor lTortuosity(tensor::zero);
            lTortuosity[0] = lTortuosityVector[0];
            lTortuosity[4] = lTortuosityVector[1];
            lTortuosity[8] = lTortuosityVector[2];
            tensor tortuosity = Rl2g & lTortuosity & Rg2l;

            forAll(zoneCellList, j)
            {
                label cellj(zoneCellList[j]);
                tortuosity_[cellj] = tortuosity;
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
    Rg2l_.correctBoundaryConditions();
    lDh_.correctBoundaryConditions();
    tortuosity_.correctBoundaryConditions();
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
    forAllIter
    (
        powerModelTable,
        powerModels_,
        iter
    )
    {
        iAact_ += iter()->iA();
    }
    iAact_.correctBoundaryConditions();

    Rl2g_.correctBoundaryConditions();
    Rg2l_.correctBoundaryConditions();

    this->constructHeatExchangers();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// ************************************************************************* //
