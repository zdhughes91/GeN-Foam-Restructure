/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
|                                                                             |
|    Built on OpenFOAM v2112                                                  |
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

#include "argList.H"
#include "Time.H"
#include "syncTools.H"
#include "faceSet.H"
#include "pointSet.H"
#include "meshTools.H"
#include "polyTopoChange.H"
#include "polyRemoveFace.H"
#include "polyModifyFace.H"
#include "indirectPrimitivePatch.H"
#include "processorPolyPatch.H"
#include "localPointRegion.H"
#include "duplicatePoints.H"
#include "ReadFields.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "processorMeshes.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void insertDuplicateMerge
(
    const polyMesh& mesh,
    const labelList& boundaryFaces,
    const labelList& duplicates,
    polyTopoChange& meshMod
)
{
    const faceList& faces = mesh.faces();
    const labelList& faceOwner = mesh.faceOwner();
    const faceZoneMesh& faceZones = mesh.faceZones();

    forAll(duplicates, bFacei)
    {
        label otherFacei = duplicates[bFacei];

        if (otherFacei != -1 && otherFacei > bFacei)
        {
            // Two duplicate faces. Merge.

            label face0 = boundaryFaces[bFacei];
            label face1 = boundaryFaces[otherFacei];

            label own0 = faceOwner[face0];
            label own1 = faceOwner[face1];

            if (own0 < own1)
            {
                // Use face0 as the new internal face.
                label zoneID = faceZones.whichZone(face0);
                bool zoneFlip = false;

                if (zoneID >= 0)
                {
                    const faceZone& fZone = faceZones[zoneID];
                    zoneFlip = fZone.flipMap()[fZone.whichFace(face0)];
                }

                meshMod.setAction(polyRemoveFace(face1));
                meshMod.setAction
                (
                    polyModifyFace
                    (
                        faces[face0],           // modified face
                        face0,                  // label of face being modified
                        own0,                   // owner
                        own1,                   // neighbour
                        false,                  // face flip
                        -1,                     // patch for face
                        false,                  // remove from zone
                        zoneID,                 // zone for face
                        zoneFlip                // face flip in zone
                    )
                );
            }
            else
            {
                // Use face1 as the new internal face.
                label zoneID = faceZones.whichZone(face1);
                bool zoneFlip = false;

                if (zoneID >= 0)
                {
                    const faceZone& fZone = faceZones[zoneID];
                    zoneFlip = fZone.flipMap()[fZone.whichFace(face1)];
                }

                meshMod.setAction(polyRemoveFace(face0));
                meshMod.setAction
                (
                    polyModifyFace
                    (
                        faces[face1],           // modified face
                        face1,                  // label of face being modified
                        own1,                   // owner
                        own0,                   // neighbour
                        false,                  // face flip
                        -1,                     // patch for face
                        false,                  // remove from zone
                        zoneID,                 // zone for face
                        zoneFlip                // face flip in zone
                    )
                );
            }
        }
    }
}


label patchSize(const polyMesh& mesh, const labelList& patchIDs)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    label sz = 0;
    forAll(patchIDs, i)
    {
        const polyPatch& pp = patches[patchIDs[i]];
        sz += pp.size();
    }
    return sz;
}


labelList patchFaces(const polyMesh& mesh, const labelList& patchIDs)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    labelList faceIDs(patchSize(mesh, patchIDs));
    label sz = 0;
    forAll(patchIDs, i)
    {
        const polyPatch& pp = patches[patchIDs[i]];

        forAll(pp, ppi)
        {
            faceIDs[sz++] = pp.start()+ppi;
        }
    }

if (faceIDs.size() != sz)
{
    FatalErrorInFunction << exit(FatalError);
}

    return faceIDs;
}


labelList findBaffles(const polyMesh& mesh, const labelList& boundaryFaces)
{
    // Get all duplicate face labels (in boundaryFaces indices!).
    labelList duplicates = localPointRegion::findDuplicateFaces
    (
        mesh,
        boundaryFaces
    );


    // Check that none are on processor patches
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    forAll(duplicates, bFacei)
    {
        if (duplicates[bFacei] != -1)
        {
            label facei = boundaryFaces[bFacei];
            label patchi = patches.whichPatch(facei);

            if (isA<processorPolyPatch>(patches[patchi]))
            {
                FatalErrorInFunction
                    << "Duplicate face " << facei
                    << " is on a processorPolyPatch."
                    << "This is not allowed." << nl
                    << "Face:" << facei
                    << " is on patch:" << patches[patchi].name()
                    << abort(FatalError);
            }
        }
    }


    // Write to faceSet for ease of post-processing.
    
    {
        faceSet duplicateSet
        (
            mesh,
            "duplicateFaces",
            (mesh.nFaces() - mesh.nInternalFaces())/256
        );

        forAll(duplicates, bFacei)
        {
            label otherFacei = duplicates[bFacei];

            if (otherFacei != -1 && otherFacei > bFacei)
            {
                duplicateSet.insert(boundaryFaces[bFacei]);
                duplicateSet.insert(boundaryFaces[otherFacei]);
            }
        }

        Info<< "Writing " << returnReduce(duplicateSet.size(), sumOp<label>())
            << " duplicate faces to faceSet " << duplicateSet.objectPath()
            << nl << endl;
        duplicateSet.write();
    }
    
    return duplicates;
}


int removeBaffles(fvMesh& mesh, Time& runTime)//(int argc, char *argv[])
{
    Info<< "Creating fluidRegionNB mesh without baffles" << nl << endl;
    //argList::addNote
    //(
    //    "Detect faces that share points (baffles).\n"
    //    "Merge them or duplicate the points."
    //);

    //#include "addOverwriteOption.H"
    //#include "addRegionOption.H"
    //#include "addDictOption.H"
    /*
    argList::addBoolOption
    (
        "detectOnly",
        "find baffles only, but do not merge or split them"
    );
    argList::addBoolOption
    (
        "split",
        "topologically split duplicate surfaces"
    );
    */
    //#include "setRootCase.H"
    //#include "createTime.H"
    //Foam::Time runTime(Foam::Time::controlDictName, args);
    runTime.functionObjects().off();
    //#include "createNamedMesh.H"

    const word oldInstance = mesh.pointsInstance();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    const bool split      = 0;//args.found("split");
    const bool overwrite  = 1;//args.found("overwrite");
    const bool detectOnly = 0;//args.found("detectOnly");

    labelList detectPatchIDs;
    labelList splitPatchIDs;
    labelList mergePatchIDs;

    {
        if (detectOnly)
        {
            detectPatchIDs = identity(patches.size());
        }
        else if (split)
        {
            splitPatchIDs = identity(patches.size());
        }
        else
        {
            mergePatchIDs = identity(patches.size());
        }
    }


    if (detectPatchIDs.size())
    {
        findBaffles(mesh, patchFaces(mesh, detectPatchIDs));

        if (detectOnly)
        {
            return 0;
        }
    }

    // Read objects in time directory
    IOobjectList objects(mesh, runTime.timeName());

    // Read vol fields.

    PtrList<volScalarField> vsFlds;
    ReadFields(mesh, objects, vsFlds);

    PtrList<volVectorField> vvFlds;
    ReadFields(mesh, objects, vvFlds);

    PtrList<volSphericalTensorField> vstFlds;
    ReadFields(mesh, objects, vstFlds);

    PtrList<volSymmTensorField> vsymtFlds;
    ReadFields(mesh, objects, vsymtFlds);

    PtrList<volTensorField> vtFlds;
    ReadFields(mesh, objects, vtFlds);

    // Read surface fields.

    PtrList<surfaceScalarField> ssFlds;
    ReadFields(mesh, objects, ssFlds);

    PtrList<surfaceVectorField> svFlds;
    ReadFields(mesh, objects, svFlds);

    PtrList<surfaceSphericalTensorField> sstFlds;
    ReadFields(mesh, objects, sstFlds);

    PtrList<surfaceSymmTensorField> ssymtFlds;
    ReadFields(mesh, objects, ssymtFlds);

    PtrList<surfaceTensorField> stFlds;
    ReadFields(mesh, objects, stFlds);



    if (mergePatchIDs.size())
    {
        Info<< "Merging duplicate faces" << nl << endl;

        // Mesh change engine
        polyTopoChange meshMod(mesh);

        const labelList boundaryFaces(patchFaces(mesh, mergePatchIDs));

        // Get all duplicate face pairs (in boundaryFaces indices!).
        labelList duplicates(findBaffles(mesh, boundaryFaces));

        // Merge into internal faces.
        insertDuplicateMerge(mesh, boundaryFaces, duplicates, meshMod);

        if (!overwrite)
        {
            runTime++;
        }

        // Change the mesh. No inflation.
        autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false);

        // Update fields
        mesh.updateMesh(map());

        // Move mesh (since morphing does not do this)
        if (map().hasMotionPoints())
        {
            mesh.movePoints(map().preMotionPoints());
        }

        if (overwrite)
        {
            mesh.setInstance(oldInstance);
        }
        Info<< "Writing mesh to time " << runTime.timeName() << endl;
        mesh.write();
    }


    if (splitPatchIDs.size())
    {
        Info<< "Topologically splitting duplicate surfaces"
            << ", i.e. duplicating points internal to duplicate surfaces"
            << nl << endl;

        // Determine points on split patches
        DynamicList<label> candidates;
        {
            label sz = 0;
            forAll(splitPatchIDs, i)
            {
                sz += patches[splitPatchIDs[i]].nPoints();
            }
            candidates.setCapacity(sz);

            bitSet isCandidate(mesh.nPoints());
            forAll(splitPatchIDs, i)
            {
                const labelList& mp = patches[splitPatchIDs[i]].meshPoints();
                forAll(mp, mpi)
                {
                    label pointi = mp[mpi];
                    if (isCandidate.set(pointi))
                    {
                        candidates.append(pointi);
                    }
                }
            }
        }


        // Analyse which points need to be duplicated
        localPointRegion regionSide(mesh, candidates);

        // Point duplication engine
        duplicatePoints pointDuplicator(mesh);

        // Mesh change engine
        polyTopoChange meshMod(mesh);

        // Insert topo changes
        pointDuplicator.setRefinement(regionSide, meshMod);

        if (!overwrite)
        {
            runTime++;
        }

        // Change the mesh. No inflation.
        autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false);

        // Update fields
        mesh.updateMesh(map());

        // Move mesh (since morphing does not do this)
        if (map().hasMotionPoints())
        {
            mesh.movePoints(map().preMotionPoints());
        }

        if (overwrite)
        {
            mesh.setInstance(oldInstance);
        }
        Info<< "Writing mesh to time " << runTime.timeName() << endl;
        mesh.write();

        topoSet::removeFiles(mesh);
        processorMeshes::removeFiles(mesh);
        
        //Dump duplicated points (if any)
        const labelList& pointMap = map().pointMap();

        labelList nDupPerPoint(map().nOldPoints(), 0);

        pointSet dupPoints(mesh, "duplicatedPoints", 100);

        forAll(pointMap, pointi)
        {
            label oldPointi = pointMap[pointi];

            nDupPerPoint[oldPointi]++;

            if (nDupPerPoint[oldPointi] > 1)
            {
                dupPoints.insert(map().reversePointMap()[oldPointi]);
                dupPoints.insert(pointi);
            }
        }

        Info<< "Writing " << returnReduce(dupPoints.size(), sumOp<label>())
            << " duplicated points to pointSet "
            << dupPoints.objectPath() << nl << endl;

        dupPoints.write();
        
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
