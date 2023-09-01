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

Class
    Foam::insertDuplicateMerge
    
Description
    Detects boundary faces that share points (baffles). Either merges them or
    duplicate the points.

Note
    - can only handle pairwise boundary faces. So three faces using
      the same points is not handled (is illegal mesh anyway)

    - surfaces consisting of duplicate faces can be topologically split
    if the points on the interior of the surface cannot walk to all the
    cells that use them in one go.

    - Parallel operation (where duplicate face is perpendicular to a coupled
    boundary) is supported but not really tested.
    (Note that coupled faces themselves are not seen as duplicate faces)

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
);


label patchSize(const polyMesh& mesh, const labelList& patchIDs);


labelList patchFaces(const polyMesh& mesh, const labelList& patchIDs);


labelList findBaffles(const polyMesh& mesh, const labelList& boundaryFaces);


int removeBaffles(fvMesh& mesh, Time& runTime);//(int argc, char *argv[])


// ************************************************************************* //
