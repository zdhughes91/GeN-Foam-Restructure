// Unfortunately, meshToMesh projection can have problems in the presence of baffles.
// As a workaround, we allow (via the keyword "removeBaffles") the creation of a 
// second thermalHydraulic mesh (fluidRegionNB) without baffles and over which
// the original fields are copied cell by cell. If present, fluidRegionNB will have
// its own set of fields, distinguished by the "NB" (No Baffle) expression in the name

const word fluidRegionName("fluidRegion");
Foam::fvMesh fluidMesh
(
    Foam::IOobject
    (
        fluidRegionName,
        runTime.timeName(),
        runTime,
        Foam::IOobject::MUST_READ
    )
);

autoPtr<fvMesh> fluidMeshNBPtr;

if(runTime.controlDict().lookupOrDefault<bool>("removeBaffles", false))
{
    if (UPstream::nProcs() > 1)
    {
        Info<< "WARNING: the removeBaffles feature is guaranteed to work only "
            << "if each processor contains complete pairs of master and slave "
            << "faces of baffe patches! Otherwise, projection artifacts can "
            << " still occur" << nl << endl;
    }
    
    //- Copy internal geometric data from fluidMesh
    pointField points(fluidMesh.points());
    faceList faces(fluidMesh.faces());
    cellList cells(fluidMesh.cells());

    //- Copy boundary geometric data from fluidMesh
    const polyBoundaryMesh& patches = fluidMesh.boundaryMesh();
    wordList patchNames(patches.names());
    List<polyPatch*> pList;
    forAll(patchNames, i)
    {
        pList.append
        (
            new polyPatch
            (
                patches[patchNames[i]],
                patches
            )
        );
    }

    //- Assemble mesh as copy from points, faces, cells with no boundary
    fluidMeshNBPtr.reset
    (
        new Foam::fvMesh
        (
            Foam::IOobject
            (
                "fluidRegionNB",
                runTime.timeName(),
                runTime,
                Foam::IOobject::NO_READ
            ),
            std::move(points),
            std::move(faces),
            std::move(cells)
        )
    );

    //- Add boundary data
    fluidMeshNBPtr().addFvPatches(pList);

    //- One of the main issues when passing the actual runTime to removeBaffles
    //  was that the functionObjects execution flags were getting reset to
    //  false. Probably other aspects of runTime were getting modified too. To
    //  avoid all these issues, let removeBaffles operate with a dummy runTime
    //  created on the fly, which is deleted afterwards
    autoPtr<Time> dummyRunTimePtr
    (
        new Time(Foam::Time::controlDictName, args)
    );
    removeBaffles(fluidMeshNBPtr(), dummyRunTimePtr());
    dummyRunTimePtr.clear();    

    //- Copy mesh zones
    List<pointZone*> pointZonesNB(0);
    List<faceZone*> faceZonesNB(0);
    List<cellZone*> cellZonesNB(0);
    forAll(fluidMesh.pointZones().names(), i)
    {
        word name(fluidMesh.pointZones().names()[i]);
        const pointZone& origZone(fluidMesh.pointZones()[name]);
        pointZonesNB.append
        (
            new pointZone
            (
                name,
                origZone,
                i,
                fluidMesh.pointZones()  //- Not sure whether I should use
                                        //  fluidMesh or fluidMeshNBPtr() here,
                                        //  both seem to result in the same
                                        //  zones getting assembled, so, not
                                        //  too big of a deal but this should
                                        //  be checked eventually
            )
        );
    }
    forAll(fluidMesh.faceZones().names(), i)
    {
        word name(fluidMesh.faceZones().names()[i]);
        const faceZone& origZone(fluidMesh.faceZones()[name]);
        faceZonesNB.append
        (
            new faceZone
            (
                name,
                origZone,
                origZone.flipMap(),
                i,
                fluidMesh.faceZones()
            )
        );
    }
    forAll(fluidMesh.cellZones().names(), i)
    {
        word name(fluidMesh.cellZones().names()[i]);
        const cellZone& origZone(fluidMesh.cellZones()[name]);
        cellZonesNB.append
        (
            new cellZone
            (
                name,
                origZone,
                i,
                fluidMesh.cellZones()
            )
        );
    }
    fluidMeshNBPtr().addZones(pointZonesNB, faceZonesNB, cellZonesNB);
}
else
{
    fluidMeshNBPtr.reset
    (
        new Foam::fvMesh
        (
            Foam::IOobject
            (
                "fluidRegion",
                runTime.timeName(),
                runTime,
                Foam::IOobject::MUST_READ
            )
        )
    );  
}

fvMesh& fluidMeshNB = fluidMeshNBPtr();

/*
Info << "Print NB zones" << endl;
forAllIter
(
    pointZoneMesh,
    fluidMeshNBPtr().pointZones(),
    iter
)
{
    Info << (*iter).name() << endl;
    Info << (*iter) << endl;
}
forAllIter
(
    faceZoneMesh,
    fluidMeshNBPtr().faceZones(),
    iter
)
{
    Info << (*iter).name() << endl;
    Info << (*iter) << endl;
}
forAllIter
(
    cellZoneMesh,
    fluidMeshNBPtr().cellZones(),
    iter
)
{
    Info << (*iter).name() << endl;
    Info << (*iter) << endl;
}
Info << "Done printing NB zones" << endl;
*/

//- Silence "unused variable" warnings, I actually do use it and have no idea
//  why the compiler does not "see" it
(void) fluidMeshNB;
