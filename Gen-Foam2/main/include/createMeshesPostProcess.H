#include "createMeshes.H"

/*
This will soon be used for a multiregion expansion
if (!fluidRegions.size() && !neutroRegions.size())
{
    FatalErrorIn(args.executable())
        << "No region meshes present" << exit(FatalError);
}

fvMesh& mesh = fluidRegions.size() ? fluidRegions[0] : neutroRegions[0];
*/
fvMesh& mesh = fluidMesh;
