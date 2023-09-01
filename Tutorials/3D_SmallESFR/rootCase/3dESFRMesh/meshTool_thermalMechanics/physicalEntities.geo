
Physical Surface("bottom") = {bottomSurfaces[]} ;
Physical Surface("bottomCR") = {bottomSurfacesCR[]} ;
//Physical Surface("walls") = {wallSurfaces[]} ;
Physical Surface("top") = {topSurfaces[]} ;
Physical Surface("topCR") = {topSurfacesCR[]} ;
Physical Surface("zeroDispl") = {zeroDisplSurfaces[]} ;

Physical Volume("innerCore") = {innerAssemblyVolumes[]} ;
Physical Volume("outerCore") = {outerAssemblyVolumes[]} ;
Physical Volume("radialReflector") = {reflectorAssemblyVolumes[]} ;

Physical Volume("follower") = {followerVolumes[]} ;

Physical Volume("controlRod") = {controlAssemblyVolumes[]} ;
Physical Volume("rest") = {upperFGPVolumes[], upperReflVolumes[], lowReflVolumes[], lowFGPVolumes[]} ;
/*
Physical Volume("rest") = {followerVolumes[]} ;

Physical Volume("rest") = {upperFGPVolumes[]} ;
Physical Volume("rest") = {upperReflVolumes[]} ;
Physical Volume("rest") = {lowReflVolumes[]} ;
Physical Volume("rest") = {lowFGPVolumes[]} ;
*/
Physical Volume("diagrid") = {diagridVolumes[]} ;
Physical Volume("softStructure") = {softStructureVolumes[]} ;


//Physical Surface("bottomGap") = {bottomGapSurfaces[]} ;
//Physical Surface("topGap") = {topGapSurfaces[]} ;
//Physical Volume("gap") = {assemblyGapVolumes[]} ;

