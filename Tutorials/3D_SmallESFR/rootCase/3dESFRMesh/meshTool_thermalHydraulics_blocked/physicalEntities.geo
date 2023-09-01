
Physical Surface("bottom") = {bottomSurfaces[]} ;
Physical Surface("blockedAssembly") = {blockedSurface[]} ;

Physical Surface("walls") = {wallSurfaces[]} ;
Physical Surface("top") = {topSurfaces[]} ;
Physical Surface("zeroDispl") = {zeroDisplSurfaces[]} ;

Physical Volume("innerCore") = {innerAssemblyVolumes[]} ;
Physical Volume("outerCore") = {outerAssemblyVolumes[]} ;
Physical Volume("radialReflector") = {reflectorAssemblyVolumes[]} ;
Physical Volume("controlRod") = {controlAssemblyVolumes[]} ;

Physical Volume("follower") = {followerVolumes[]} ;

Physical Volume("axialReflector") = {upperFGPVolumes[], upperReflVolumes[], lowReflVolumes[], lowFGPVolumes[]} ;
/*
Physical Volume("axialReflector") = {upperReflVolumes[]} ;
Physical Volume("axialReflector") = {lowReflVolumes[]} ;
Physical Volume("axialReflector") = {lowFGPVolumes[]} ;
*/
Physical Volume("diagrid") = {diagridVolumes[]} ;
Physical Volume("softStructure") = {softStructureVolumes[]} ;


Physical Surface("bottomGap") = {bottomGapSurfaces[]} ;
Physical Surface("topGap") = {topGapSurfaces[]} ;
Physical Volume("gap") = {assemblyGapVolumes[]} ;

