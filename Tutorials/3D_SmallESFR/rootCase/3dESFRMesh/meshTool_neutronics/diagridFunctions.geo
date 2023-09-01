
Function insertDiagrid
  surface[] = Translate {X, Y, 1} { Duplicata{ Surface{1}; } };
  If(!forTM)
    bottomSurfaces +=  surface[0];
  EndIf
  assembly[] = Extrude {0, 0, diagridHeight} { 
  Surface{surface[0]};  Layers{diagridNodes} ; Recombine;
  };
  If(!forNeutronics)
    wallSurfaces += assembly[2];
    wallSurfaces += assembly[3];
    wallSurfaces += assembly[4];
    wallSurfaces += assembly[5];
    wallSurfaces += assembly[6];
    wallSurfaces += assembly[7];

    surfaceGap[] = Translate {X, Y, 1} { Duplicata{ Surface{2}; } };
    If(!forTM)
      bottomGapSurfaces +=  surfaceGap[0];
    EndIf
    assemblyGap[] = Extrude {0, 0, diagridHeight} { 
    Surface{surfaceGap[0]};  Layers{diagridNodes} ; Recombine;
    };
    assemblyGapVolumes += assemblyGap[1];
  EndIf

  diagridVolumes += assembly[1];
  
Return
