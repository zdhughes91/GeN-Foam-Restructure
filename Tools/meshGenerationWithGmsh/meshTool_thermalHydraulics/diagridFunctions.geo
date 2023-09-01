
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


    If(gap != 0.0)
    	surfaceGap[] = Translate {X, Y, 1+lowReflHeight+lowFGPHeight+diagridHeight+coreHeight+upperFGPHeight} { Duplicata{ Surface{2}; } };
    	If(!forTM)
      		bottomGapSurfaces +=  surfaceGap[0];
    	EndIf 
   	assemblyGap[] = Extrude {0, 0, upperReflHeight} { 
    	Surface{surfaceGap[0]};  Layers{upperReflNodes} ; Recombine;
    	};
    	assemblyGapVolumes += assemblyGap[1];
    EndIf

  EndIf

  diagridVolumes += assembly[1];
  
Return
