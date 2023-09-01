
Function insertUpperFGP
  surface[] = Translate {X, Y, 1+lowReflHeight+lowFGPHeight+diagridHeight+coreHeight} { Duplicata{ Surface{1}; } };

  assembly[] = Extrude {0, 0, upperFGPHeight} { 
  Surface{surface[0]};  Layers{upperFGPNodes} ; Recombine;
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
    	assemblyGap[] = Extrude {0, 0, upperReflHeight} { 
    	Surface{surfaceGap[0]};  Layers{upperReflNodes} ; Recombine;
    	};
    	assemblyGapVolumes += assemblyGap[1];
    EndIf
  EndIf


  isCR = 0 ; 
  For l In {0 : (#CRposX[] - 1)}
    If(((X-CRposX[l])*(X-CRposX[l]))<0.01 && ((Y-CRposY[l])*(Y-CRposY[l]))<0.01)
      isCR = 1;
    EndIf
  EndFor

  isO = 0;
  For l In {0 : (#OposX[] - 1)}
    If(((X-OposX[l])*(X-OposX[l]))<0.01 && ((Y-OposY[l])*(Y-OposY[l]))<0.01)
      isO = 1;
    EndIf
  EndFor

  isR = 0;
  For l In {0 : (#RposX[] - 1)}
    If(((X-RposX[l])*(X-RposX[l]))<0.01 && ((Y-RposY[l])*(Y-RposY[l]))<0.01)
      isR = 1;
    EndIf
  EndFor
    If(j>=(innerRows+outerRows) && !isO)
      isR = 1;
    EndIf
  
  If(isCR)
    controlAssemblyVolumes += assembly[1];
  EndIf

  If(isR)
    reflectorAssemblyVolumes += assembly[1];
  EndIf

  If(!isCR  && !isR)
    upperFGPVolumes += assembly[1];
  EndIf

Return
