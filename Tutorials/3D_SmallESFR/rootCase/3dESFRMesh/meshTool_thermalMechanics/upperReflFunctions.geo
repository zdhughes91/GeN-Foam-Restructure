
Function insertUpperRefl
  surface[] = Translate {X, Y, 1+lowReflHeight+lowFGPHeight+diagridHeight+coreHeight+upperFGPHeight} { Duplicata{ Surface{1}; } };
  assembly[] = Extrude {0, 0, upperReflHeight} { 
  Surface{surface[0]};  Layers{upperReflNodes} ; Recombine;
  };

  isCR = 0 ; 
  For l In {0 : (#CRposX[] - 1)}
    If(((X-CRposX[l])*(X-CRposX[l]))<0.01 && ((Y-CRposY[l])*(Y-CRposY[l]))<0.01)
      isCR = 1;
    EndIf
  EndFor


  If(isCR)
    topSurfacesCR +=  assembly[0];
  EndIf
  If(!isCR)
    topSurfaces +=  assembly[0];
  EndIf  



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
    upperReflVolumes += assembly[1];
  EndIf

Return
