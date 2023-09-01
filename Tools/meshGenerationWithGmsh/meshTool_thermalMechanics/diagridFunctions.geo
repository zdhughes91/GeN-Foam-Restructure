
Function insertDiagrid
  surface[] = Translate {X, Y, 1} { Duplicata{ Surface{1}; } };
  If(!forTM)
    bottomSurfaces +=  surface[0];
  EndIf
  assembly[] = Extrude {0, 0, diagridHeight} { 
  Surface{surface[0]};  Layers{diagridNodes} ; Recombine;
  };

  isCR = 0 ; 
  For l In {0 : (#CRposX[] - 1)}
    If(((X-CRposX[l])*(X-CRposX[l]))<0.01 && ((Y-CRposY[l])*(Y-CRposY[l]))<0.01)
      isCR = 1;
    EndIf
  EndFor

  If(isCR)
    bottomSurfacesCR +=  surface[0];
  EndIf
  If(!isCR)
    bottomSurfaces +=  surface[0];
  EndIf


  diagridVolumes += assembly[1];
  
Return
