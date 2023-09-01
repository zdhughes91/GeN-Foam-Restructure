  If(forTM)
    surface[] = Translate {0, 0, 1-0.2} { Duplicata{ Surface{1}; } };
    zeroDisplSurfaces +=  surface[0];
    assembly[] = Extrude {0, 0, 0.2} { 
    Surface{surface[0]};  Layers{3} ; Recombine;
    };
    softStructureVolumes += assembly[1];
  EndIf

