
forNeutronics = 1;
forTM = 1;
Sym60 = 1;

lc = 0.5;
edge = 0.12171 ; //0.11391
gap = 0; //0.009

coreHeight = 1;
lowReflHeight = 0.3;
lowFGPHeight = 0.3;
diagridHeight = 0.3;
upperFGPHeight = 0.11;
upperReflHeight = 0.3;

coreNodes = 10;
lowReflNodes = 3;
lowFGPNodes = 3;
diagridNodes = 3;
upperFGPNodes = 1;
upperReflNodes = 3;

innerRows = 5;
outerRows = 2;
reflectorRows = 1;
assemblyBalance = 0;

//CR position
CRposxn = {};
CRposyn = {};

//CRposxn += {0};
//CRposyn += {0};

CRposxn += {3};
CRposyn += {0};

CRposxn += {5};
CRposyn += {1};

CRposxn += {6};
CRposyn += {-1};

CRposxn += {7};
CRposyn += {2};

CRposxn += {9};
CRposyn += {-2};

CRposxn += {5};
CRposyn += {5};

//Inner position
Iposxn = {};
Iposyn = {};

//Iposxn += {5};
//Iposyn += {4};

//Iposxn += {6};
//Iposyn += {3};

//Iposxn += {9};
//Iposyn += {-3};

//Iposxn += {9};
//Iposyn += {-4};

//Outer position
Oposxn = {};
Oposyn = {};

//Oposxn += {7};
//Oposyn += {6};

//Oposxn += {8};
//Oposyn += {5};

//Oposxn += {13};
//Oposyn += {-5};

//Oposxn += {13};
//Oposyn += {-6};


//Reflector position
Rposxn = {};
Rposyn = {};

//Rposxn += {12};
//Rposyn += {0};

totalHeight =  lowFGPHeight + lowReflHeight + coreHeight + upperFGPHeight + upperReflHeight  ;

 
halfEdge = edge/2;
flatToFlat = edge * Sqrt(3);
halfFlatToFlat = flatToFlat/2;
pitch = flatToFlat + gap;

edgeExt = edge + gap/Sqrt(3);
halfEdgeExt = edgeExt/2;
flatToFlatExt = edgeExt * Sqrt(3);
halfFlatToFlatExt = flatToFlatExt/2;

totalRows = innerRows + outerRows;
totalAssemblyNumber = 1;
For j In {1 : (totalRows-1)}
  totalAssemblyNumber += 6*j;
EndFor
totalAssemblyNumber += assemblyBalance;

Iposx = {};
Iposy = {};
IposX = {};
IposY = {};
For i In {0 : (#Iposxn[]-1)}
  Iposx += {Iposxn[i]*pitch};
  Iposy += {Iposyn[i]*pitch};
  IposX += {Iposx[i] * Sqrt(3)/2};
  IposY += {Iposx[i] * 1/2 + Iposy[i]};

  If(Sym60)
    For kk In {1:5}
      IposX0 = 0.5*IposX[#IposX[]-1] - Sin(Pi/3)*IposY[#IposX[]-1];
      IposY0 = Sin(Pi/3)*IposX[#IposX[]-1] + 0.5*IposY[#IposX[]-1];
      IposX += IposX0;
      IposY += IposY0;
    EndFor    
  EndIf
EndFor

Oposx = {};
Oposy = {};
OposX = {};
OposY = {};
For i In {0 : (#Oposxn[]-1)}
  Oposx += {Oposxn[i]*pitch};
  Oposy += {Oposyn[i]*pitch};
  OposX += {Oposx[i] * Sqrt(3)/2};
  OposY += {Oposx[i] * 1/2 + Oposy[i]};

  If(Sym60)
    For kk In {1:5}
      OposX0 = 0.5*OposX[#OposX[]-1] - Sin(Pi/3)*OposY[#OposX[]-1];
      OposY0 = Sin(Pi/3)*OposX[#OposX[]-1] + 0.5*OposY[#OposX[]-1];
      OposX += OposX0;
      OposY += OposY0;
    EndFor    
  EndIf
EndFor

Rposx = {};
Rposy = {};
RposX = {};
RposY = {};
For i In {0 : (#Rposxn[]-1)}
  Rposx += {Rposxn[i]*pitch};
  Rposy += {Rposyn[i]*pitch};
  RposX += {Rposx[i] * Sqrt(3)/2};
  RposY += {Rposx[i] * 1/2 + Rposy[i]};

  If(Sym60)
    For kk In {1:5}
      RposX0 = 0.5*RposX[#RposX[]-1] - Sin(Pi/3)*RposY[#RposX[]-1];
      RposY0 = Sin(Pi/3)*RposX[#RposX[]-1] + 0.5*RposY[#RposX[]-1];
      RposX += RposX0;
      RposY += RposY0;
    EndFor    
  EndIf
EndFor

CRposx = {};
CRposy = {};
CRposX = {};
CRposY = {};
For i In {0 : (#CRposxn[]-1)}
  CRposx += {CRposxn[i]*pitch};
  CRposy += {CRposyn[i]*pitch};
  CRposX += {CRposx[i] * Sqrt(3)/2};
  CRposY += {CRposx[i] * 1/2 + CRposy[i]};

  If(Sym60)
    For kk In {1:5}
      CRposX0 = 0.5*CRposX[#CRposX[]-1] - Sin(Pi/3)*CRposY[#CRposX[]-1];
      CRposY0 = Sin(Pi/3)*CRposX[#CRposX[]-1] + 0.5*CRposY[#CRposX[]-1];
      CRposX += CRposX0;
      CRposY += CRposY0;
    EndFor    
  EndIf
EndFor




