


For j In {0 : (innerRows-1)}
  x = 0;
  y = j*pitch;
  X = x * Sqrt(3)/2;
  Y = x * 1/2 + y;
  Call insertInnerAssembly;
  If (j>0)
  For k In {1:5}
    X0 = 0.5*X - Sin(Pi/3)*Y;
    Y0 = Sin(Pi/3)*X + 0.5*Y;
    X = X0;
    Y = Y0;
    Call insertInnerAssembly;      
  EndFor
  
  For i In {1:(j-1)}
    x = i*pitch;
    y -= pitch;
    X = x * Sqrt(3)/2;
    Y = x * 1/2 + y;
    Call insertInnerAssembly;
    For k In {1:5}
      X0 = 0.5*X - Sin(Pi/3)*Y;
      Y0 = Sin(Pi/3)*X + 0.5*Y;
      X = X0;
      Y = Y0;
      Call insertInnerAssembly;      
    EndFor
  EndFor
  EndIf

EndFor


For j In {innerRows : (innerRows+outerRows-1)}
  x = 0;
  y = j*pitch;
  X = x * Sqrt(3)/2;
  Y = x * 1/2 + y;
  Call insertOuterAssembly;
  For k In {1:5}
    X0 = 0.5*X - Sin(Pi/3)*Y;
    Y0 = Sin(Pi/3)*X + 0.5*Y;
    X = X0;
    Y = Y0;
    Call insertOuterAssembly;      
  EndFor
  For i In {1:(j-1)}
    x = i*pitch;
    y -= pitch;
    X = x * Sqrt(3)/2;
    Y = x * 1/2 + y;
    Call insertOuterAssembly;
    For k In {1:5}
      X0 = 0.5*X - Sin(Pi/3)*Y;
      Y0 = Sin(Pi/3)*X + 0.5*Y;
      X = X0;
      Y = Y0;
      Call insertOuterAssembly;      
    EndFor
  EndFor
EndFor


For j In {(innerRows+outerRows) : (innerRows+outerRows+reflectorRows-1)}
  x = 0;
  y = j*pitch;
  X = x * Sqrt(3)/2;
  Y = x * 1/2 + y;
  Call insertReflectorAssembly;
  For k In {1:5}
    X0 = 0.5*X - Sin(Pi/3)*Y;
    Y0 = Sin(Pi/3)*X + 0.5*Y;
    X = X0;
    Y = Y0;
    Call insertReflectorAssembly;      
  EndFor
  For i In {1:(j-1)}
    x = i*pitch;
    y -= pitch;
    X = x * Sqrt(3)/2;
    Y = x * 1/2 + y;
    Call insertReflectorAssembly;
    For k In {1:5}
      X0 = 0.5*X - Sin(Pi/3)*Y;
      Y0 = Sin(Pi/3)*X + 0.5*Y;
      X = X0;
      Y = Y0;
      Call insertReflectorAssembly;      
    EndFor
  EndFor
EndFor






