


For j In {0 : (innerRows+outerRows+reflectorRows-1)}
  x = 0;
  y = j*pitch;
  X = x * Sqrt(3)/2;
  Y = x * 1/2 + y;
  Call insertUpperFGP;
  If (j>0)
  For k In {1:5}
    X0 = 0.5*X - Sin(Pi/3)*Y;
    Y0 = Sin(Pi/3)*X + 0.5*Y;
    X = X0;
    Y = Y0;
    Call insertUpperFGP;      
  EndFor
  
  For i In {1:(j-1)}
    x = i*pitch;
    y -= pitch;
    X = x * Sqrt(3)/2;
    Y = x * 1/2 + y;
    Call insertUpperFGP;
    For k In {1:5}
      X0 = 0.5*X - Sin(Pi/3)*Y;
      Y0 = Sin(Pi/3)*X + 0.5*Y;
      X = X0;
      Y = Y0;
      Call insertUpperFGP;      
    EndFor
  EndFor
  EndIf

EndFor





