
Point(1) = {-halfEdge, -halfFlatToFlat,  -1, lc} ;
Point(2) = {-edge, 0, -1, lc} ;
Point(3) = {-halfEdge, halfFlatToFlat,  -1, lc} ;
Point(4) = {halfEdge, halfFlatToFlat,  -1, lc} ;
Point(5) = {edge, 0, -1, lc} ;
Point(6) = {halfEdge, -halfFlatToFlat,  -1, lc} ;

Line(1) = {1,2} ;
Line(2) = {2,3} ;
Line(3) = {3,4} ;
Line(4) = {4,5} ;
Line(5) = {5,6} ;
Line(6) = {6,1} ;

Line Loop(7) = {1,2,3,4,5,6} ;

Plane Surface(1) = {7} ;


If(!forNeutronics)
  Point(7)  = {-halfEdgeExt, -halfFlatToFlatExt,  -1, lc} ;
  Point(8)  = {-edgeExt, 0, -1, lc} ;
  Point(9)  = {-halfEdgeExt, halfFlatToFlatExt,  -1, lc} ;
  Point(10) = {halfEdgeExt, halfFlatToFlatExt,  -1, lc} ;
  Point(11) = {edgeExt, 0, -1, lc} ;
  Point(12) = {halfEdgeExt, -halfFlatToFlatExt,  -1, lc} ;

  Line(8) = {7,8} ;
  Line(9) = {8,9} ;
  Line(10) = {9,10} ;
  Line(11) = {10,11} ;
  Line(12) = {11,12} ;
  Line(13) = {12,7} ;

  Line Loop(14) = {8,9,10,11,12,13} ;

Plane Surface(2) = {7,-14} ;
EndIf

