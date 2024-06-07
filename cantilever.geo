// Gmsh project created on Wed May 15 10:19:39 2024
//+

//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Point(3) = {1, .05, 0, 1.0};
//+
Point(4) = {0, .05, 0, 1.0};
//+
Point(5) = {0, .1, 0, 1.0};
//+
Point(6) = {1, .1, 0, 1.0};
//+
Line(1) = {5, 4};
//+
Line(2) = {4, 1};
//+
Line(3) = {1, 2};
//+
Line(4) = {2, 3};
//+
Line(5) = {3, 6};
//+
Line(6) = {6, 5};
//+
Line(7) = {4, 3};
//+
Curve Loop(1) = {6, 1, 7, 5};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {2, 3, 4, -7};
//+
Plane Surface(2) = {-2};
//+
Physical Point("ponto engaste", 8) = {4};
//+
Physical Point("ponto forca", 9) = {3};
//+
Physical Curve("engaste", 10) = {1, 2};
//+
Physical Surface("dominio", 11) = {1, 2};
//+
Transfinite Curve {1, 2, 5, 4} = 2 Using Progression 1;
//+
Transfinite Curve {6, 7, 3} = 51 Using Progression 1;
//+
Recombine Surface {1, 2};
//+
Transfinite Surface {1};
//+
Transfinite Surface {2};
