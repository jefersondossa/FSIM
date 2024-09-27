// Gmsh project created on Mon Aug  5 15:43:37 2024
//+
Point(1) = { 0.5*Sqrt(Pi/4),  0.5*Sqrt(Pi/4), 0, 1.0};
Point(2) = {-0.5*Sqrt(Pi/4),  0.5*Sqrt(Pi/4), 0, 1.0};
Point(3) = {-0.5*Sqrt(Pi/4), -0.5*Sqrt(Pi/4), 0, 1.0};
Point(4) = { 0.5*Sqrt(Pi/4), -0.5*Sqrt(Pi/4), 0, 1.0};
Point(5) = {0., 0., 0, 1.0};
//+
Circle(1) = {1, 5, 2};
//+
Circle(2) = {2, 5, 3};
//+
Circle(3) = {3, 5, 4};
//+
Circle(4) = {4, 5, 1};

in=15;
ou=40;
up=15;
Point(6) = {-in, up, 0, 1.0};
Point(7) = {-in, -up, 0, 1.0};
Point(8) = { ou, up, 0, 1.0};
Point(9) = { ou, -up, 0, 1.0};
//+
Line(5) = {6, 7};
//+
Line(6) = {7, 3};
//+
Line(7) = {4, 9};
//+
Line(8) = {9, 7};
//+
Line(9) = {1, 8};
//+
Line(10) = {8, 9};
//+
Line(11) = {2, 6};
//+
Line(12) = {6, 8};
//+
Curve Loop(1) = {11, 5, 6, -2};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {6, 3, 7, 8};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {7, -10, -9, -4};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {9, -12, -11, -1};
//+
Plane Surface(4) = {4};

//+
Extrude {0, 0, 1} {
  Surface{4}; Surface{1}; Surface{2}; Surface{3}; 
}
//+
Transfinite Curve {24, 28, 19, 46, 72, 68, 50, 20} = 2 Using Progression 1;
//+
cyl = 10;
diag = 30;
transval=1.05;
Transfinite Curve {17, 1, 2, 39, 59, 83, 4, 3, 61, 8, 81, 10, 12, 15, 5, 37} = cyl Using Progression 1;
//+

Transfinite Curve {11} = diag Using Progression transval;
Transfinite Curve {16} = diag Using Progression 1/transval;
Transfinite Curve {6, 38} = diag Using Progression 1/transval;
//+
Transfinite Curve {9, 14} = diag Using Progression transval;
Transfinite Curve {7, 60} = diag Using Progression transval;
//+
Transfinite Surface {21};
//+
Transfinite Surface {25};
//+
Transfinite Surface {29};
//+
Transfinite Surface {33};
//+
Transfinite Surface {4};
//+
Transfinite Surface {34};
//+
Transfinite Surface {56};
//+
Transfinite Surface {1};
//+
Transfinite Surface {55};
//+
Transfinite Surface {47};
//+
Transfinite Surface {51};
//+
Transfinite Surface {77};
//+
Transfinite Surface {73};
//+
Transfinite Surface {91};
//+
Transfinite Surface {100};
//+
Transfinite Surface {3};
//+
Transfinite Surface {78};
//+
Transfinite Surface {2};
//+
Transfinite Surface {69};
//+
Transfinite Surface {99};

Recombine Surface {99, 69, 2, 78, 3, 100, 91, 73, 77, 51, 47, 55, 1, 56, 34, 4, 33, 29 ,25, 21};
//+
Transfinite Volume{2};
//+
Transfinite Volume{1};
//+
Transfinite Volume{4};
//+
Transfinite Volume{3};
//+
Physical Surface("cylinder", 101) = {33, 55, 69, 99};
//+
Physical Surface("inlet", 102) = {47};
//+
Physical Surface("outlet", 103) = {91};
//+
Physical Surface("symmetry", 104) = {77, 25, 1, 56, 2, 78, 3, 100, 4, 34};
//+
Physical Volume("fluid", 105) = {2, 1, 4, 3};
