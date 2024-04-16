// Gmsh project created on Thu Nov 23 08:39:57 2023
SetFactory("OpenCASCADE");
a = 0.015;
//+a
Point(1) = {0.0, 0, 0, a};
//+
Point(2) = {5, 0, 0, a};
//+
Point(3) = {5, 5, 0, a};
//+
Point(4) = {0.0, 5, 0, a*0.05};
//+
Point(5) = {.5, .0, 0, a*0.01};
//+
Point(6) = {.5, 5, 0, a*0.05};

//+
Line(1) = {4, 6};
//+
Line(2) = {6, 3};
//+
Line(3) = {3, 2};
//+
Line(4) = {2, 5};
//+
Line(5) = {5, 1};
//+
Line(6) = {1, 4};
//+
Line(7) = {6, 5};

Transfinite Curve {7,3} = 21 Using Progression 1.3;
Transfinite Curve {6} = 21 Using Progression 1/1.3;
Transfinite Curve {5} = 11 Using Progression 1.3;
Transfinite Curve {1} = 11 Using Progression 1/	1.3;
Transfinite Curve {2} = 16 Using Progression 1.40;
Transfinite Curve {4} = 16 Using Progression 1/1.40;
//+
Curve Loop(1) = {6, 1, 7, 5};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {2, 3, 4, -7}; 
//+
Plane Surface(2) = {2};

Transfinite Surface{1,2};
Recombine Surface {1,2};
//+
Physical Curve("Load", 8) = {1};
//+
Physical Curve("Symmetry", 9) = {6};
//+
Physical Curve("Bottom", 10) = {5, 4};
//+
Physical Curve("Right", 11) = {3};
//+
Physical Surface("Domain", 12) = {1, 2};
//+
Physical Curve("Top", 13) = {2};
