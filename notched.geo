// Gmsh project created on Sun Nov 26 21:56:34 2023
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 15, 0, 2.0};
//+
Point(2) = {5, 15, 0, 2.0};
//+
Point(3) = {5, 0, 0, 2.0};
//+
Point(4) = {0, 0, 0, .05};
Point(5) = {0.5, 0, 0, 0.005};
//+
//+
Line(1) = {4, 5};
Line(2) = {5, 3};
//+
Line(3) = {3, 2};
//+
Line(4) = {2, 1};
//+
Line(5) = {4, 1};
//+
Curve Loop(1) = {-5, 4, 3, 2,1};
//+
Plane Surface(1) = {-1};
//+
Physical Curve("Bottom", 6) = {1};
Physical Curve("Notch", 7) = {2};
//+
Physical Curve("Left", 8) = {5};
//+
Physical Curve("Right", 9) = {3};
//+
Physical Curve("Top", 10) = {4};
//+
Physical Surface("Domain", 11) = {1};
