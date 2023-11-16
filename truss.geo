// Gmsh project created on Sun Nov  5 15:11:54 2023
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Physical Point("Left", 2) = {1};
//+
Physical Point("Right", 3) = {2};
//+
Physical Curve("Domain", 4) = {1};
