// Gmsh project created on Sun Nov  5 15:11:54 2023
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1000.0};
//+
Point(2) = {-173.205080757/100, 0, 0, 1000.0};
Point(3) = {173.205080757/100, 0, 0, 1000.0};
Point(4) = {0, -100./100, 0, 1000.0};

//+
//Line(1) = {2, 4};
//+
Line(2) = {4, 1};
//+
//Line(3) = {4, 3};
//+
Physical Curve("Domain", 4) = {2};
//+
Physical Point("Fixed", 5) = {1};
//+
Physical Point("Disp", 6) = {4};
