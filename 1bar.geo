// Gmsh project created on Wed Jul 31 10:11:55 2024
//+
Point(1) = {-1, 0, 0, 2.0};
//+
Point(2) = {1, 0, 0, 2.0};
//+
Line(1) = {1, 2};
//+
Physical Point("engaste", 3) = {1};
//+
Physical Point("forca", 4) = {2};
//+
Physical Curve("viga", 5) = {1};
