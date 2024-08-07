// Gmsh project created on Wed Aug  7 09:33:58 2024
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {0, 1, 0, 1.0};
//+
Point(3) = {1, 1, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Physical Point("apoios", 3) = {1, 3};
//+
Physical Curve("vertical", 4) = {1};
//+
Physical Curve("horizontal", 5) = {2};
