// Gmsh project created on Wed Jul 31 10:11:55 2024
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Physical Point("Apoioesquerda", 2) = {1};
//+
Physical Point("Apoiodireita", 3) = {2};
//+
Physical Curve("viga", 4) = {1};
