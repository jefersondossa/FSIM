// Gmsh project created on Thu Nov 23 08:39:57 2023
SetFactory("OpenCASCADE");
a = 0.2;
//+
Point(1) = {0, 0, 0, a};
//+
Point(2) = {1, 0, 0, a};
//+
Point(3) = {1, 1, 0, a};
//+
Point(4) = {0, 1, 0, a};
//+
Point(5) = {.3, 1, 0, a};
//+
Point(6) = {.7, 1, 0, a};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 6};
//+
Line(4) = {6, 5};
//+
Line(5) = {5, 4};
//+
Line(6) = {4, 1};
//+
Curve Loop(1) = {5, 6, 1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Physical Curve("Fundo", 7) = {1};
//+
Physical Curve("Laterais", 8) = {6, 2};
//+
Physical Curve("Topo", 9) = {5, 3};
//+
Physical Curve("Carga", 10) = {4};
//+
Physical Surface("Dominio", 11) = {1};
