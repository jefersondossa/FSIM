// Gmsh project created on Thu Nov 23 08:39:57 2023
SetFactory("OpenCASCADE");
a = 0.01;
alturageogrelha = 0.1;
//+
Point(1) = {0, 0, 0, a};
//+
Point(2) = {.15, 0, 0, a};
//+
Point(3) = {.15, .15, 0, a};
//+
Point(4) = {0, .15, 0, a};
//+
Point(5) = {.055, .15, 0, a};
//+
Point(6) = {.095, .15, 0, a};

b = 0.005;
Point(7) = {.0, alturageogrelha, 0, a};
Point(8) = {.025, alturageogrelha, 0, b};
Point(9) = {.125, alturageogrelha, 0, b};
Point(10) = {.15, alturageogrelha, 0, a};

//+
Line(1) = {4, 5};
//+
Line(2) = {5, 6};
//+
Line(3) = {6, 3};
//+
Line(4) = {3, 10};
//+
Line(5) = {10, 2};
//+
Line(6) = {2, 1};
//+
Line(7) = {1, 7};
//+
Line(8) = {7, 4};
//+
Line(9) = {7, 8};
//+
Line(10) = {8, 9};
//+
Line(11) = {9, 10};
//+
Curve Loop(1) = {8, 1, 2, 3, 4, -11, -10, -9};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {7, 9, 10, 11, 5, 6};
//+
Plane Surface(2) = {2};
//+
Physical Curve("Fundo", 12) = {6};
//+
Physical Curve("Laterais", 13) = {8, 7, 4, 5};
//+
Physical Curve("Topo", 14) = {1, 3};
//+
Physical Curve("Carga", 15) = {2};
//+
Physical Curve("Geogrelha", 16) = {10};
//+
Physical Surface("Solo", 17) = {1, 2};
