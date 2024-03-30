// Gmsh project created on Thu Nov 23 08:39:57 2023
SetFactory("OpenCASCADE");
a = 0.015;
alturageogrelha = 0.01;
k=0.4;
//+
Point(1) = {0.0, 0-k, 0, a};
//+
Point(2) = {0.075+k, 0-k, 0, a};
//+
Point(3) = {0.075+k, 0.10, 0, a};
//+
Point(4) = {0.0, 0.10, 0, a*0.05};
//+
//Point(5) = {.055, .1, 0, a*0.01};
//+
Point(6) = {.02, 0.10, 0, a*0.05};

b = 0.001;
Point(7) = {.0, 0.10-alturageogrelha, 0, b};
//Point(8) = {.025, alturageogrelha, 0, b};
Point(9) = {0.05, 0.10-alturageogrelha, 0, b};
Point(10) = {0.075+k, 0.10-alturageogrelha, 0, a};

//+
Line(1) = {4, 6};
//+
//Line(2) = {5, 6};
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
Line(9) = {7, 9};
//+
//Line(10) = {8, 9};
//+
Line(11) = {9, 10};
//+
Curve Loop(1) = {8, 1, 3, 4, -11, -9};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {7, 9, 11, 5, 6};
//+
Plane Surface(2) = {2};
//+
Physical Curve("Fundo", 12) = {6};
//+
Physical Curve("Laterais", 13) = {8, 7, 4, 5};
//+
Physical Curve("Topo", 14) = { 3};
//+
Physical Curve("Carga", 15) = {1};
//+
Physical Curve("Geogrelha", 16) = {9};
//+
Physical Surface("Domain", 17) = {1, 2};
