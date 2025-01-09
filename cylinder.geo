// Gmsh project created on Tue Jan 30 18:59:12 2024
SetFactory("OpenCASCADE");
//+
Rectangle(1) = {0, 0, 0, 30, 10, 0};
//+
//Circle(5) = {2.5, 2.5, 0.0, 0.25, 0, 2*Pi};

//+
//Curve Loop(2) = {3, 4, 1, 2};
//+
//Curve Loop(3) = {5};
//+
//Plane Surface(2) = {2, 3};
//+
Physical Curve("Inflow", 6) = {4};
//+
Physical Curve("Outflow", 7) = {2};
//+
Physical Curve("Symmetry", 8) = {3, 1};
//+
//Physical Curve("Wall", 9) = {5};
//+
Physical Surface("Fluid", 10) = {1};
