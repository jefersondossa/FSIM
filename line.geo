// Gmsh project created on Mon Jan 29 15:34:55 2024
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 5.10};
//+
Point(2) = {1, 0, 0, 5.10};
//+
Line(1) = {1, 2};
//+
Transfinite Curve {1} = 21 Using Progression 1;
//+
Physical Point("Left", 2) = {1};
//+
Physical Point("Right", 3) = {2};
//+
Physical Curve("Line", 4) = {1};

