// Gmsh project created on Mon Oct 30 10:12:05 2023
SetFactory("OpenCASCADE");
//+
Rectangle(1) = {0, 0, 0, 1, 1, 0};
//+
Physical Curve("Left", 5) = {4};
//+
Physical Curve("Right", 6) = {2};
//+
Physical Curve("TB", 7) = {3, 1};
//+
Physical Surface("Domain", 8) = {1};
//+
Transfinite Curve {3, 1} = 3 Using Progression 1;
Transfinite Curve {4,2} = 2 Using Progression 1;
Transfinite Surface {1};
