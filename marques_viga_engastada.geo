@ -0,0 +1,43 @@
// Gmsh project created on Thu Sep 28 15:16:20 2023
SetFactory("OpenCASCADE");
Point(1) = {0, 0, 0, 1.0};
//+
Point(3) = {1., 0, 0, 1.0};
//+
Point(4) = {1., 1, 0, 1.0};
//+
Point(6) = {0, 1, 0, 1.0};


//+
Line(1) = {6, 1};
//+
Line(3) = {4, 3};
//+
Line(4) = {3, 1};
//+
//+
Line(6) = {6, 4};

//+
Curve Loop(1) = {1, -4, -3, -6};
//+
Plane Surface(1) = {1};
//+
Transfinite Curve {1,3} = 2 Using Progression 1;
//+
Transfinite Curve {4,6} = 2 Using Progression 1;
//+
Transfinite Surface {1};

//+
Physical Curve("Left", 9) = {4};

Physical Curve("Right", 12) = {6};
//+
Physical Surface("Domain", 10) = {1};
//+
//Physical Curve("Bottom", 11) = {1};

//+
Physical Curve("Laterais", 13) = {1, 3};

