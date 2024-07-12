// Gmsh project created on Sat Sep 30 14:10:51 2023
//Geograde

m = 1;
n = 0.5;

Point(10) = {0,0,0,m};
Point(2) = {15,0,0,m};
Point(3) = {15,15,0,m};
Point(4) = {0,15,0,m};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Curve Loop(1) = {1,2,3,4};

Point(5) = {1,12,0,n};
Point(6) = {14,12,0,n};

Line(5) = {5,6};

Line Loop(2) = {5,-5};

Plane Surface(1) = {1,2};

Physical Surface("Domain", 6) = {1};
//+
Physical Curve("Right", 7) = {5};
//+
Physical Curve("Left", 8) = {4, 2};
//+
Physical Curve("Top", 9) = {3};
//+
Physical Curve("Bottom", 10) = {1};
