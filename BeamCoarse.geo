//Gmsh project created on Thu Apr 25 09:06:35 2024
//Rectangle 6 X 0.6 (Lenght Suport = 0.75)

m = 0.3;

Point(1) = {0,0,0,m};
Point(2) = {0,0.6,0,m};
Point(3) = {6,0.6,0,m};
Point(4) = {6,0,0,m};
Point(5) = {5.25,0,0,m};
Point(6) = {0.75,0,0,m};
Point(7) = {0.75,0.6,0,m};
Point(8) = {5.25,0.6,0,m};

Line(1) = {1,2};
Line(2) = {2,7};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,1};
Line(7) = {6,7};
Line(8) = {5,8};
Line(9) = {7,8};
Line(10) = {8,3};

//+
Curve Loop(1) = {1, 2, -7, 6};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {9, -8, 5, 7};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {10, 3, 4, 8};
//+
Plane Surface(3) = {3};
//+
Transfinite Curve {1, 7, 8, 3} = 10 Using Progression 1;
//+
Transfinite Curve {9,5} = 20 Using Progression 1;
//+
Transfinite Curve {2, 10, 4, 6} = 4.45 Using Progression 1;
//+
Transfinite Surface {1};
//+
Transfinite Surface {2};
//+
Transfinite Surface {3};
//+
Recombine Surface {1, 2, 3};
//+
Physical Curve("Load", 11) = {2, 9, 10};
//+
Physical Curve("SupportLeft", 12) = {6};
//+
Physical Curve("SupportRight", 13) = {4};
//+
Physical Surface("DomainCoarse", 14) = {1, 2, 3};
