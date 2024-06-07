// Gmsh project created on Tue Sep 26 19:10:50 2023

m = 0.08;
Point(1) = {0,0,0,m};
Point(2) = {0,-4,0,m};
Point(3) = {0,4,0,m};

Circle(1) = {2,1,3};
Circle(2) = {3,1,2};

//--------------------

Point(4) = {-1.75,1,0,m};
Point(5) = {-2.75,1,0,m};
Point(6) = {-0.75,1,0,m};

Circle(3) = {5,4,6};
Circle(4) = {6,4,5};

Curve Loop(2) = {3,4};
//Plane Surface(2) = {2};


Point(7) = {1.75,1,0,m};
Point(8) = {2.75,1,0,m};
Point(9) = {0.75,1,0,m};

Circle(5) = {8,7,9};
Circle(6) = {9,7,8};

Curve Loop(3) = {5,6};
//Plane Surface(3) = {3};


//------------------

Point(10) = {-1.75, -1,0,m};
Point(11) = {-1.75, -0.7,0,m};
Point(12) = {-1.75, -1.3,0,m};
Point(13) = {1.75, -1,0,m};
Point(14) = {1.75, -0.7,0,m};
Point(15) = {1.75, -1.3,0,m};

Circle(7) = {11,10,12};
Circle(8) = {15,13,14};

//Plane Surface(5) = {1,2,3,4};

//-------------------------
Point(16) = {0,-0.7,0,m};
Point(17) = {0,-1.3,0,m};
Line(9) = {11,16};
Line(11) = {16,14};
Line(12) = {3,1};
Line(13) = {1,16};
Line(14) = {17,2};
Line(15) = {17,15};
Line(10) = {12,17};

//------------------------

Curve Loop(4) = {-2,12,13,-9,7,10,14};
Curve Loop(5) = {-1,-14,15,8,-11,-13,-12};

Plane Surface(1) = {4,2};
Plane Surface(2) = {5,3};

//+
Physical Curve("Circle", 19) = {4, 3, 5, 6};
//+
Physical Curve("Left", 17) += {2, 12, 13, 9, 7, 10, 14};
//+
Physical Curve("Right", 16) += {1, 12, 13, 11, 8, 15, 14};
//+
Physical Surface("Domain", 20) = {1};
//+
Physical Surface("Domain2", 21) = {2};
//+
Physical Curve("Bottom", 22) = {9, 7, 10, 15, 11, 8};
