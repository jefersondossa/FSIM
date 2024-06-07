// Gmsh project created on Thu Sep 21 15:43:17 2023

m = 0.8;
Point(1) = {0,0,0,m};
Point(2) = {-4,0,0,m};
Point(3) = {4,0,0,m};

Circle(1) = {2,1,3};
Circle(2) = {3,1,2};

Curve Loop(1) = {1,2};
//Plane Surface(1) = {1};

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
Line(9) = {14,11};
Line (10) = {12,15};

Curve Loop(4) = {7,8,9,10};
//Plane Surface(4) = {4};

Plane Surface(5) = {1,2,3,4};

//------------------------

//+
Physical Curve("Left", 11) = {4, 3, 5, 6};
//+
Physical Curve("Right", 12) = {9, 7, 10, 8};
//+
Physical Surface("Domain", 13) = {5};
//+
Physical Curve("Top", 14) = {2, 1};
