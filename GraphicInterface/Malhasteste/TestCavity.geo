m = 0.1;

//Building
//Inside Points

Point(1) = {0,0,0,m};
Point(2) = {3,0,0,m};
Point(3) = {3,3,0,m};
Point(4) = {0,3,0,m};

Point(5) = {0,0,1,m};
Point(6) = {3,0,1,m};
Point(7) = {3,3,1,m};
Point(8) = {0,3,1,m};
//+
Line(1) = {5, 8};
//+
Line(2) = {8, 4};
//+
Line(3) = {4, 1};
//+
Line(4) = {1, 5};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 2};
//+
Line(7) = {2, 1};
//+
Line(8) = {6, 7};
//+
Line(9) = {7, 3};
//+
Line(10) = {3, 2};
//+
Line(11) = {7, 8};
//+
Line(12) = {4, 3};

//+
Curve Loop(1) = {3, 4, 1, 2};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {12, -9, 11, 2};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {10, -6, 8, 9};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {7, 4, 5, 6};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {3, -7, -10, -12};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {1, -11, -8, -5};
//+
Plane Surface(6) = {6};

//+
Surface Loop(1) = {1, 5, 4, 6, 2, 3};
//+
Volume(1) = {1};

//+
//Physical Curve("Curve inlet", 13) = {3, 2, 1, 4};
//+
//Physical Curve("Curve outlet", 14) = {10, 9, 8, 6};
//+
//Physical Curve("Curve top", 15) = {2, 12, 9, 11};
//+
//Physical Curve("Curve botton", 16) = {4, 7, 6, 5};
//+
//Physical Curve("Curve front", 17) = {1, 11, 8, 5};
//+
//Physical Curve("Curve back", 18) = {7, 3, 12, 10};


//+
Physical Surface("movingWall", 19) = {2};
//+
Physical Surface("fixedWalls", 20) = {1, 4, 3};
//+
Physical Surface("frontAndBack", 21) = {6, 5};


//+
Physical Volume("Volume", 22) = {1};
