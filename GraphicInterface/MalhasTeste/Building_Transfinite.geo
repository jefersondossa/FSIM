h = 1;
//+
Point(1) = {0, -50, 0, h};
//+
Point(2) = {50, -50, 0, h};
//+
Point(3) = {50, 140, 0, h};
//+
Point(4) = {0, 140, 0, h};
//+
Point(5) = {20, 20, 0, h};
//+
Point(6) = {30, 20, 0, h};
//+
Point(7) = {30, 50, 0, h};
//+
Point(8) = {20, 50, 0, h};

Point(9) = {20, 20, 30, h};
//+
Point(10) = {30, 20, 30, h};
//+
Point(11) = {30, 50, 30, h};
//+
Point(12) = {20, 50, 30, h};

Point(13) = {0, -50, 50, h};
//+
Point(14) = {50, -50, 50, h};
//+
Point(15) = {50, 140, 50, h};
//+
Point(16) = {0, 140, 50, h};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Line(5) = {5, 8};
//+
Line(6) = {8, 7};
//+
Line(7) = {7, 6};
//+
Line(8) = {6, 5};
//+
Line(9) = {9, 10};
//+
Line(10) = {10, 11};
//+
Line(11) = {11, 12};
//+
Line(12) = {12, 9};
//+
Line(13) = {10, 6};
//+
Line(14) = {9, 5};
//+
Line(15) = {11, 7};
//+
Line(16) = {8, 12};
//+
Line(17) = {13, 14};
//+
Line(18) = {14, 15};
//+
Line(19) = {15, 16};
//+
Line(20) = {16, 13};
//+
Line(21) = {13, 1};
//+
Line(22) = {16, 4};
//+
Line(23) = {15, 3};
//+
Line(24) = {11, 15};
//+
Line(25) = {7, 3};
//+
Line(26) = {6, 2};
//+
Line(27) = {5, 1};
//+
Line(28) = {9, 13};
//+
Line(29) = {12, 16};
//+
Line(30) = {8, 4};
//+
Line(31) = {10, 14};
//+
Line(32) = {14, 2};
//+
Curve Loop(1) = {6, 25, 3, -30};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {25, -2, -26, -7};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {1, -26, 8, 27};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {27, -4, -30, -5};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {5, 16, 12, 14};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {16, -11, 15, -6};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {12, 9, 10, 11};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {8, -14, 9, 13};
//+
Plane Surface(8) = {8};
//+
Curve Loop(9) = {7, -13, 10, 15};
//+
Plane Surface(9) = {9};
//+
Curve Loop(10) = {15, 25, -23, -24};
//+
Plane Surface(10) = {10};
//+
Curve Loop(11) = {30, -22, -29, -16};
//+
Plane Surface(11) = {11};
//+
Curve Loop(12) = {14, 27, -21, -28};
//+
Plane Surface(12) = {12};
//+
Curve Loop(13) = {26, -32, -31, 13};
//+
Plane Surface(13) = {13};
//+
Curve Loop(14) = {24, 19, -29, -11};
//+
Plane Surface(14) = {14};
//+
Curve Loop(15) = {29, 20, -28, -12};
//+
Plane Surface(15) = {15};
//+
Curve Loop(16) = {28, 17, -31, -9};
//+
Plane Surface(16) = {16};
//+
Curve Loop(17) = {31, 18, -24, -10};
//+
Plane Surface(17) = {17};
//+
Curve Loop(18) = {23, 3, -22, -19};
//+
Plane Surface(18) = {18};
//+
Curve Loop(19) = {22, 4, -21, -20};
//+
Plane Surface(19) = {19};
//+
Curve Loop(20) = {21, 1, -32, -17};
//+
Plane Surface(20) = {20};
//+
Curve Loop(21) = {32, 2, -23, -18};
//+
Plane Surface(21) = {21};
//+
Curve Loop(22) = {20, 17, 18, 19};
//+
Plane Surface(22) = {22};
//+
Surface Loop(1) = {21, 2, 9, 13, 10, 17};
//+
Volume(1) = {1};
//+
Surface Loop(2) = {22, 15, 16, 14, 17, 7};
//+
Volume(2) = {2};
//+
Surface Loop(3) = {3, 20, 8, 13, 16, 12};
//+
Volume(3) = {3};
//+
Surface Loop(4) = {19, 4, 5, 12, 15, 11};
//+
Volume(4) = {4};
//+
Surface Loop(5) = {14, 6, 1, 18, 10, 11};
//+
Volume(5) = {5};

Physical Surface("Building", 33) = {5,6,8,7,9};
Physical Surface("Inlet", 34) = {20};
Physical Surface("Outlet", 35) = {18};
Physical Surface("frontAndBack", 36) = {21, 19};
Physical Surface("Up", 37) = {22};
Physical Surface("Down", 38) = {1, 2, 3, 4};
Physical Volume("Fluid", 39) = {1, 2, 3, 4, 5};

//+
Transfinite Curve {4, 20, 3, 19, 18, 18, 17, 21, 1, 32, 22, 23, 2, 5, 7, 6, 14, 8, 9, 13, 15, 16, 10, 12, 11} = 10 Using Progression 1;
//+
Transfinite Curve {26, 31, 28, 27, 30, 25, 24, 29} = 15 Using Progression 2;
//+
Transfinite Surface {21};
//+
Transfinite Surface {18};
//+
Transfinite Surface {20};
//+
Transfinite Surface {22};
//+
Transfinite Surface {17};
//+
Transfinite Surface {10};
//+
Transfinite Surface {14};
//+
Transfinite Surface {1};
//+
Transfinite Surface {2};
//+
Transfinite Surface {4};
//+
Transfinite Surface {15};
//+
Transfinite Surface {12};
//+
Transfinite Surface {6};
//+
Transfinite Surface {3};
//+
Transfinite Surface {16};
//+
Transfinite Surface {19};
//+
Transfinite Surface {11};
//+
Transfinite Surface {7};
//+
Transfinite Surface {9};
//+
Transfinite Surface {8};
//+
Transfinite Surface {5};
//+
Transfinite Surface {13};
//+
Recombine Surface {15, 22, 19, 14, 18, 11, 10, 5, 9, 4, 21, 12, 13, 8, 7, 17, 16, 20, 3, 2, 1, 6};
//+
Transfinite Volume{5};
//+
Transfinite Volume{1};
//+
Transfinite Volume{2};
//+
Transfinite Volume{4};
//+
Transfinite Volume{3};
