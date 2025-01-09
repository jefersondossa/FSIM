// Gmsh project created on Tue Oct 22 22:15:47 2024
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {.5*Sin(Pi/4), .5*Sin(Pi/4), 0, 1.0};
//+
Point(3) = {-.5*Sin(Pi/4), .5*Sin(Pi/4), 0, 1.0};
//+
Point(4) = {-.5*Sin(Pi/4), -.5*Sin(Pi/4), 0, 1.0};
//+
Point(5) = {.5*Sin(Pi/4), -.5*Sin(Pi/4), 0, 1.0};

altura = 15;
Point(10) = {0, 0, altura, 1.0};
//+
Point(11) = {.5*Sin(Pi/4), .5*Sin(Pi/4), altura, 1.0};
//+
Point(12) = {-.5*Sin(Pi/4), .5*Sin(Pi/4), altura, 1.0};
//+
Point(13) = {-.5*Sin(Pi/4), -.5*Sin(Pi/4), altura, 1.0};
//+
Point(14) = {.5*Sin(Pi/4), -.5*Sin(Pi/4), altura, 1.0};

//+
Point(6) = {-10, -10, 0, 1.0};
//+
Point(7) = {-10, 10, 0, 1.0};
//+
Point(8) = {30, 10, 0, 1.0};
//+
Point(9) = {30, -10, 0, 1.0};

hvento=20;
Point(15) = {-10, -10, hvento, 1.0};
//+
Point(16) = {-10, 10, hvento, 1.0};
//+
Point(17) = {30, 10, hvento, 1.0};
//+
Point(18) = {30, -10, hvento, 1.0};
//+
Circle(1) = {2, 1, 3};
//+
Circle(2) = {3, 1, 4};
//+
Circle(3) = {4, 1, 5};
//+
Circle(4) = {5, 1, 2};
//+
Line(5) = {4, 6};
//+
Line(6) = {6, 9};
//+
Line(7) = {9, 5};
//+
Line(8) = {7, 3};
//+
Line(9) = {6, 7};
//+
Line(10) = {7, 8};
//+
Line(11) = {8, 2};
//+
Line(12) = {9, 8};
//+
Circle(13) = {11, 10, 12};
//+
Circle(14) = {12, 10, 13};
//+
Circle(15) = {13, 10, 14};
//+
Circle(16) = {14, 10, 11};
//+
Line(17) = {4, 13};
//+
Line(18) = {12, 3};
//+
Line(19) = {5, 14};
//+
Line(20) = {11, 2};
//+
Line(21) = {13, 15};
//+
Line(22) = {12, 16};
//+
Line(23) = {14, 18};
//+
Line(24) = {11, 17};
//+
Line(25) = {17, 18};
//+
Line(26) = {18, 9};
//+
Line(27) = {8, 17};
//+
Line(28) = {17, 16};
//+
Line(29) = {16, 7};
//+
Line(30) = {15, 6};
//+
Line(31) = {15, 18};
//+
Line(32) = {15, 16};
//+
Curve Loop(1) = {20, -4, 19, 16};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {19, -15, -17, 3};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {17, -14, 18, 2};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {1, -18, -13, 20};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {14, 15, 16, 13};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {7, -3, 5, 6};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {5, 9, 8, 2};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {1, -8, 10, 11};
//+
Plane Surface(8) = {8};
//+
Curve Loop(9) = {12, 11, -4, -7};
//+
Plane Surface(9) = {9};
//+
Curve Loop(10) = {5, -30, -21, -17};
//+
Plane Surface(10) = {10};
//+
Curve Loop(11) = {7, 19, 23, 26};
//+
Plane Surface(11) = {11};
//+
Curve Loop(12) = {11, -20, 24, -27};
//+
Plane Surface(12) = {12};
//+
Curve Loop(13) = {22, 29, 8, -18};
//+
Plane Surface(13) = {13};
//+
Curve Loop(14) = {22, -32, -21, -14};
//+
Plane Surface(14) = {14};
//+
Curve Loop(15) = {21, 31, -23, -15};
//+
Plane Surface(15) = {15};
//+
Curve Loop(16) = {16, 24, 25, -23};
//+
Plane Surface(16) = {16};
//+
Curve Loop(17) = {24, 28, -22, -13};
//+
Plane Surface(17) = {17};
//+
Curve Loop(18) = {25, -31, 32, -28};
//+
Plane Surface(18) = {18};
//+
Curve Loop(19) = {29, -9, -30, 32};
//+
Plane Surface(19) = {19};
//+
Curve Loop(20) = {30, 6, -26, -31};
//+
Plane Surface(20) = {20};
//+
Curve Loop(21) = {26, 12, 27, 25};
//+
Plane Surface(21) = {21};
//+
Curve Loop(22) = {10, 27, 28, 29};
//+
Plane Surface(22) = {22};
//+
Surface Loop(1) = {6, 2, 20, 15, 10, 11};
//+
Volume(1) = {1};
//+
Surface Loop(2) = {7, 19, 3, 14, 13, 10};
//+
Volume(2) = {2};
//+
Surface Loop(3) = {13, 22, 8, 4, 12, 17};
//+
Volume(3) = {3};
//+
Surface Loop(4) = {9, 21, 1, 16, 12, 11};
//+
Volume(4) = {4};
//+
Surface Loop(5) = {18, 15, 16, 17, 14, 5};
//+
Volume(5) = {5};
//+
Transfinite Curve {8, 11, 7, 5, 24, 23, 22, 21} = 10 Using Progression 1;
//+
Transfinite Curve {18, 17, 20, 19, 30, 29, 26, 27} = 10 Using Progression 1;
//+
Transfinite Curve {32, 9, 28, 10, 31, 6, 12, 25, 15, 16, 13, 14, 1, 2, 3, 4} = 10 Using Progression 1;
//+
Transfinite Surface {1};
Transfinite Surface {2};
Transfinite Surface {3};
Transfinite Surface {4};
Transfinite Surface {8};
//+
Transfinite Surface {9};
//+
Transfinite Surface {7};
//+
Transfinite Surface {6};
Transfinite Surface {5};
//+
Transfinite Surface {10};
//+
Transfinite Surface {13};
//+
Transfinite Surface {12};
//+
Transfinite Surface {11};
//+
Transfinite Surface {18};
//+
Transfinite Surface {16};
//+
Transfinite Surface {15};
//+
Transfinite Surface {14};
//+
Transfinite Surface {17};
//+
Transfinite Surface {20};
//+
Transfinite Surface {19};
//+
Transfinite Surface {22};
//+
Transfinite Surface {21};
//+
Transfinite Volume{3};
//+
Transfinite Volume{4};
//+
Transfinite Volume{5};
//+
Transfinite Volume{1};
//+
Transfinite Volume{2};
//+

//+
Physical Surface("Inflow", 33) = {19};
//+
Physical Surface("Outflow", 34) = {21};
//+
Physical Surface("Floor", 35) = {9, 6, 7, 8};
//+
Physical Surface("Cylinder", 36) = {3, 4, 1, 2, 5};
//+
Physical Surface("NormalZ", 37) = {18};
//+
Physical Surface("NormalY", 38) = {22, 20};
//+
Physical Volume("Fluid", 39) = {1, 4, 3, 2, 5};
