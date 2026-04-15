SetFactory("OpenCASCADE");

Mesh.RecombineAll = 1;

// Bottom points, z = 0.
Point(1) = {0, 0, 0};
Point(2) = {0.5, 0, 0};
Point(3) = {1, 0, 0};
Point(4) = {0, 0.5, 0};
Point(5) = {0.5, 0.5, 0};
Point(6) = {1, 0.5, 0};
Point(7) = {0, 1, 0};
Point(8) = {0.5, 1, 0};
Point(9) = {1, 1, 0};

// Top points, z = 1.
Point(10) = {0, 0, 1};
Point(11) = {0.5, 0, 1};
Point(12) = {1, 0, 1};
Point(13) = {0, 0.5, 1};
Point(14) = {0.5, 0.5, 1};
Point(15) = {1, 0.5, 1};
Point(16) = {0, 1, 1};
Point(17) = {0.5, 1, 1};
Point(18) = {1, 1, 1};

// Bottom grid curves.
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {4, 5};
Line(4) = {5, 6};
Line(5) = {7, 8};
Line(6) = {8, 9};
Line(7) = {1, 4};
Line(8) = {4, 7};
Line(9) = {2, 5};
Line(10) = {5, 8};
Line(11) = {3, 6};
Line(12) = {6, 9};

// Top grid curves.
Line(13) = {10, 11};
Line(14) = {11, 12};
Line(15) = {13, 14};
Line(16) = {14, 15};
Line(17) = {16, 17};
Line(18) = {17, 18};
Line(19) = {10, 13};
Line(20) = {13, 16};
Line(21) = {11, 14};
Line(22) = {14, 17};
Line(23) = {12, 15};
Line(24) = {15, 18};

// Vertical curves.
Line(25) = {1, 10};
Line(26) = {2, 11};
Line(27) = {3, 12};
Line(28) = {4, 13};
Line(29) = {5, 14};
Line(30) = {6, 15};
Line(31) = {7, 16};
Line(32) = {8, 17};
Line(33) = {9, 18};

// Quadrant 1 surfaces: x=[0,0.5], y=[0,0.5].
Curve Loop(1) = {1, 9, -3, -7};
Plane Surface(1) = {1};
Curve Loop(2) = {13, 21, -15, -19};
Plane Surface(2) = {2};
Curve Loop(3) = {1, 26, -13, -25};
Plane Surface(3) = {3};
Curve Loop(4) = {9, 29, -21, -26};
Plane Surface(4) = {4};
Curve Loop(5) = {3, 29, -15, -28};
Plane Surface(5) = {5};
Curve Loop(6) = {7, 28, -19, -25};
Plane Surface(6) = {6};

// Quadrant 2 surfaces: x=[0.5,1], y=[0,0.5].
Curve Loop(7) = {2, 11, -4, -9};
Plane Surface(7) = {7};
Curve Loop(8) = {14, 23, -16, -21};
Plane Surface(8) = {8};
Curve Loop(9) = {2, 27, -14, -26};
Plane Surface(9) = {9};
Curve Loop(10) = {11, 30, -23, -27};
Plane Surface(10) = {10};
Curve Loop(11) = {4, 30, -16, -29};
Plane Surface(11) = {11};

// Quadrant 3 surfaces: x=[0.5,1], y=[0.5,1].
Curve Loop(12) = {4, 12, -6, -10};
Plane Surface(12) = {12};
Curve Loop(13) = {16, 24, -18, -22};
Plane Surface(13) = {13};
Curve Loop(14) = {12, 33, -24, -30};
Plane Surface(14) = {14};
Curve Loop(15) = {6, 33, -18, -32};
Plane Surface(15) = {15};
Curve Loop(16) = {10, 32, -22, -29};
Plane Surface(16) = {16};

// Quadrant 4 surfaces: x=[0,0.5], y=[0.5,1].
Curve Loop(17) = {3, 10, -5, -8};
Plane Surface(17) = {17};
Curve Loop(18) = {15, 22, -17, -20};
Plane Surface(18) = {18};
Curve Loop(19) = {5, 32, -17, -31};
Plane Surface(19) = {19};
Curve Loop(20) = {8, 31, -20, -28};
Plane Surface(20) = {20};

Surface Loop(1) = {1, 2, 3, 4, 5, 6};
Volume(1) = {1};
Surface Loop(2) = {7, 8, 9, 10, 11, 4};
Volume(2) = {2};
Surface Loop(3) = {12, 13, 11, 14, 15, 16};
Volume(3) = {3};
Surface Loop(4) = {17, 18, 5, 16, 19, 20};
Volume(4) = {4};

// Five elements in each half of x/y, ten elements through z.
Transfinite Curve {1:24} = 8 Using Progression 1;
Transfinite Curve {25:33} = 15 Using Progression 1;
Transfinite Surface {1:20};
Recombine Surface {1:20};
Transfinite Volume {1:4};
Recombine Volume {1:4};

Physical Volume("Internal", 15) = {1, 2, 3, 4};

// Bench support and load points.
Physical Point("Supports", 16) = {1, 3, 7, 9};
Physical Point("Load", 17) = {14};
