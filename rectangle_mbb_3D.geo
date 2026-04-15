SetFactory("OpenCASCADE");

// Keep the original rectangle point IDs and coordinates on z = 0.
Point(1) = {0, 0, 0};
Point(2) = {1, 0, 0};
Point(3) = {1, 1, 0};
Point(4) = {0, 1, 0};

Point(5) = {0, 0, 1};
Point(6) = {1, 0, 1};
Point(7) = {1, 1, 1};
Point(8) = {0, 1, 1};

// Original rectangle curves on z = 0:
// 1 = Bottom, 2 = Right, 3 = Top, 4 = Left.
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};

Line(9) = {1, 5};
Line(10) = {2, 6};
Line(11) = {3, 7};
Line(12) = {4, 8};

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Curve Loop(2) = {5, 6, 7, 8};
Plane Surface(2) = {2};

Curve Loop(3) = {1, 10, -5, -9};
Plane Surface(3) = {3};

Curve Loop(4) = {2, 11, -6, -10};
Plane Surface(4) = {4};

Curve Loop(5) = {3, 12, -7, -11};
Plane Surface(5) = {5};

Curve Loop(6) = {4, 9, -8, -12};
Plane Surface(6) = {6};

Surface Loop(1) = {1, 2, 3, 4, 5, 6};
Volume(1) = {1};

// Beam Top-Bottom-Left-Right
// LEFT-RIGHT
Transfinite Curve {4, 2, 8, 6} = 11 Using Progression 1;
// TOP-BOTTOM
Transfinite Curve {3, 1, 7, 5} = 11 Using Progression 1;
// DEPTH
Transfinite Curve {9, 10, 11, 12} = 11 Using Progression 1;

Transfinite Surface {1, 2, 3, 4, 5, 6};
Recombine Surface {1, 2, 3, 4, 5, 6};
Transfinite Volume {1};
Recombine Volume {1};

Physical Volume("Internal", 15) = {1};

// MBB
Physical Surface("Left", 16) = {6};
Physical Curve("RightBottom", 17) = {10};
Physical Surface("Top", 18) = {5};
Physical Surface("Bottom", 19) = {3};
Physical Curve("LeftTop", 23) = {12};
