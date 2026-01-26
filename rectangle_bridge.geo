// =========================================
// Manual geometry (built-in kernel)
// Rectangle: width = 2, height = 1
// Bottom split into 2 lines using midpoint
// =========================================

// Corner points
Point(1) = {0, 0, 0};
Point(2) = {2, 0, 0};
Point(3) = {2, 1, 0};
Point(4) = {0, 1, 0};

// Mid-bottom point
Point(10) = {1, 0, 0};

// Bottom edge split into two lines
Line(101) = {1, 10};   // left half
Line(102) = {10, 2};   // right half

// Other edges
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Surface
Line Loop(1) = {101, 102, 2, 3, 4};
Plane Surface(1) = {1};

// Mesh control
Transfinite Curve {101, 102, 2, 4} = 41 Using Progression 1;
Transfinite Curve {101, 102, 2, 4} = 41 Using Progression 1;
Transfinite Curve {3} = 81 Using Progression 1;
Transfinite Surface{1} = {1, 2, 3, 4};
Recombine Surface{1};

// Physical groups
Physical Surface("Internal", 15) = {1};

Physical Point("Bottom_Left",  24) = {1};
Physical Point("Bottom_Right", 17) = {2};

//Physical Curve("Right", 203) = {2};
//Physical Curve("Top",   18) = {3};
//Physical Curve("Left",  16) = {4};

Physical Point("BottomMiddle", 23) = {10};
