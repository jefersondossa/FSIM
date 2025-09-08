SetFactory("OpenCASCADE");

Rectangle(1) = {0, 0, 0, 3, 1};

// Beam Top-Bottom-Left-Right
// LEFT-RIGHT
Transfinite Curve {4, 2} = 101 Using Progression 1;
// TOP-BOTTOM
Transfinite Curve {3, 1} = 301 Using Progression 1;
Transfinite Surface{1,2,3,4};
Recombine Surface{1,2,3,4};

Physical Surface("Internal", 15) = {1};

Physical Point("RightBottom", 17) = {2};
Physical Point("LeftTop", 23) = {4};
//Physical Point("RightTop", 21) = {3};
//Physical Point("LeftPoints", 22) = {4,1};
//Physical Curve("Right", 17) = {2};
Physical Curve("Bottom", 19) = {1};
Physical Curve("Top", 18) = {3};
Physical Curve("Left", 16) = {4};
