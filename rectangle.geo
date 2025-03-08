SetFactory("OpenCASCADE");

Rectangle(1) = {0, 0, 0, 2, 1};

// Beam Top-Bottom-Left-Right
// LEFT-RIGHT
Transfinite Curve {4, 2} = 21 Using Progression 1;
// TOP-BOTTOM
Transfinite Curve {3, 1} = 41 Using Progression 1;
Transfinite Surface{1,2,3,4};
Recombine Surface{1,2,3,4};

Physical Surface("Internal", 15) = {1};

Physical Curve("Bottom", 19) = {1};
Physical Curve("Right", 17) = {2};
Physical Curve("Top", 18) = {3};
Physical Curve("Left", 16) = {4};
