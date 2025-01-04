// Cantilever 3D beam.
SetFactory("OpenCASCADE");

// Box(1) = {0, 0, 0, 0.1, 1, 5};
Rectangle(1) = {0, 0, 0, 2, 1};

// Beam Top-Bottom-Left-Right
// LEFT-RIGHT
Transfinite Curve {4, 2} = 21 Using Progression 1;
// TOP-BOTTOM
Transfinite Curve {3, 1} = 41 Using Progression 1;

//Recombine Surface{1,2,3,4};

Physical Surface("Internal", 15) = {1};
Physical Curve("Load", 16) = {2};
Physical Curve("Engaste", 17) = {4};
Physical Curve("Livre", 18) = {1, 3};
