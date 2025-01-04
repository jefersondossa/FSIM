// Cantilever 3D beam.
SetFactory("OpenCASCADE");

Box(1) = {0, 0, 0, 0.1, 1, 5};

// Beam Top-Bottom-Left-Right
Transfinite Curve {9, 10, 11, 12} = 1 Using Progression 1;
Transfinite Curve {4, 2, 8, 6} = 21 Using Progression 1;
// Beam Ends Front-Back
Transfinite Curve {3, 7, 5, 1} = 101 Using Progression 1;

Transfinite Surface {1, 2, 3, 4, 5, 6};

Recombine Surface{1,2,3,4,5,6};
Transfinite Volume{1};

Physical Volume("Volume", 15) = {1};
Physical Surface("Load", 16) = {5};
Physical Surface("Engaste", 17) = {6};
Physical Surface("Livre", 18) = {1, 2, 3, 4};
