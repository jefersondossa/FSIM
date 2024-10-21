// Cantilever 3D beam.
SetFactory("OpenCASCADE");

Box(1) = {0, 0, 0, 1, 1, 10};

// Beam Top-Bottom-Left-Right
Transfinite Curve {11, 4, 9, 2, 10, 8, 12, 6} = 6 Using Progression 1;
// Beam Ends Front-Back
Transfinite Curve {3, 7, 5, 1} = 51 Using Progression 1;

Transfinite Surface {1, 2, 3, 4, 5, 6};

Recombine Surface{1,2,3,4,5,6};
Transfinite Volume{1};

Physical Volume("Volume", 15) = {1};
Physical Surface("Load", 16) = {5};
Physical Surface("Engaste", 17) = {6};
Physical Surface("Livre", 18) = {1, 2, 3, 4};
