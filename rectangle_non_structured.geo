SetFactory("OpenCASCADE");

// Geometry
Rectangle(1) = {0, 0, 0, 3, 1};

// --- Structured mesh commands removed ---
// Transfinite Curve {4, 2} = 51 Using Progression 1;
// Transfinite Curve {3, 1} = 151 Using Progression 1;
// Transfinite Surface{1,2,3,4};
Recombine Surface{1,2,3,4};

// --- Mesh settings ---
Mesh.Algorithm = 6; // Frontal-Delaunay (robust for 2D)
Mesh.CharacteristicLengthMin = 0.05;
Mesh.CharacteristicLengthMax = 0.1;

// Physical groups
Physical Surface("Internal", 15) = {1};

Physical Point("RightBottom", 17) = {2};
Physical Point("LeftTop", 23) = {4};
Physical Curve("Left", 16) = {4};
Physical Curve("Top", 18) = {3};
Physical Curve("Bottom", 19) = {1};
