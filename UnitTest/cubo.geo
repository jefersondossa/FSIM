// Gmsh project created on Fri Sep 18 14:26:14 2026
SetFactory("OpenCASCADE");
N = 1;
Box(1) = {0, 0, 0, 1, 1, 1};
Physical Surface("face", 13) = {6, 2, 5, 1, 4, 3};
Physical Volume("volume", 14) = {1};
Transfinite Curve {3, 4, 2, 1, 5, 6, 7, 8, 9, 10, 11, 12} = N + 1 Using Progression 1;
Transfinite Surface {1,2,3,4,5,6};
//Recombine Surface {1,2,3,4,5,6};
Transfinite Volume {1};
