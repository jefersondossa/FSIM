// Gmsh project created on Fri Dec 15 15:06:46 2023
SetFactory("OpenCASCADE");
//+
Box(1) = {0, 0, 0, 1, 1, 1};

//+
Transfinite Curve {3, 11, 4, 9, 1, 2, 10, 8, 12, 7, 6, 5} = 2 Using Progression 1;
//+
Transfinite Surface {4};
//+
Transfinite Surface {2};
//+
Transfinite Surface {6};
//+
Transfinite Surface {5};
//+
Transfinite Surface {1};
//+
Transfinite Surface {3};
//+
//Recombine Surface{1,2,3,4,5,6};
Transfinite Volume{1};
//Recombine Volume{1};

//+
//+
Physical Volume("Volume", 15) = {1};
//+
Physical Surface("Base", 16) = {1};
//+
Physical Surface("Load", 17) = {2};
//+
//Physical Surface("Free", 18) = {1, 5, 2, 6};
//+
//Physical Point("FixPoint", 19) = {1};
//+
//Physical Point("PointX", 20) = {2};
//+
Physical Surface("Xfixo", 18) = {3,4};
//+
Physical Surface("zfixo", 19) = {6, 5};
