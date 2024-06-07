Point(1) = {0.000000, 0.000000, 0.000000, 1.000000}; 
//
Point(2) = {1.000000, 0.000000, 0.000000, 1.000000}; 
//
Point(3) = {1.000000, 0.300000, 0.000000, 1.000000}; 
//
Point(4) = {0.000000, 0.300000, 0.000000, 1.000000}; 
//
Line(1) = {1, 2}; 
//
Line(2) = {2, 3}; 
//
Line(3) = {3, 4}; 
//
Line(4) = {4, 1};
//
Line Loop(1) = {1, 2, 3, 4};
//
Plane Surface(1) = {1};
//
Transfinite Line {1} = 37 Using Progression 1;
//
Transfinite Line {2} = 11 Using Progression 1;
//
Transfinite Line {3} = 37 Using Progression 1;
//
Transfinite Line {4} = 11 Using Progression 1;
//+
Physical Curve("Engaste", 11) = {4};
//+
Physical Curve("Forca", 12) = {3};
//+
Physical Surface("Dominio", 13) = {1};
