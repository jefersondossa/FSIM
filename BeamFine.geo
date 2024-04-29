// Gmsh project created on Thu Apr 25 09:32:35 2024
//Rectangle 2 X 0.6
//Crack on points 5,6 e 7

m = 0.5;
gluezone = 0.1;
altura = 0.6;
alturatrinca = altura/10;
postrinca = 3;
largura = 1;

Point(1) = {postrinca-largura,0,0,m};
Point(2) = {postrinca-largura,altura,0,m};
Point(3) = {postrinca+largura,altura,0,m};
Point(4) = {postrinca+largura,0,0,m};
Point(5) = {postrinca,0,0,m};
Point(6) = {postrinca,alturatrinca,0,m};
Point(7) = {postrinca,0,0,m};

Point(8) = {postrinca-largura,alturatrinca,0,m};
Point(9) = {postrinca+largura,alturatrinca,0,m};
Point(10) = {postrinca,altura,0,m};


Line(1) = {1,8};
Line(2) = {8,2};
Line(3) = {2,10};
Line(12) = {10,3};
Line(4) = {3,9};
Line(5) = {4,5};
Line(6) = {5,6};
Line(7) = {6,7};
Line(8) = {7,1};
Line(9) = {8,6};
Line(10) = {6,9};
Line(11) = {9,4};
Line(13) = {6,10};

//+
Curve Loop(1) = {2, 3, -13, -9};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {12, 4, -10, 13};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {8, 1, 9, 7};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {10, 11, 5, 6};
//+
Plane Surface(4) = {4};
//+
//+
Transfinite Curve {2, 13, 4} = 10 Using Progression 1;
//+
Transfinite Curve {3, 12, 10, 9, 8, 5} = 20 Using Progression 1;
//+
Transfinite Curve {1, 7, 6, 11} = 5 Using Progression 1;

//+
Transfinite Surface {1};
//+
Transfinite Surface {2};
//+
Transfinite Surface {3};
//+
Transfinite Surface {4};
Transfinite Surface {5};
Transfinite Surface {6};
Transfinite Surface {7};
Transfinite Surface {8};
//+
//Recombine Surface {3, 4, 2, 1, 5, 6, 7, 8};
//+
Physical Curve("Glue", 19) = {2, 1, 4, 11};
//+
Physical Curve("Load", 15) = {12, 3};
//+
Physical Surface("Domain", 16) = {2, 4, 3, 1};


