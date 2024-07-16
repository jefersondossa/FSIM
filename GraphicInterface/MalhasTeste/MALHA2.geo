// Gmsh project created on Tue Jul 09 21:41:58 2024
SetFactory("OpenCASCADE");
hmaior = 2;
hmenor = 0.2;
//pontos face inferior
Point(1) = {0, 0, 0, hmaior};
Point(2) = {50, 0, 0, hmaior};
Point(3) = {50, 60, 0, hmaior};
Point(4) = {0, 60, 0, hmaior};
Point(5) = {20, 20, 0, hmenor};
Point(6) = {30, 20, 0, hmenor};
Point(7) = {30, 40, 0, hmenor};
Point(8) = {20, 40, 0, hmenor};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
//Line Loop(1) = {1,2,3,4};
//Plane Surface (1) = {1};
//face inferior maior
Line(5) = {5,6 };
Line(6) = {6,7 };
Line(7) = {7,8 };
Line(8) = {8,5 };
//Line Loop(2) = {5,6,7,8};
//Plane Surface (2) = {2};
//face inferior menor
//Plane Surface (1) = {1,2}; remover 1 e 2 e deixar esse (maior menos menor)

//pontos face superior
Point(9) = {0, 0, 20, hmaior};
Point(10) = {50, 0, 20, hmaior};
Point(11) = {50, 60, 20, hmaior};
Point(12) = {0, 60, 20, hmaior};
Point(13) = {20, 20, 10, hmenor};
Point(14) = {30, 20, 10, hmenor};
Point(15) = {30, 40, 10, hmenor};
Point(16) = {20, 40, 10, hmenor};

Line(9) = {9,10};
Line(10) = {10, 11};
Line(11) = {11, 12};
Line(12) = {12, 9};
//Line Loop(3) = {9,10,11,12};
//Plane Surface (3) = {3};
//face superior maior

Line(13) = {13,14 };
Line(14) = {14,15 };
Line(15) = {15,16 };
Line(16) = {16,13 };
//Line Loop(4) = {13,14,15,16};
//Plane Surface (4) = {4};
//face superior menor

Line(17) = {1, 9};
Line(18) = {2, 10};
Line(19) = {3, 11};
Line(20) = {4, 12};
//Line Loop(5) = {17,9,18,1};
//Plane Surface (5) = {5};
//face frente maior

Line(21) = {5, 13};
Line(22) = {6, 14};
Line(23) = {7, 15};
Line(24) = {8, 16};
//Line Loop(6) = {21,13,-22,-5};
//Plane Surface (6) = {6};
//face frente menor





//+
Curve Loop(1) = {13, -22, -5, 21};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {16, -21, -8, 24};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {14, -23, -6, 22};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {15, -24, -7, 23};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {16, 13, 14, 15};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {20, -11, -19, 3};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {19, -10, -18, 2};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {18, -9, -17, 1};
//+
Plane Surface(8) = {8};
//+
Curve Loop(9) = {17, -12, -20, 4};
//+
Plane Surface(9) = {9};
//+
Curve Loop(10) = {9, 10, 11, 12};
//+
Plane Surface(10) = {10};
//+
Curve Loop(11) = {2, 3, 4, 1};
//+
Curve Loop(12) = {6, 7, 8, 5};
//+
Plane Surface(11) = {11, 12};
//+
Show "*";
//+
Show "*";
//+
Surface Loop(1) = {10, 8, 7, 6, 9, 11, 1, 5, 2, 4, 3};
//+
Volume(1) = {1};
//+
Surface Loop(2) = {7, 6, 9, 8, 10, 11, 1, 5, 2, 4, 3};
//+
Physical Surface("edificio", 25) = {4, 2, 5, 3, 1};
//+
Physical Surface("entrada", 26) = {7};
//+
Physical Surface("saida", 27) = {9};
//+
Physical Surface("simetria", 28) = {6, 8};
//+
Physical Surface("topo", 29) = {10};
//+
Physical Surface("base", 30) = {11};
//+
Physical Volume("fluido", 31) = {1};
