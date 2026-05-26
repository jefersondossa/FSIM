//Pontos:
Point(1) = {1.414213562, 0, 0, 0};
Point(2) = {0, 1.414213562, 0, 0};
Point(3) = {-0.707106781, 0.707106781, 0, 0};
Point(4) = {0, 0, 0, 0};
Point(5) = {-0.707106781, -0.707106781, 0, 0};
Point(6) = {0, -1.414213562, 0, 0};
Point(7) = {0.707106781, 0.707106781, 0, 0};
Point(8) = {0.707106781, -0.707106781, 0, 0};

//Linhas:
Line(1) = {1, 7};
Line(2) = {7, 2};
Line(3) = {2, 3};
Line(4) = {3, 4};
Line(5) = {4, 5};
Line(6) = {5, 6};
Line(7) = {6, 8};
Line(8) = {8, 1};
Line(9) = {7, 4};
Line(10) = {8, 4};

Line Loop(1) = {8, 1, 9, -10};
Line Loop(2) = {-9, 2, 3, 4};
Line Loop(3) = {7, 10, 5, 6};
Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};

// Physical groups
Physical Surface("Dominio", 1) = {1, 2, 3};
Physical Point("Apoio movel horizontal", 2) = {2, 6};
Physical Point("Apoio movel vertical", 3) = {1};
Physical Line("C. C. AB", 4) = {1, 2}; 
Physical Line("C. C. BC", 5) = {3}; 
Physical Line("C. C. EF", 6) = {6}; 
Physical Line("C. C. FA", 7) = {7, 8}; 

Transfinite Curve {1, 2, 3, 4, 5, 6, 7, 8, 9, 10} = 3; //Using Progression 1;

Transfinite Surface {1, 2, 3};
Recombine Surface {1, 2, 3};
