//Pontos
Point(1) = {0.707106791, 0, 0, 0};
Point(2) = {0, 0.707106791, 0, 0};
Point(3) = {-0.353355339, 0.353355339, 0, 0};
Point(4) = {0, 0, 0, 0};
Point(5) = {-0.353355339, -0.353355339, 0, 0};
Point(6) = {0, -0.707106781, 0, 0};
Point(7) = {0.353355339, 0.353355339, 0, 0};
Point(8) = {0.353355339, -0.353355339, 0, 0};

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

Line Loop(1) = {6, 7, 10, 5};
Line Loop(2) = {-10, 8, 1, 9};
Line Loop(3) = {-9, 2, 3, 4};

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};

// Physical groups
Physical Surface("Dominio Local", 1) = {1, 2, 3};
Physical Line("Boundary overlapping", 2) = {1, 2, 3, 6, 7, 8};
//Physical Point("Apoio fixo", 3) = {4};

Transfinite Curve {1, 2, 3, 4, 5, 6, 7, 8, 9, 10} = 5; //Using Progression 1;

Transfinite Surface {1, 2, 3};
Recombine Surface {1, 2, 3};

