//Pontos
Point(1) = {0, 0, 0, 1.0};
Point(2) = {1, 0, 0, 1.0};
Point(3) = {1, 1, 0, 1.0};
Point(4) = {0, 1, 0, 1.0};

//Linhas
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(1) = {1, 2, 3, 4};

Plane Surface(1) = {1};

//Physical groups
Physical Surface("Dominio", 1) = {1};
Physical Line("Desloc. prescrito", 2) = {3};
Physical Line("Apoio móvel vertical", 3) = {4};
Physical Line("Apoio horizontal", 4) = {1};

Transfinite Curve {1, 2, 3, 4} = 3;

Transfinite Surface {1};
Recombine Surface {1};