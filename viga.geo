//Pontos
Point(1) = {0, 0, 0, 1.0};
Point(2) = {100, 0, 0, 1.0};
Point(3) = {100, 10, 0, 1.0};
Point(4) = {0, 10, 0, 1.0};

//Linhas
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(1) = {1, 2, 3, 4};

Plane Surface(1) = {1};

// Physical groups
Physical Surface("Dominio Global", 1) = {1};
Physical Point("Apoio fixo", 2) = {1};
Physical Point("Apoio móvel horizontal", 3) = {4};
Physical Line("Momento", 4) = {2};

Transfinite Curve {1, 2, 3, 4} = 2;

Transfinite Surface {1};
Recombine Surface {1};


