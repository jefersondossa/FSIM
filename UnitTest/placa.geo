// PLACA
N = 8;

// PONTOS
Point(1) = {0, 0, 0, 1.0/N};
Point(2) = {1, 0, 0, 1.0/N};
Point(3) = {1, 1, 0, 1.0/N};
Point(4) = {0, 1, 0, 1.0/N};

// LINHAS
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// PHYSICAL GROUPS
Physical Curve("linhas", 5) = {3, 4, 2, 1};
Curve Loop(1) = {4, 1, 2, 3};
Plane Surface(1) = {1};
Physical Surface("superficie", 6) = {1};

// OUTROS
Transfinite Curve {3, 4, 2, 1} = N + 1 Using Progression 1;
Transfinite Surface {1};
Recombine Surface {1};
