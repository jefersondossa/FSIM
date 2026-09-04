// PLACA

// PONTOS
Point(1) = {0, 0, 0, 10.0};
Point(2) = {1, 0, 0, 10.0};
Point(3) = {1, 1, 0, 10.0};
Point(4) = {0, 1, 0, 10.0};

// LINHAS
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// PHYSICAL GROUPS
Physical Curve("linhas", 5) = {3, 4, 2, 1};
