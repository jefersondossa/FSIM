// Gmsh project created on Mon Jun 08 09:57:51 2026
SetFactory("OpenCASCADE");

//PONTOS
Point(1) = {1, 0, 0, 10.0};
Point(2) = {4, 0, 0, 10.0};
Point(3) = {7, 0, 0, 10.0};

//RETAS
Line(1) = {1, 2};
Line(2) = {2, 3};

//GRUPOS
Physical Point("apoio_fixo_esquerda", 3) = {1};
Physical Point("apoio_movel_centro", 4) = {2};
Physical Point("apoio_movel_direita", 5) = {3};
Physical Curve("barras", 6) = {1, 2};
