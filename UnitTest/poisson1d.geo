// Gmsh project created on Mon Jun 08 09:57:51 2026
SetFactory("OpenCASCADE");

//PONTOS
Point(1) = {0, 0, 0, 0.5};
Point(2) = {1, 0, 0, 0.5};

//RETAS
Line(1) = {1, 2};

//GRUPOS
Physical Point("apoio_fixo_esquerda", 3) = {1};
Physical Point("apoio_movel_direita", 4) = {2};
Physical Curve("barras", 6) = {1};
