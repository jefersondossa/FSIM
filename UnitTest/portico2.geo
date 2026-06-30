// Gmsh project created on Mon Jun 29 12:07:16 2026
SetFactory("OpenCASCADE");

//PONTOS
Point(1) = {0, 0, 0, 10.0};
Point(2) = {0, 10, 0, 10.0};

//LINHAS
Line(1) = {1, 2};

//GRUPOS
Physical Point("engaste", 1) = {1};
Physical Curve("coluna", 2) = {1};
Physical Point("carga_concentrada", 3) = {2};
