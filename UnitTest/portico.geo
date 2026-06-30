// Gmsh project created on Mon Jun 08 09:46:53 2026
SetFactory("OpenCASCADE");

//PONTOS
Point(1) = {1, 0, 0, 10.0};
Point(2) = {1, 4, 0, 10.0};
Point(3) = {1, 7, 0, 10.0};
Point(4) = {1, 10, 0, 10.0};
Point(5) = {6, 0, 0, 10.0};
Point(6) = {6, 7, 0, 10.0};
Point(8) = {6, 4, 0, 10.0};
Point(9) = {6, 10, 0, 10.0};

//RETAS
Line(1) = {2, 8};
Line(2) = {3, 6};
Line(3) = {4, 9};
Line(4) = {1, 2};
Line(5) = {2, 3};
Line(6) = {3, 4};
Line(7) = {5, 8};
Line(8) = {8, 6};
Line(9) = {6, 9};

//GRUPOS
Physical Point("engaste_esquerdo", 10) = {1};
Physical Point("engaste_direito", 11) = {5};
Physical Point("forca_40knewton", 12) = {9};
Physical Curve("carga_distribuida_10kn", 13) = {1, 2, 3};
Physical Curve("barras", 14) = {4, 5, 6, 9, 8, 7};
