// Gmsh project created on Wed May 06 09:33:18 2026
SetFactory("OpenCASCADE");

// TRELIÇA

Point(1) = {0, 0, 0, 10.0};
Point(2) = {1, 0, 0, 10.0};
Point(3) = {1, 1.5, 0, 10.0};
Point(4) = {0, 1.5, 0, 10.0};
Point(5) = {2, 1.5, 0, 10.0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {1, 3};
Line(6) = {3, 5};
Line(7) = {5, 2};

Physical Point("apoio_fixo", 8) = {4};
Physical Point("apoio_movel", 9) = {1};
Physical Point("forca", 10) = {5};
Physical Curve("barras", 11) = {1, 2, 3, 4, 5, 6, 7};
