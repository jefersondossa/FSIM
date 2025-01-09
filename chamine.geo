// Gmsh project created on Wed Dec 18 10:50:43 2024
SetFactory("OpenCASCADE");
//+
Circle(1) = {0, 0, 0, 0.218, 0, 2*Pi};
Extrude {0, 0, 20} {
  Curve{1}; Layers {20};
}
//+
Circle(4) = {0, 0, 20, 0.190, 0, 2*Pi};

Curve Loop(2) = {3};
//+
Curve Loop(3) = {4};
//+
Plane Surface(2) = {2, 3};
Extrude {0, 0, 93} {
  Curve{4};  Layers {93};
}
//+
Curve Loop(5) = {6};
//+
Plane Surface(4) = {5};
//+
Curve Loop(6) = {1};
//+
Plane Surface(5) = {6};
