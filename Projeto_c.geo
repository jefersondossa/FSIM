// Gmsh project created on Thu Dec 05 13:10:42 2024
SetFactory("OpenCASCADE");
//+
Box(1) = {0, 0, 0, 18.40, 16.90, 45};
//+
Box(2) = {0, 0, 0, 1.15, 5.40, 45};
//+
Box(3) = {17.25, 0, 0, 1.15, 5.40, 45};
//+
Box(4) = {4.75, 0, 0, 1.45, 3, 45};
//+
Box(5) = {12.20, 0, 0, 1.45, 3, 45};
//+
Box(6) = {7, 4.05, 0, 4.40, 7.45, 45};
//+
Recursive Delete {
  Volume{6}; 
}
//+
Box(6) = {7, 9.45, 0, 4.40, 7.45, 45};
//+
BooleanDifference{ Volume{1}; Delete; }{ Volume{6}; Delete; }
//+
BooleanDifference{ Volume{1}; Delete; }{ Volume{2}; Volume{4}; Volume{5}; Volume{3}; Delete; }
