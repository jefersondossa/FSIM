//Pontos:
Point(1) = {1.414213562, 0, 0, 0};
Point(2) = {1.060660171, 0.35355339, 0, 0};
Point(3) = {0.707106781, 0.707106781, 0, 0};
Point(4) = {0.35355339, 1.060660171, 0, 0};
Point(5) = {0, 1.414213562, 0, 0};
Point(6) = {-0.35355339, 1.060660171, 0, 0};
Point(7) = {-0.707106781, 0.707106781, 0, 0};
Point(8) = {-0.35355339, 0.35355339, 0, 0};
Point(9) = {0, 0, 0, 0};
Point(10) = {-0.35355339, -0.35355339, 0, 0};
Point(11) = {-0.707106781, -0.707106781, 0, 0};
Point(12) = {-0.35355339, -1.060660171, 0, 0};
Point(13) = {0, -1.414213562, 0, 0};
Point(14) = {0.35355339, -1.060660171, 0, 0};
Point(15) = {0.707106781, -0.707106781, 0, 0};
Point(16) = {1.060660171, -0.35355339, 0, 0};
Point(17) = {0.707106781, 0, 0, 0};
Point(18) = {0, 0.707106781, 0, 0};
Point(19) = {0, -0.707106781, 0, 0};
Point(20) = {0.35355339, 0.35355339, 0, 0};
Point(21) = {0.35355339, -0.35355339, 0, 0};

//Linhas:
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 9};
Line(9) = {9, 10};
Line(10) = {10, 11};
Line(11) = {11, 12};
Line(12) = {12, 13};
Line(13) = {13, 14};
Line(14) = {14, 15};
Line(15) = {15, 16};
Line(16) = {16, 1};

//Boundary overlapping
Line(17) = {10, 19};
Line(18) = {19, 21};
Line(19) = {21, 17};
Line(20) = {17, 20};
Line(21) = {20, 18};
Line(22) = {18, 8};

Line(23) = {12, 19};
Line(24) = {14, 19};
Line(25) = {15, 21};
Line(26) = {16, 17};
Line(27) = {2, 17};
Line(28) = {3, 20};
Line(29) = {4, 18};
Line(30) = {6, 18};
Line(31) = {21, 9};
Line(32) = {20, 9};

Line Loop(1) = {10, 11, 23, -17};
Line Loop(2) = {12, 13, 24, -23};
Line Loop(3) = {-24, 14, 25, -18};
Line Loop(4) = {-25, 15, 26, -19};
Line Loop(5) = {-26, 16, 1, 27};
Line Loop(6) = {-27, 2, 28, -20};
Line Loop(7) = {-28, 3, 29, -21};
Line Loop(8) = {-29, 4, 5, 30};
Line Loop(9) = {7, -22, -30, 6};
Line Loop(10) = {17, 18, 31, 9};
Line Loop(11) = {-31, 19, 20, 32};
Line Loop(12) = {8, -32, 21, 22};

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};
Plane Surface(6) = {6};
Plane Surface(7) = {7};
Plane Surface(8) = {8};
Plane Surface(9) = {9};
Plane Surface(10) = {10};
Plane Surface(11) = {11};
Plane Surface(12) = {12};

// Physical groups
Physical Surface("Dominio Global", 1) = {1, 2, 3, 4, 5, 6, 7, 8, 9};
Physical Surface("Dominio Local", 2) = {10, 11, 12};
Physical Point("Apoio horizontal", 3) = {5, 13};
//Physical Point("Apoio fixo", 3) = {9};
Physical Point("Apoio vertical", 4) = {1};
Physical Line("C. C. AB", 5) = {1, 2, 3, 4}; 
Physical Line("C. C. BC", 6) = {5, 6}; 
Physical Line("C. C. EF", 7) = {11, 12}; 
Physical Line("C. C. FA", 8) = {13, 14, 15, 16}; 

Transfinite Curve {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32} = 3; //Using Progression 1;

Transfinite Surface {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
Recombine Surface {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
