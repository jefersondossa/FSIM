//Points
Point(1) = {0, 0, 0, 1.0};
Point(2) = {4, 0, 0, 1.0};
Point(3) = {4, 2, 0, 1.0};
Point(4) = {0, 2, 0, 1.0};

//Lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(5) = {1, 2, 3,4};
Plane Surface(6) = {5};

// Physical groups
Physical Surface("Dominio", 1) = {6};
Physical Line("Apoio fixo", 2) = {2};  
Physical Line("Apoio movel",   3) = {1,3};
Physical Line("Carregamento",  4) = {4}; 

Transfinite Curve {4, 2} = 2; //Using Progression 1;
Transfinite Curve {1, 3} = 3; //Using Progression 1;

Transfinite Surface {6};
Recombine Surface {6};

