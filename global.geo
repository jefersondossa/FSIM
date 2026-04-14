//Points
Point(1) = {0, 0, 0, 1.0};
Point(2) = {2, 0, 0, 1.0};
Point(3) = {4, 0, 0, 1.0};
Point(4) = {4, 2, 0, 1.0};
Point(5) = {2, 2, 0, 1.0};
Point(6) = {0, 2, 0, 1.0};

//Lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};
Line(7) = {2, 5};

Line Loop(1) = {2, 3, 4, -7};
Line Loop(2) = {1, 7, 5, 6};

Plane Surface(1) = {1};
Plane Surface(2) = {2};

// Physical groups
Physical Surface("Dominio Global", 1) = {1};
Physical Surface("Dominio Local", 2) = {2};
Physical Line("Apoio fixo", 3) = {3};  
Physical Line("Apoio movel", 4) = {1,2,4,5};
Physical Line("Carregamento", 5) = {6}; 

Transfinite Curve {2, 3, 4, 6, 7} = 2; //Using Progression 1;
Transfinite Curve {1, 5} = 2; //Using Progression 1;

Transfinite Surface {1};
Recombine Surface {1};

Transfinite Surface {2};
Recombine Surface {2};
