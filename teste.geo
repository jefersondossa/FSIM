tm = 1; //tamanho da malha
tmr = 0.1; //tamanho refinado
r = 3; //raio da circunsferência

//Pontos
Point(1) = {0, 0, 0, tm};
Point(2) = {40, 0, 0, tm};
Point(3) = {40, 20, 0, tm};
Point(4) = {0, 20, 0, tm};

Point(5) = {20, 10, 0};
Point(6) = {20+r, 10, 0, tmr};
Point(7) = {20-r, 10, 0, tmr};

//Linhas
Line(1) = {1, 2}; //Transfinite Line {1} = 40 Using Progression 1.5;
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Circle(5) = {6, 5, 7}; //{inicio, centro, fim} 
Circle(6) = {7, 5, 6};

//Contornos
Curve Loop(1) = {1, 2, 3, 4}; //contorno retangular
Curve Loop(2) = {5, 6}; //contorno circular

//Superficie
Plane Surface(1) = {1, 2}; //{contorno exterior, contorno interior}

//Grupos físicos
Physical Curve("Borda direita", 1) = {2}; //borda direita (Neumann) = lado 2
Physical Curve("Borda esquerda", 2) = {4}; //borda esquerda (Dirichlet) = lado 4
Physical Surface("Superficie", 3) = {1}; //{plane surface}

Mesh 2; //dimensão da malha

//Mesh.SurfaceFaces = 1;
//Mesh.Points = 1;

//Transfinite Surface {1}; //malha estruturada {plane surface}
//Recombine Surface {1}; //malha com elementos quadrangulares {plane surface}
//Transfinite Curve {1, 2, 3, 4} = 5;//divide os 4 lados em 5 nós 

