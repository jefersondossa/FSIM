// Gmsh Script (.geo) - Tetraedro Cúbico (20 nós / Ordem 3)

// 1. Definição dos Vértices do Tetraedro no Espaço de Referência
// Point(tag) = {x, y, z, tamanho_caracteristico_da_malha};
Point(1) = {0.0, 0.0, 1.0, 10.0};
Point(2) = {1.0, 0.0, 0.0, 10.0};
Point(3) = {0.0, 0.0, 0.0, 10.0};
Point(4) = {0.0, 1.0, 0.0, 10.0};

// 2. Definição das Arestas (Linhas)
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 1};
Line(4) = {1, 4};
Line(5) = {2, 4};
Line(6) = {3, 4};

// 3. Definição dos Contornos de Face (Curve Loops) e Superfícies
Curve Loop(1) = {1, 2, 3};
Plane Surface(1) = {1}; // Base z=0

Curve Loop(2) = {1, 5, -4};
Plane Surface(2) = {2}; // Face y=0

Curve Loop(3) = {2, 6, -5};
Plane Surface(3) = {3}; // Face inclinada

Curve Loop(4) = {3, 4, -6};
Plane Surface(4) = {4}; // Face x=0

// 4. Definição do Volume do Tetraedro
Surface Loop(1) = {1, 2, 3, 4};
Volume(1) = {1};

// Opcional: Define grupo físico para o volume
Physical Surface("face", 13) = {1,2,3,4};
Physical Volume("TetrahedronDomain", 14) = {1};
