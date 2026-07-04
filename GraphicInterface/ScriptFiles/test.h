// Defines the problem dimension
const int dimension = 2;

{
auto yieldFunction = [](const REAL &plast, REAL &sigma_y, REAL &hardening){
hardening = 500.;
sigma_y =.5+hardening*plast;
};

CompMesh* cmesh = new CompMesh();

//Physical Group Right
L2Projection * Right = new L2Projection("matid", "dimension", Neumann, 0, -9, 0);

cmesh->InsertMaterial(Right);

//Physical Group Domain
Elasticity2D * Domain = new Elasticity2D("matid", 0, 0, false);

cmesh->InsertMaterial(Domain);

};
