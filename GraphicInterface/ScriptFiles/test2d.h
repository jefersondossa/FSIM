// Defines the problem dimension
const int dimension = 2;

{
auto yieldFunction = [](const double &plast, double &sigma_y, double &hardening){
hardening = 500.;
sigma_y =.5+hardening*plast;
};

CompMesh* cmesh = new CompMesh();

//Physical Group Right
Elasticity2D * Right = new Elasticity2D("matid", 0, 0, false);

cmesh->InsertMaterial(Right);

//Physical Group Top
NavierStokes * Top = new NavierStokes("matid", "dimension", 0, 0);

cmesh->InsertMaterial(Top);

//Physical Group Left
L2Projection * Left = new L2Projection("matid", "dimension", Neumann, 0, 0, 0);

cmesh->InsertMaterial(Left);

};
