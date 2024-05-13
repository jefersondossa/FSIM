// Defines the problem dimension
const int dimension = d

{
auto yieldFunction = [](const double &plast, double &sigma_y, double &hardening){
hardening = 500.;
sigma_y =.5+hardening*plast;
};

CompMesh* cmesh = new CompMesh();

//Physical Group 

//Physical Group Right

Elasticity2D * Right = new Elasticity2D("matid", 0, 0, 0);

cmesh->InsertMaterial((Right);
//Physical Group Left

ElasticityPositional2D * Left = new ElasticityPositional2D("matid", 0, 0, 0);

cmesh->InsertMaterial((Left);
//Physical Group Top

Stokes * Top = new Stokes("matid", "dimension", 0, 0);

cmesh->InsertMaterial((Top);

//Physical Group Domain

Stokes * Domain = new Stokes("matid", "dimension", 0, 0);

cmesh->InsertMaterial((Domain);

};
