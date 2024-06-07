// Defines the problem dimension
const int dimension = d

{
auto yieldFunction = [](const double &plast, double &sigma_y, double &hardening){
hardening = 500.;
sigma_y =.5+hardening*plast;
};

CompMesh* cmesh = new CompMesh();

//Physical Group Left

Stokes * Left = new Stokes("matid", "dimension", 0, 0);

cmesh->InsertMaterial((Left);

//Physical Group 

//Physical Group Right

