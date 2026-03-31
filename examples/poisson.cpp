#include "GeoMesh.h"
#include "LinearAnalysis.h"
#include "MixedCompMesh.h"
#include "GmshTools.h"
#include "VTUGenerator.h"
#include "L2Projection.h"
#include "Poisson.h"

auto exactSol = [](const VecDouble &coord, VecDouble &u, MatrixDouble &gradU){
    const auto &x=coord[0];
    const auto &y=coord[1];

    u[0] = .5*(x*x-1);
    gradU(0,0) = x;
    gradU(1,0) = 0.;
};

// Defines the problem dimension
const int dimension = 1;
auto forcing = [](const VecDouble &coord, VecDouble &force){
    const auto &x=coord[0];
    const auto &y=coord[1];
    force[0] = -1;
};

int main(int argc, char **args) { 

    //Geometric Mesh
    GeoMesh * gmesh = new GeoMesh();
    GmshTools::Read(*gmesh,"../hierarquic2d.msh");
    gmesh->Print("gmesh.txt");

    CompMesh* cmesh = new CompMesh(gmesh,ApproxType::EIsoparametric);


    Poisson * matelas = new Poisson(8,1,1);
    matelas->SetForcingFunction(forcing);
    matelas->SetExactSolution(exactSol);
    cmesh->InsertMaterial(matelas);

    //BC
    MatrixDouble val1(1,1);
    val1.setZero();
    VecDouble val2(1);
    val2.setZero();
    L2Projection * matbc1 = new L2Projection(5,0,BoundaryConditionType::kDirichlet,val1,val2);
    L2Projection * matbc2 = new L2Projection(6,0,BoundaryConditionType::kDirichlet,val1,val2);
    L2Projection * matbc3 = new L2Projection(7,0,BoundaryConditionType::kDirichlet,val1,val2);
    matbc1->SetExactSolution(exactSol);
    matbc2->SetExactSolution(exactSol);
    matbc3->SetExactSolution(exactSol);

    cmesh->InsertMaterial(matbc1);
    cmesh->InsertMaterial(matbc2);
    cmesh->InsertMaterial(matbc3);
 
    cmesh->AutoBuild();
    cmesh->Print("cmesh.txt");


    LinearAnalysis an(cmesh,SolverType::ELU);
   
    std::vector<std::string> ScalarNames, VectorNames;
    // VectorNames = {"Solution"};
    ScalarNames = {"Solution","ExactSolution"};
    an.Run();


    VecDouble errors(3);
    an.PostProcessError(errors);

    VTUGenerator::PrintResults(cmesh,"poisson",ScalarNames,VectorNames);    
        
}           