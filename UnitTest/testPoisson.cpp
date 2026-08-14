#include <fstream>
#include <LinearAnalysis.h>
#include <GmshTools.h>
#include <VTUGenerator.h>
#include <L2Projection.h>
#include <Poisson.h>
#include <memory>

auto forcing = [](const VecDouble &coord, VecDouble &force){
    const auto &x=coord[0];
    const auto &y=coord[1];
    force[0] = std::sin(M_PI*x/1);
};

auto exactSol = [](const VecDouble &coord, VecDouble &u, MatrixDouble &gradU){
    const auto &x=coord[0];
    const auto &y=coord[1];

    u[0] = std::sin(M_PI*x/1);
    gradU(0,0) = std::cos(M_PI*x/1)*1*1/(M_PI * M_PI);
};

int main()
{
    //Geometric Mesh
    GeoMesh * gmesh = new GeoMesh();
    GmshTools::Read(*gmesh,"../../UnitTest/poisson1d.msh");
    gmesh->Print("gmesh.txt");

    CompMesh* cmesh = new CompMesh(gmesh,ApproxType::EIsoparametric);
    cmesh->Dimension() = 1;

    Poisson * mat = new Poisson(6,1,1);
    mat->SetForcingFunction(forcing);
    mat->SetExactSolution(exactSol);
    cmesh->InsertMaterial(mat);

    //BC
    MatrixDouble val1(1,1);
    val1.setZero();
    VecDouble val2(1);
    val2.setZero();
    
    //Apoio fixo
    L2Projection * matbc1 = new L2Projection(3,0,BoundaryConditionType::kDirichlet,val1,val2);

    //Apoio fixo
    L2Projection * matbc2 = new L2Projection(4,0,BoundaryConditionType::kDirichlet,val1,val2);

    cmesh->InsertMaterial(matbc1);
    cmesh->InsertMaterial(matbc2);
 
    cmesh->AutoBuild();
    cmesh->Print("cmesh.txt");

    LinearAnalysis an(cmesh,SolverType::ELDLt);
   
    std::vector<std::string> ScalarNames, VectorNames;
    ScalarNames = {"Solution"};
    VectorNames = {"Derivative"};
    an.Run();

    //Put a check criterion here.

    VTUGenerator::PrintResults(cmesh,"resultPoisson",ScalarNames,VectorNames); 

    VecDouble errors(3);
    an.PostProcessError(errors);

    return 0;
}