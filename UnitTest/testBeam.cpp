#include <fstream>
#include <LinearAnalysis.h>
#include <GmshTools.h>
#include <VTUGenerator.h>
#include <L2Projection.h>
#include <LinearBeam.h>
#include <memory>

auto forcing = [](const VecDouble &coord, VecDouble &force){
    const auto &x=coord[0];
    const auto &y=coord[1];
    force[0] = -30000;
};

int main()
{
    //Geometric Mesh
    GeoMesh * gmesh = new GeoMesh();
    GmshTools::Read(*gmesh,"../../UnitTest/viga.msh");
    gmesh->Print("gmesh.txt");

    CompMesh* cmesh = new CompMesh(gmesh,ApproxType::EIsoparametric);
    cmesh->Dimension() = 1;

    LinearBeam * matelas = new LinearBeam(6,1e11,8.333333e-02);
    matelas->SetForcingFunction(forcing);
    cmesh->InsertMaterial(matelas);

    //BC
    MatrixDouble val1(2,2);
    val1.setZero();
    VecDouble val2(2);
    val2.setZero();
    
    //Engaste
    L2Projection * matbc1 = new L2Projection(3,0,BoundaryConditionType::kDirichlet,val1,val2);

    //Apoio móvel
    val2[0]=1.;
    L2Projection * matbc2 = new L2Projection(4,0,BoundaryConditionType::kDirectionalHomogeneousDirichlet,val1,val2);
    L2Projection * matbc3 = new L2Projection(5,0,BoundaryConditionType::kDirectionalHomogeneousDirichlet,val1,val2);

    cmesh->InsertMaterial(matbc1);
    cmesh->InsertMaterial(matbc2);
    cmesh->InsertMaterial(matbc3);
 
    cmesh->AutoBuild();
    cmesh->Print("cmesh.txt");

    LinearAnalysis an(cmesh,SolverType::ELDLt);
   
    std::vector<std::string> ScalarNames, VectorNames;
    ScalarNames = {"Displacement"};
    // ScalarNames = {"Stress","Strain","AxialForce"};
    an.Run();

    an.PrintGlobalMatrix();
    an.PrintGlobalRhs();
    an.PrintSolution();

    //Put a check criterion here.

    VTUGenerator::PrintResultsHierarquic(cmesh,"resultBeam",ScalarNames,VectorNames); 

    return 0;
}