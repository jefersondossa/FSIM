#include <fstream>
#include <LinearAnalysis.h>
#include <GmshTools.h>
#include <VTUGenerator.h>
#include <L2Projection.h>
#include <ElasticTruss.h>
#include <memory>

auto forcing = [](const VecDouble &coord, VecDouble &force){
    const auto &x=coord[0];
    const auto &y=coord[1];
    force[1] = -5000.;
    force[0] = 10000;
};

int main()
{
    //Geometric Mesh
    GeoMesh * gmesh = new GeoMesh();
    GmshTools::Read(*gmesh,"../../UnitTest/trelica1.msh");
    gmesh->Print("gmesh.txt");

    CompMesh* cmesh = new CompMesh(gmesh,ApproxType::EIsoparametric);
    cmesh->Dimension() = 1;

    ElasticTruss * matelas = new ElasticTruss(11,2,100000e6,7.8540e-05);
    cmesh->InsertMaterial(matelas);

    //BC
    MatrixDouble val1(2,2);
    val1.setZero();
    VecDouble val2(2);
    val2.setZero();
    //Apoio fixo
    L2Projection * matbc1 = new L2Projection(8,0,BoundaryConditionType::kDirichlet,val1,val2);

    //Apoio móvel
    val2[0]=1.;
    L2Projection * matbc2 = new L2Projection(9,0,BoundaryConditionType::kDirectionalHomogeneousDirichlet,val1,val2);
    
    //Forca
    L2Projection * matbc3 = new L2Projection(10,0,BoundaryConditionType::kNeumann,val1,val2);
    matbc3->SetForcingFunction(forcing);

    cmesh->InsertMaterial(matbc1);
    cmesh->InsertMaterial(matbc2);
    cmesh->InsertMaterial(matbc3);
 
    cmesh->AutoBuild();
    cmesh->Print("cmesh.txt");

    LinearAnalysis an(cmesh,SolverType::ELDLt);
   
    std::vector<std::string> ScalarNames, VectorNames;
    VectorNames = {"Displacement"};
    ScalarNames = {"Stress","Strain","AxialForce"};
    an.Run();

    an.PrintGlobalMatrix();
    an.PrintGlobalRhs();
    an.PrintSolution();

    VTUGenerator::PrintResultsHierarquic(cmesh,"resultTruss",ScalarNames,VectorNames); 

    return 0;
}