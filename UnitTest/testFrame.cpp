#include <fstream>
#include <LinearAnalysis.h>
#include <GmshTools.h>
#include <VTUGenerator.h>
#include <L2Projection.h>
#include <LinearFrame.h>
#include <memory>

auto forcing = [](const VecDouble &coord, VecDouble &force){
    const auto &x=coord[0];
    const auto &y=coord[1];
    //force[1] = -10000;
};

auto forcing_concentrate = [](const VecDouble &coord, VecDouble &force){
    const auto &x=coord[0];
    const auto &y=coord[1];
    force[0] = 40000;
};

int main()
{
    //Malha geométrica (arquivo .msh)
    GeoMesh * gmesh = new GeoMesh();
    GmshTools::Read(*gmesh,"../../UnitTest/portico.msh");
    gmesh->Print("gmesh.txt");

    CompMesh* cmesh = new CompMesh(gmesh,ApproxType::EIsoparametric);
    cmesh->Dimension() = 1;

    LinearFrame * matelas = new LinearFrame(13,1e11,8.333333e-02, 1.);
    matelas->SetForcingFunction(forcing);
    cmesh->InsertMaterial(matelas);

    LinearFrame * matelas2 = new LinearFrame(14,1e11,8.333333e-02, 1.);
    cmesh->InsertMaterial(matelas2);

    // BC
    MatrixDouble val1(3,3);
    val1.setZero();
    VecDouble val2(3);
    val2.setZero();

    // Engaste Esquerdo
    L2Projection * matbc1 = new L2Projection(10,0,BoundaryConditionType::kDirichlet,val1,val2);

    // Engaste Direito
    L2Projection * matbc2 = new L2Projection(11,0,BoundaryConditionType::kDirichlet,val1,val2);

    // Força concentrada
    L2Projection * matbc3 = new L2Projection(12,0,BoundaryConditionType::kNeumann,val1,val2);
    matbc3->SetForcingFunction(forcing_concentrate);

    cmesh->InsertMaterial(matbc1);
    cmesh->InsertMaterial(matbc2);
    cmesh->InsertMaterial(matbc3);
 
    cmesh->AutoBuild();
    cmesh->Print("cmesh.txt");

    LinearAnalysis an(cmesh,SolverType::ELDLt);
   
    std::vector<std::string> ScalarNames, VectorNames;
    VectorNames = {"Displacement"};
    // ScalarNames = {"Stress","Strain","AxialForce"};
    an.Run();

    an.PrintGlobalMatrix();
    an.PrintGlobalRhs();
    an.PrintSolution();

    //Put a check criterion here.

    VTUGenerator::PrintResultsHierarquic(cmesh,"resultFrame",ScalarNames,VectorNames); 

    return 0;
}