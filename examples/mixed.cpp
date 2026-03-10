#include "GeoMesh.h"
#include "LinearAnalysis.h"
#include "MixedCompMesh.h"
#include "GmshTools.h"
#include "VTUGenerator.h"
#include "L2Projection.h"
#include "MixedElasticity.h"

CompMesh* CreateDisplacementMesh(GeoMesh *gmesh);
CompMesh* CreatePressureMesh(GeoMesh *gmesh);
MixedCompMesh* CreateMixedMesh(std::vector<CompMesh *> &meshvector);

int main(int argc, char **args) {   
    //Geometric Mesh
    GeoMesh * gmesh = new GeoMesh();
    GmshTools::Read(*gmesh,"../hierarquic2d.msh");
    gmesh->Print("gmesh.txt");

    //Create displacement mesh
    CompMesh* cmeshdisp = CreateDisplacementMesh(gmesh);
    
    //Create pressure mesh 
    CompMesh* cmeshpressure = CreatePressureMesh(gmesh);

    std::vector<CompMesh *> meshvector = {cmeshdisp, cmeshpressure};

    MixedCompMesh *cmesh = CreateMixedMesh(meshvector);

    LinearAnalysis an(cmesh,SolverType::ELDLt);
    an.Run();
    an.PrintGlobalMatrix();
    an.PrintGlobalRhs();
    an.PrintSolution();

    std::vector<std::string> ScalarNames, VectorNames;
    // ScalarNames = {"Displacement"};
    VectorNames = {"Displacement"};

    VTUGenerator::PrintResults(cmesh,"mixed",ScalarNames,VectorNames);

    // VecDouble errors;
    // an.PostProcessError(errors);

} 


CompMesh* CreateDisplacementMesh(GeoMesh *gmesh){

    CompMesh* cmesh = new CompMesh(gmesh, ApproxType::EHierarquic); 
    cmesh->SetDefaultOrder(2);

    int nstate = 2;
    WeakForm * mat = new WeakForm(8,nstate);
    cmesh->InsertMaterial(mat);

    //BC
    MatrixDouble val1(nstate,nstate);
    val1.setZero();
    VecDouble val2(nstate);
    VecDouble val3(nstate);
    val2.setZero();val3.setZero();
    // val2[0] = 1.0;
    L2Projection * matbc3 = new L2Projection(6,1,BoundaryConditionType::kDirichlet,val1,val2);
    val2.setZero();
    L2Projection * matbc1 = new L2Projection(5,1,BoundaryConditionType::kDirichlet,val1,val2);
    L2Projection * matbc2 = new L2Projection(7,1,BoundaryConditionType::kDirichlet,val1,val3);

    cmesh->InsertMaterial(matbc1);
    cmesh->InsertMaterial(matbc2);
    cmesh->InsertMaterial(matbc3);

    cmesh->AutoBuild();
    
    cmesh->Print("cmesh_disp.txt");

    return cmesh;

}



CompMesh* CreatePressureMesh(GeoMesh *gmesh){

    CompMesh* cmesh = new CompMesh(gmesh, ApproxType::EHierarquic); 
    cmesh->SetDefaultOrder(1);

    int nstate = 1;
    WeakForm * mat = new WeakForm(8,nstate);
    cmesh->InsertMaterial(mat);

    //BC
    MatrixDouble val1(nstate,nstate);
    val1.setZero();
    VecDouble val2(nstate);
    val2.setZero();
    L2Projection * matbc3 = new L2Projection(6,1,BoundaryConditionType::kNeumann,val1,val2);
    L2Projection * matbc1 = new L2Projection(5,1,BoundaryConditionType::kNeumann,val1,val2);
    L2Projection * matbc2 = new L2Projection(7,1,BoundaryConditionType::kNeumann,val1,val2);

    cmesh->InsertMaterial(matbc1);
    cmesh->InsertMaterial(matbc2);
    cmesh->InsertMaterial(matbc3);

    cmesh->AutoBuild();
    
    cmesh->Print("cmesh_pressure.txt");

    return cmesh;
}

MixedCompMesh* CreateMixedMesh(std::vector<CompMesh *> &meshvector){

    MixedCompMesh* cmesh = new MixedCompMesh(meshvector);

    MixedElasticity * mat = new MixedElasticity(8,2,1.0,0.0);
    cmesh->InsertMaterial(mat);

    cmesh->AutoBuild();

    cmesh->Print("cmesh_mixed.txt");

    return cmesh;
};
