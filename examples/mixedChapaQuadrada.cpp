#include "GeoMesh.h"
#include "LinearAnalysis.h"
#include "MixedCompMesh.h"
#include "GmshTools.h"
#include "VTUGenerator.h"
#include "L2Projection.h"
#include "EigenSpMatrix.h"
#include "MixedElasticity.h"

CompMesh* CreateDisplacementMesh(GeoMesh *gmesh);
CompMesh* CreatePressureMesh(GeoMesh *gmesh);
MixedCompMesh* CreateMixedMesh(std::vector<CompMesh *> &meshvector);

// Defines the problem dimension
const int dimension = 2;

void SolveProblem(CompMesh *cmesh);

int main(int argc, char **args) {   
    //Geometric Mesh
    GeoMesh * gmesh = new GeoMesh();
    GmshTools::Read(*gmesh,"../chapaQuadrada.msh");
    gmesh->Print("gmesh.txt");

    //Create displacement mesh
    CompMesh* cmeshdisp = CreateDisplacementMesh(gmesh);
    
    //Create pressure mesh 
    CompMesh* cmeshpressure = CreatePressureMesh(gmesh);

    std::vector<CompMesh *> meshvector = {cmeshdisp, cmeshpressure};

    MixedCompMesh *cmesh = CreateMixedMesh(meshvector);

    SolveProblem(cmesh);

    //cmesh->Print("cmesh_mixed2.txt");
    //cmeshpressure->Print("cmesh_press2.txt");
    //cmeshdisp->Print("cmesh_disp2.txt");
} 

void SolveProblem(CompMesh *cmesh){

    LinearAnalysis an(cmesh,SolverType::ELDLt);
    an.Run();
    // an.PrintGlobalMatrix();
    // an.PrintGlobalRhs();
    // an.PrintSolution();

    EigenSpMatrix *spMat = dynamic_cast<EigenSpMatrix *>(an.GlobalMatrix());
    if (!spMat) {
        std::cerr << "Error: GlobalMatrix is not of type EigenSpMatrix." << std::endl;
        return;
    }
    
    VecDouble sol = spMat->Solution();
    VecDouble rhs = spMat->Rhs();
    double strainEnergy = (sol.dot(rhs))/2;
    double strainEnergy2 = 0.0;
    for (int i = 0; i < sol.size(); i++){
        if (fabs(sol[i])>1.e3 || fabs(rhs[i])>1.e3)continue;
        strainEnergy2 += sol[i]*rhs[i]/2.;
    }
    std::cout << std::fixed << std::setprecision(10) << "Strain Energy: "<< strainEnergy << std::endl;
    std::cout << std::fixed << std::setprecision(10) << "Strain Energy 2: "<< strainEnergy2 << std::endl;

    std::vector<std::string> ScalarNames, VectorNames;
    ScalarNames = {"Pressure"};
    VectorNames = {"Displacement"};

    VTUGenerator::PrintResults(cmesh,"mixed",ScalarNames,VectorNames);

    // VecDouble errors;
    // an.PostProcessError(errors);
}

CompMesh* CreateDisplacementMesh(GeoMesh *gmesh){

    CompMesh* cmesh = new CompMesh(gmesh, ApproxType::EHierarquic); 
    cmesh->SetDefaultOrder(2);

    int nstate = 2;
    WeakForm * mat = new WeakForm(1, nstate);
    cmesh->InsertMaterial(mat);

    //BC
    MatrixDouble val1(nstate,nstate);
    val1.setZero();

    VecDouble val2(nstate);
    val2.setZero();
    
    val2[0] = 0.001;
    L2Projection * matbc1 = new L2Projection(2,dimension-1,BoundaryConditionType::kDirectionalNonHomogeneousDirichlet,val1,val2);
    val2[1] = 1.0;
    L2Projection * matbc2 = new L2Projection(3,dimension-1,BoundaryConditionType::kDirectionalHomogeneousDirichlet,val1,val2);
    val2.setZero();
    val2[0] = 1.0;
    L2Projection * matbc3 = new L2Projection(4,dimension-2,BoundaryConditionType::kDirectionalHomogeneousDirichlet,val1,val2);

    cmesh->InsertMaterial(matbc1);
    cmesh->InsertMaterial(matbc2);
    cmesh->InsertMaterial(matbc3);

    cmesh->AutoBuild(); 

    return cmesh;
}

CompMesh* CreatePressureMesh(GeoMesh *gmesh){

    CompMesh* cmesh = new CompMesh(gmesh, ApproxType::EHierarquic); 
    cmesh->SetDefaultOrder(1);

    int nstate = 1;
    WeakForm * mat = new WeakForm(1, nstate);
    cmesh->InsertMaterial(mat);

    //BC
    MatrixDouble val1(nstate,nstate);
    val1.setZero();
    VecDouble val2(nstate);
    val2.setZero();
    L2Projection * matbc1 = new L2Projection(2,1,BoundaryConditionType::kNeumann,val1,val2);
    L2Projection * matbc2 = new L2Projection(3,1,BoundaryConditionType::kNeumann,val1,val2);
    L2Projection * matbc3 = new L2Projection(4,1,BoundaryConditionType::kNeumann,val1,val2);

    cmesh->InsertMaterial(matbc1);
    cmesh->InsertMaterial(matbc2);
    cmesh->InsertMaterial(matbc3);

    //cmesh->CreateDisconnectedElements();
    cmesh->AutoBuild();
    
    cmesh->Print("cmesh_pressure.txt");

    return cmesh;
}

MixedCompMesh* CreateMixedMesh(std::vector<CompMesh *> &meshvector){

    MixedCompMesh* cmesh = new MixedCompMesh(meshvector);

    MixedElasticity * mat = new MixedElasticity(1, 2, 1.0, 0.3);
    cmesh->InsertMaterial(mat);

    cmesh->AutoBuild();

    cmesh->Print("cmesh_mixed.txt");

    return cmesh;
};
