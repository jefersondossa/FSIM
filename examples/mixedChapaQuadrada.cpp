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

void SolveProblem(CompMesh *cmesh, VecDouble &solution);
void CreateDAux(CompMesh *cmesh);
void CreatePAux(CompMesh *cmesh);

REAL ModElasticity = 1.e6;
REAL PoissonRatio = 0.3;

int main(int argc, char **args) {   
    //Geometric Mesh
    GeoMesh * gmesh = new GeoMesh();
    GmshTools::Read(*gmesh,"../chapaQuadrada.msh");
    gmesh->Print("gmesh.txt");

    //Create displacement mesh
    CompMesh* cmeshdisp = CreateDisplacementMesh(gmesh);
    //cmeshdisp->Print("cmesh.txt");

    //Create aux displacement mesh
    CompMesh *cmeshDaux = new CompMesh(gmesh,ApproxType::EHierarquic);
    cmeshDaux->SetDefaultOrder(2);
    CreateDAux(cmeshDaux);
    //cmeshDaux->Print("cmeshAux.txt");
    
    //Create pressure mesh 
    CompMesh* cmeshpressure = CreatePressureMesh(gmesh);

    //Create aux pressure mesh
    CompMesh *cmeshPaux = new CompMesh(gmesh,ApproxType::EHierarquic);
    cmeshPaux->SetDefaultOrder(0);
    CreatePAux(cmeshPaux);
    //cmeshPaux->Print("cmeshPAux.txt");

    std::vector<CompMesh *> meshvector = {cmeshdisp, cmeshpressure};
    std::vector<CompMesh *> meshvectorAux = {cmeshDaux, cmeshPaux};

    MixedCompMesh *cmesh = CreateMixedMesh(meshvector);
    MixedCompMesh *cmeshaux = CreateMixedMesh(meshvectorAux);

    VecDouble solution;
    SolveProblem(cmesh, solution);

    //Compute strain energy
    LinearAnalysis an(cmeshaux,SolverType::ELDLt);
    an.Compute();   
    //an.PrintGlobalMatrix();

    EigenSpMatrix *spMat = dynamic_cast<EigenSpMatrix *>(an.GlobalMatrix());
    if (!spMat) {
        std::cerr << "Error: GlobalMatrix is not of type EigenSpMatrix." << std::endl;
    }
    auto globalMat = spMat->Matrix();

    // std::cout << "Global Matrix: \n" << globalMat << std::endl;

    VecDouble Force = globalMat * solution;
    //std::cout << "Force: \n" << Force << std::endl;
    //std::cout << "solution: \n" << solution << std::endl;
    REAL strainEnergy = (solution.dot(Force))/2;
    std::cout << std::fixed << std::setprecision(10) << "Strain Energy: "<< strainEnergy << std::endl;
    
} 

void SolveProblem(CompMesh *cmesh, VecDouble &solution){

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
    solution = sol;

    std::vector<std::string> ScalarNames, VectorNames;
    ScalarNames = {"Pressure"};
    VectorNames = {"Displacement"};

    VTUGenerator::PrintResults(cmesh,"mixed",ScalarNames,VectorNames);

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
    
    //Chapa quadrada cisalhamento
    val2[0] = 0.001;
    L2Projection * matbc1 = new L2Projection(2,dimension-1,BoundaryConditionType::kDirectionalNonHomogeneousDirichlet,val1,val2);
    val2.setZero();
    val2[1] = 1.0;
    L2Projection * matbc2 = new L2Projection(3,dimension-1,BoundaryConditionType::kDirectionalHomogeneousDirichlet,val1,val2);
    //L2Projection * matbc2 = new L2Projection(3,dimension-1,BoundaryConditionType::kNeumann,val1,val2);
    val2.setZero();
    L2Projection * matbc3 = new L2Projection(4,dimension-1,BoundaryConditionType::kDirichlet,val1,val2);

    cmesh->InsertMaterial(matbc1);
    cmesh->InsertMaterial(matbc2);
    cmesh->InsertMaterial(matbc3);

    cmesh->AutoBuild(); 

    return cmesh;
}

void CreateDAux(CompMesh *cmesh){

    int nstate = 2;
    WeakForm * mat = new WeakForm(1, nstate);
    cmesh->InsertMaterial(mat);
    
    cmesh->AutoBuild();
}

CompMesh* CreatePressureMesh(GeoMesh *gmesh){

    CompMesh* cmesh = new CompMesh(gmesh, ApproxType::EHierarquic); 
    cmesh->SetDefaultOrder(0);

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

void CreatePAux(CompMesh *cmesh){

    int nstate = 1;
    WeakForm * mat = new WeakForm(1, nstate);
    cmesh->InsertMaterial(mat);
    
    cmesh->AutoBuild();
}

MixedCompMesh* CreateMixedMesh(std::vector<CompMesh *> &meshvector){

    MixedCompMesh* cmesh = new MixedCompMesh(meshvector);

    MixedElasticity * mat = new MixedElasticity(1, dimension, ModElasticity, PoissonRatio);
    cmesh->InsertMaterial(mat);

    cmesh->AutoBuild();

    cmesh->Print("cmesh_mixed.txt");

    return cmesh;
};
