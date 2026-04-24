#include "GeoMesh.h"
#include "LinearAnalysis.h"
#include "MixedCompMesh.h"
#include "GmshTools.h"
#include "VTUGenerator.h"
#include "L2Projection.h"
#include "Elasticity2D.h"
#include "CompMeshTools.h"
#include "InterpolatedBC.h"
#include "EigenSpMatrix.h"
#include "ElementEnriched.h"
#include "Element.h"
#include "GlobalLocalEnrichment.h"
#include "MEFGGlobalLocalTools.h"

//Comentar erro ao debugar após incluir a weakform GlobalLocalEnrichment e o element ElementEnriched. 

// Defines the problem dimension
const int dimension = 2;
auto forcingFunction = [](const VecDouble &coord, VecDouble &force){
    const auto &x=coord[0];
    const auto &y=coord[1];
    //force[0] = 0;
    //force[1] = -1;
};

auto exactSol = [](const VecDouble &coord, VecDouble &u, MatrixDouble &gradU){
    const auto &x=coord[0];
    const auto &y=coord[1];

};

void CreateModel(CompMesh *cmesh);
void CreateModel2(CompMesh *cmesh);
void SolveProblem(CompMesh *cmesh,VecDouble &Solution);

double ModElasticity = 1.e6;
double PoissonRatio = 0.49999;

int main(int argc, char **args) { 

    //Create Model
    GeoMesh *gmesh = new GeoMesh();
    GmshTools::Read(*gmesh,"../chapaQuadrada.msh");
    CompMesh *cmesh = new CompMesh(gmesh,ApproxType::EIsoparametric);
    CreateModel(cmesh);
    cmesh->Print("cmesh.txt");

    CompMesh *cmeshaux = new CompMesh(gmesh,ApproxType::EIsoparametric);
    CreateModel2(cmeshaux);
    //gmesh->Print("gmesh.txt");
    cmeshaux->Print("cmeshAux.txt");
    
    VecDouble Solution;
    //Solve Problem
    SolveProblem(cmesh,Solution);

    //Compute strain energy
    LinearAnalysis an(cmeshaux,SolverType::ELDLt);
    an.Compute();   
    // an.PrintGlobalMatrix();

    EigenSpMatrix *spMat = dynamic_cast<EigenSpMatrix *>(an.GlobalMatrix());
    if (!spMat) {
        std::cerr << "Error: GlobalMatrix is not of type EigenSpMatrix." << std::endl;
    }
    auto globalMat = spMat->Matrix();

    // std::cout << "Global Matrix: \n" << globalMat << std::endl;

    VecDouble Force = globalMat * Solution;
    // std::cout << "Force: \n" << Force << std::endl;
    double strainEnergy = (Solution.dot(Force))/2;
    std::cout << std::fixed << std::setprecision(10) << "Strain Energy: "<< strainEnergy << std::endl;
    
}   

void CreateModel(CompMesh *cmesh){

    Elasticity2D *matelasticity = new Elasticity2D(1, ModElasticity, PoissonRatio, false, 1.0); //domínio global

    cmesh->InsertMaterial(matelasticity);
    
    //BC 
    MatrixDouble val1(2,2);
    val1.setZero();
    VecDouble val2(2);
    val2.setZero();

    //Chapa quadrada cisalhamento
    val2[0] = 0.001;
    L2Projection * matbc1 = new L2Projection(2,dimension-1,BoundaryConditionType::kDirectionalNonHomogeneousDirichlet,val1,val2);
    val2.setZero();
    val2[1] = 1.0;
    L2Projection * matbc2 = new L2Projection(3,dimension-1,BoundaryConditionType::kDirectionalHomogeneousDirichlet,val1,val2);
    val2.setZero();
    //val2[0] = 1.0;
    L2Projection * matbc3 = new L2Projection(4,dimension-1,BoundaryConditionType::kDirichlet,val1,val2);

    cmesh->InsertMaterial(matbc1);
    cmesh->InsertMaterial(matbc2);
    cmesh->InsertMaterial(matbc3);
    
    cmesh->AutoBuild();
}

void CreateModel2(CompMesh *cmesh){

    Elasticity2D *matelasticity = new Elasticity2D(1, ModElasticity, PoissonRatio, false, 1.0); //domínio global

    cmesh->InsertMaterial(matelasticity);
    
    // //BC 
    // MatrixDouble val1(2,2);
    // val1.setZero();
    // VecDouble val2(2);
    // val2.setZero();

    // //Chapa quadrada cisalhamento
    // val2[0] = 0.001;
    // L2Projection * matbc1 = new L2Projection(2,dimension-1,BoundaryConditionType::kDirichlet,val1,val2);
    // val2.setZero();
    // val2[1] = 1.0;
    // L2Projection * matbc2 = new L2Projection(3,dimension-1,BoundaryConditionType::kDirectionalHomogeneousDirichlet,val1,val2);
    // val2.setZero();
    // //val2[0] = 1.0;
    // L2Projection * matbc3 = new L2Projection(4,dimension-1,BoundaryConditionType::kDirichlet,val1,val2);

    // cmesh->InsertMaterial(matbc1);
    // cmesh->InsertMaterial(matbc2);
    // cmesh->InsertMaterial(matbc3);
    
    cmesh->AutoBuild();
}

void SolveProblem(CompMesh *cmesh, VecDouble &Solution){

    LinearAnalysis an(cmesh,SolverType::ELDLt);
       
    std::vector<std::string> ScalarNames, VectorNames;
    ScalarNames = {"SigmaX","SigmaY","TauXY","StrainEnergy"};
    VectorNames = {"Displacement"}; //, "ExactDisplacement"}

    an.Run();
    // an.PrintGlobalMatrix();
    // an.PrintSolution();
    // an.PrintGlobalRhs();

    //VecDouble errors(4);
    //anG.PostProcessError(errors);

    EigenSpMatrix *spMat = dynamic_cast<EigenSpMatrix *>(an.GlobalMatrix());
    if (!spMat) {
        std::cerr << "Error: GlobalMatrix is not of type EigenSpMatrix." << std::endl;
        return;
    }

    VecDouble sol = spMat->Solution();
    VecDouble rhs = spMat->Rhs();
    Solution = sol;
    

    //double strainEnergy = (sol.dot(rhs))/2;
    //std::cout << std::fixed << std::setprecision(10) << "Strain Energy: "<< strainEnergy << std::endl;

    VTUGenerator::PrintResults(cmesh,"mefResult",ScalarNames,VectorNames);
}
