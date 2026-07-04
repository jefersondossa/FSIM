#include "DataTypes.h"
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
auto forcingFunction1 = [](const VecDouble &coord, VecDouble &force){
    const auto &x=coord[0];
    const auto &y=coord[1];

    force[0] = 24*y - 120;
    force[1] = 0;
};
auto forcingFunction2 = [](const VecDouble &coord, VecDouble &force){
    const auto &x=coord[0];
    const auto &y=coord[1];

    force[0] = -24*y + 120;
    force[1] = 0;
};
auto forcingFunctionAB = [](const VecDouble &coord, VecDouble &force){
    const auto &x=coord[0];
    const auto &y=coord[1];
    
    REAL A = 1.0;
    REAL Q = 0.543075579;
    REAL lambda = 0.544483737;
    REAL r = sqrt(x*x + y*y);
    REAL theta = atan2(y, x);

    MatrixDouble stress(2,2);
    //stress x
    stress(0,0) = A*lambda*pow(r, lambda-1)*((2 - Q*(lambda+1))*cos((lambda-1)*theta) - (lambda-1)*cos((lambda-3)*theta));
    //stress y
    stress(1,1) = A*lambda*pow(r, lambda-1)*((2 + Q*(lambda+1))*cos((lambda-1)*theta) + (lambda-1)*cos((lambda-3)*theta));
    //stress xy
    stress(0,1) = A*lambda*pow(r, lambda-1)*((lambda-1)*sin((lambda-3)*theta) + Q*(lambda+1)*sin((lambda-1)*theta));
    stress(1,0) = stress(0,1);

    VecDouble n(2);
    n[0] = sqrt(2)/2;
    n[1] = sqrt(2)/2;

    force = stress * n;
};

auto forcingFunctionBC = [](const VecDouble &coord, VecDouble &force){
    const auto &x=coord[0];
    const auto &y=coord[1];
    
    REAL A = 1.0;
    REAL Q = 0.543075579;
    REAL lambda = 0.544483737;
    REAL r = sqrt(x*x + y*y);
    REAL theta = atan2(y, x);

    MatrixDouble stress(2,2);
    //stress x
    stress(0,0) = A*lambda*pow(r, lambda-1)*((2 - Q*(lambda+1))*cos((lambda-1)*theta) - (lambda-1)*cos((lambda-3)*theta));
    //stress y
    stress(1,1) = A*lambda*pow(r, lambda-1)*((2 + Q*(lambda+1))*cos((lambda-1)*theta) + (lambda-1)*cos((lambda-3)*theta));
    //stress xy
    stress(0,1) = A*lambda*pow(r, lambda-1)*((lambda-1)*sin((lambda-3)*theta) + Q*(lambda+1)*sin((lambda-1)*theta));
    stress(1,0) = stress(0,1);

    VecDouble n(2);
    n[0] = - sqrt(2)/2;
    n[1] = sqrt(2)/2;

    force = stress * n;
};

auto forcingFunctionEF = [](const VecDouble &coord, VecDouble &force){
    const auto &x=coord[0];
    const auto &y=coord[1];
    
    REAL A = 1.0;
    REAL Q = 0.543075579;
    REAL lambda = 0.544483737;
    REAL r = sqrt(x*x + y*y);
    REAL theta = atan2(y, x);

    MatrixDouble stress(2,2);
    //stress x
    stress(0,0) = A*lambda*pow(r, lambda-1)*((2 - Q*(lambda+1))*cos((lambda-1)*theta) - (lambda-1)*cos((lambda-3)*theta));
    //stress y
    stress(1,1) = A*lambda*pow(r, lambda-1)*((2 + Q*(lambda+1))*cos((lambda-1)*theta) + (lambda-1)*cos((lambda-3)*theta));
    //stress xy
    stress(0,1) = A*lambda*pow(r, lambda-1)*((lambda-1)*sin((lambda-3)*theta) + Q*(lambda+1)*sin((lambda-1)*theta));
    stress(1,0) = stress(0,1);

    VecDouble n(2);
    n[0] = -sqrt(2)/2;
    n[1] = -sqrt(2)/2;

    force = stress * n;
};

auto forcingFunctionFA = [](const VecDouble &coord, VecDouble &force){
    const auto &x=coord[0];
    const auto &y=coord[1];
    
    REAL A = 1.0;
    REAL Q = 0.543075579;
    REAL lambda = 0.544483737;
    REAL r = sqrt(x*x + y*y);
    REAL theta = atan2(y, x);

    MatrixDouble stress(2,2);
    //stress x
    stress(0,0) = A*lambda*pow(r, lambda-1)*((2 - Q*(lambda+1))*cos((lambda-1)*theta) - (lambda-1)*cos((lambda-3)*theta));
    //stress y
    stress(1,1) = A*lambda*pow(r, lambda-1)*((2 + Q*(lambda+1))*cos((lambda-1)*theta) + (lambda-1)*cos((lambda-3)*theta));
    //stress xy
    stress(0,1) = A*lambda*pow(r, lambda-1)*((lambda-1)*sin((lambda-3)*theta) + Q*(lambda+1)*sin((lambda-1)*theta));
    stress(1,0) = stress(0,1);

    VecDouble n(2);
    n[0] = sqrt(2)/2;
    n[1] = -sqrt(2)/2;

    force = stress * n;
};

auto exactSol = [](const VecDouble &coord, VecDouble &u, MatrixDouble &gradU){
    const auto &x=coord[0];
    const auto &y=coord[1];

    REAL E = 1.e7;

    u[0] = (1/E) * (24*x*y - 120*x);
    u[1] = (1/E) * (-12*x*x - 3.6*y*y + 36*y);

};

void CreateModel(CompMesh *cmesh);
// void CreateModel2(CompMesh *cmesh);
void SolveProblem(CompMesh *cmesh);

REAL ModElasticity = 1.;
REAL PoissonRatio = 0.3;

int main(int argc, char **args) { 

#ifdef HAS_PETSC
    // Starts main program invoking PETSc
    PetscInitialize(&argc, &args, (char*)0, (char*)0);
#endif

    //Create Model
    GeoMesh *gmesh = new GeoMesh();
    GmshTools::Read(*gmesh,"../chapaLRef1.msh");
    CompMesh *cmesh = new CompMesh(gmesh,ApproxType::EIsoparametric);
    CreateModel(cmesh);
    cmesh->Print("cmesh.txt");
    gmesh->Print("gmesh.txt");

    // CompMesh *cmeshaux = new CompMesh(gmesh,ApproxType::EIsoparametric);
    // CreateModel2(cmeshaux);
    // //gmesh->Print("gmesh.txt");
    // cmeshaux->Print("cmeshAux.txt");
    
    // VecDouble Solution;
    //Solve Problem
    SolveProblem(cmesh);

    //Compute strain energy
    //LinearAnalysis an(cmeshaux,SolverType::ELDLt);
    //an.Compute();   
    // an.PrintGlobalMatrix();

    // EigenSpMatrix *spMat = dynamic_cast<EigenSpMatrix *>(an.GlobalMatrix());
    // if (!spMat) {
    //     std::cerr << "Error: GlobalMatrix is not of type EigenSpMatrix." << std::endl;
    // }
    // auto globalMat = spMat->Matrix();

    // std::cout << "Global Matrix: \n" << globalMat << std::endl;

    // VecDouble Force = globalMat * Solution;
    // std::cout << "Force: \n" << Force << std::endl;
    // REAL strainEnergy = (Solution.dot(Force))/2;
    // std::cout << std::fixed << std::setprecision(10) << "Strain Energy: "<< strainEnergy << std::endl;
    

#ifdef HAS_PETSC
    //Finalize main program   
    PetscFinalize();
#endif


}   

void CreateModel(CompMesh *cmesh){

    Elasticity2D *matelasticity = new Elasticity2D(1, ModElasticity, PoissonRatio, false, 1.0); //domínio global
    //matelasticity->SetExactSolution(exactSol);

    cmesh->InsertMaterial(matelasticity);
    
    //BC 
    MatrixDouble val1(2,2);
    val1.setZero();
    VecDouble val2(2);
    val2.setZero();

    //Chapa L: apoio horizontal no pontos B e F, apoio vertical no ponto A
    val2[0] = 1.0;
    L2Projection * matbc1 = new L2Projection(2,dimension-2,BoundaryConditionType::kDirectionalHomogeneousDirichlet,val1,val2);
    val2.setZero();
    val2[1] = 1.0;
    L2Projection * matbc2 = new L2Projection(3,dimension-2,BoundaryConditionType::kDirectionalHomogeneousDirichlet,val1,val2);
    
    val2.setZero();
    L2Projection * matbc3 = new L2Projection(4,dimension-1,BoundaryConditionType::kNeumann,val1,val2);
    matbc3->SetForcingFunction(forcingFunctionAB);
    L2Projection * matbc4 = new L2Projection(5,dimension-1,BoundaryConditionType::kNeumann,val1,val2);
    matbc4->SetForcingFunction(forcingFunctionBC);
    L2Projection * matbc5 = new L2Projection(6,dimension-1,BoundaryConditionType::kNeumann,val1,val2);
    matbc5->SetForcingFunction(forcingFunctionEF);
    L2Projection * matbc6 = new L2Projection(7,dimension-1,BoundaryConditionType::kNeumann,val1,val2);
    matbc6->SetForcingFunction(forcingFunctionFA);

    //Chapa quadrada cisalhamento
    // val2[0] = 0.001;
    // L2Projection * matbc1 = new L2Projection(2,dimension-1,BoundaryConditionType::kDirectionalNonHomogeneousDirichlet,val1,val2);
    // val2.setZero();
    // val2[1] = 1.0;
    // L2Projection * matbc2 = new L2Projection(3,dimension-1,BoundaryConditionType::kDirectionalHomogeneousDirichlet,val1,val2);
    // val2.setZero();
    // //val2[0] = 1.0;
    // L2Projection * matbc3 = new L2Projection(4,dimension-1,BoundaryConditionType::kDirichlet,val1,val2);

    //Viga flexão pura
    //L2Projection * matbc1 = new L2Projection(2,dimension-2,BoundaryConditionType::kDirichlet,val1,val2);
    //val2[0] = 1.0;
    //L2Projection * matbc2 = new L2Projection(3,dimension-2,BoundaryConditionType::kDirectionalHomogeneousDirichlet,val1,val2);
    //L2Projection * matbc3 = new L2Projection(4,dimension-1,BoundaryConditionType::kNeumann,val1,val2);
    //matbc3->SetForcingFunction(forcingFunction1);
    //L2Projection * matbc4 = new L2Projection(5,dimension-1,BoundaryConditionType::kNeumann,val1,val2);
    //matbc4->SetForcingFunction(forcingFunction2);

    cmesh->InsertMaterial(matbc1);
    cmesh->InsertMaterial(matbc2);
    cmesh->InsertMaterial(matbc3);
    cmesh->InsertMaterial(matbc4);
    cmesh->InsertMaterial(matbc5);
    cmesh->InsertMaterial(matbc6);
    
    cmesh->AutoBuild();
}

/* void CreateModel2(CompMesh *cmesh){

    Elasticity2D *matelasticity = new Elasticity2D(1, ModElasticity, PoissonRatio, false, 1.0); //domínio global
    cmesh->InsertMaterial(matelasticity);
    
    cmesh->AutoBuild();
} */

void SolveProblem(CompMesh *cmesh){

    LinearAnalysis an(cmesh,SolverType::ELU);
    
    std::vector<std::string> ScalarNames, VectorNames;
    ScalarNames = {"SigmaX","SigmaY","TauXY","StrainEnergy"};
    VectorNames = {"Displacement","ExactDisplacement"}; //, "ExactDisplacement"}

    an.Run();
    // an.PrintGlobalMatrix();
    // an.PrintSolution();
    // an.PrintGlobalRhs();

    //VecDouble errors(4);
    //an.PostProcessError(errors);

#ifdef HAS_PETSC
    PETScMatrix *spMat = dynamic_cast<PETScMatrix *>(an.GlobalMatrix());
    if (!spMat) {
        std::cerr << "Error: GlobalMatrix is not of type PETSc." << std::endl;
        return;
    }

    Vec sol = spMat->Solution();
    Vec rhs = spMat->Rhs();   
    REAL strainEnergy;
    VecDot(sol, rhs, &strainEnergy);
    std::cout << std::fixed << std::setprecision(10) << "Strain Energy: "<< strainEnergy/2 << std::endl;



#else
    EigenSpMatrix *spMat = dynamic_cast<EigenSpMatrix *>(an.GlobalMatrix());
    if (!spMat) {
        std::cerr << "Error: GlobalMatrix is not of type EigenSpMatrix." << std::endl;
        return;
    }

    VecDouble sol = spMat->Solution();
    VecDouble rhs = spMat->Rhs();   

    REAL strainEnergy = (sol.dot(rhs))/2;
    std::cout << std::fixed << std::setprecision(30) << "Strain Energy: "<< strainEnergy << std::endl;


#endif
    
   
    VTUGenerator::PrintResults(cmesh,"mefResult",ScalarNames,VectorNames);
}
