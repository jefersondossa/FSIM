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

auto forcingFunction = [](const VecDouble &coord, VecDouble &force){
    const auto &x=coord[0];
    const auto &y=coord[1];
    // force[0] = -0.5+1.0*y;
    force[0] = -0.5+1.0*y;
};

auto forcingFunctionAB = [](const VecDouble &coord, VecDouble &force){
    const auto &x=coord[0];
    const auto &y=coord[1];
    
    double A = 1.0;
    double Q = 0.543075579;
    double lambda = 0.544483737;
    double r = sqrt(x*x + y*y);
    double theta = atan2(y, x);

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
    
    double A = 1.0;
    double Q = 0.543075579;
    double lambda = 0.544483737;
    double r = sqrt(x*x + y*y);
    double theta = atan2(y, x);

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
    
    double A = 1.0;
    double Q = 0.543075579;
    double lambda = 0.544483737;
    double r = sqrt(x*x + y*y);
    double theta = atan2(y, x);

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
    
    double A = 1.0;
    double Q = 0.543075579;
    double lambda = 0.544483737;
    double r = sqrt(x*x + y*y);
    double theta = atan2(y, x);

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

void SolveProblem(CompMesh *cmesh);

int main(int argc, char **args) {   
    //Geometric Mesh
    GeoMesh * gmesh = new GeoMesh();
    GmshTools::Read(*gmesh,"../chapaLRef1.msh");
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
    std::cout << "Strain Energy: "<< strainEnergy << std::endl;

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
    VecDouble val3(nstate);

    val2.setZero();
    val3.setZero();
    
    val2[0] = 1.0;
    L2Projection * matbc1 = new L2Projection(2,dimension-2,BoundaryConditionType::kDirectionalHomogeneousDirichlet,val1,val2);
    val2.setZero();
    val2[1] = 1.0;
    L2Projection * matbc2 = new L2Projection(3,dimension-2,BoundaryConditionType::kDirectionalHomogeneousDirichlet,val1,val2);
    
    L2Projection * matbc3 = new L2Projection(4,dimension-1,BoundaryConditionType::kNeumann,val1,val3);
    matbc3->SetForcingFunction(forcingFunctionAB);
    L2Projection * matbc4 = new L2Projection(5,dimension-1,BoundaryConditionType::kNeumann,val1,val3);
    matbc4->SetForcingFunction(forcingFunctionBC);
    L2Projection * matbc5 = new L2Projection(6,dimension-1,BoundaryConditionType::kNeumann,val1,val3);
    matbc5->SetForcingFunction(forcingFunctionEF);
    L2Projection * matbc6 = new L2Projection(7,dimension-1,BoundaryConditionType::kNeumann,val1,val3);
    matbc6->SetForcingFunction(forcingFunctionFA);

    cmesh->InsertMaterial(matbc1);
    cmesh->InsertMaterial(matbc2);
    cmesh->InsertMaterial(matbc3);
    cmesh->InsertMaterial(matbc4);
    cmesh->InsertMaterial(matbc5);
    cmesh->InsertMaterial(matbc6);

    cmesh->AutoBuild(); 

    return cmesh;
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
    L2Projection * matbc1 = new L2Projection(4,1,BoundaryConditionType::kNeumann,val1,val2);
    L2Projection * matbc2 = new L2Projection(5,1,BoundaryConditionType::kNeumann,val1,val2);
    L2Projection * matbc3 = new L2Projection(6,1,BoundaryConditionType::kNeumann,val1,val2);
    L2Projection * matbc4 = new L2Projection(7,1,BoundaryConditionType::kNeumann,val1,val2);

    cmesh->InsertMaterial(matbc1);
    cmesh->InsertMaterial(matbc2);
    cmesh->InsertMaterial(matbc3);
    cmesh->InsertMaterial(matbc4);

    cmesh->AutoBuild();
    
    cmesh->Print("cmesh_pressure.txt");

    return cmesh;
}

MixedCompMesh* CreateMixedMesh(std::vector<CompMesh *> &meshvector){

    MixedCompMesh* cmesh = new MixedCompMesh(meshvector);

    MixedElasticity * mat = new MixedElasticity(1, 2, 1.0, 0.);
    cmesh->InsertMaterial(mat);

    cmesh->AutoBuild();

    cmesh->Print("cmesh_mixed.txt");

    return cmesh;
};
