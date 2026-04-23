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

    double A = 1.0;
    double Q = 0.543075579;
    double lambda = 0.544483737;
    double r = sqrt(x*x + y*y);
    double theta = atan2(y, x);

    double young = 1.0;
    double nu = 0.3;
    double kappa = 3 - 4*nu;
    double G = young / (2*(1+nu));

    u[0] = (A/(2*G)) * pow(r, lambda) * ((kappa - Q * (lambda+1))*cos(lambda*theta) - lambda*cos((lambda-2)*theta));
    u[1] = (A/(2*G)) * pow(r, lambda) * ((kappa + Q * (lambda+1))*sin(lambda*theta) + lambda*sin((lambda-2)*theta));
    
    double dUxdr = (A/(2*G)) * lambda * pow(r, (lambda-1)) * ((kappa - Q*(lambda+1))*cos(lambda*theta) - lambda*cos((lambda-2)*theta));
    double dUydr = (A/(2*G)) * lambda * pow(r, (lambda-1)) * ((kappa + Q*(lambda+1))*cos(lambda*theta) + lambda*cos((lambda-2)*theta));

    double dUxdtheta = (A/(2*G)) * pow(r, lambda) * (-(kappa - Q*(lambda+1)) * lambda * sin(lambda*theta) + lambda * (lambda-2) * sin((lambda-2)*theta));
    double dUydtheta = (A/(2*G)) * pow(r, lambda) * ((kappa + Q*(lambda+1)) * lambda * cos(lambda*theta) + lambda * (lambda-2) * cos((lambda-2)*theta));

    gradU(0,0) = dUxdr*(x/r) + dUxdtheta*(-y/pow(r, 2));
    gradU(0,1) = dUxdr*(y/r) + dUxdtheta*(x/pow(r, 2));
    gradU(1,0) = dUydr*(x/r) + dUydtheta*(-y/pow(r, 2));
    gradU(1,1) = dUydr*(y/r) + dUxdtheta*(x/pow(r, 2));
};

void CreateModel(CompMesh *cmesh);
void SolveProblem(CompMesh *cmesh);

double ModElasticity = 1.;
double PoissonRatio = 0.3;

int main(int argc, char **args) { 

    //Create Model
    GeoMesh *gmesh = new GeoMesh();
    GmshTools::Read(*gmesh,"../chapaQuadrada.msh");
    CompMesh *cmesh = new CompMesh(gmesh,ApproxType::EIsoparametric);
    CreateModel(cmesh);
    //gmesh->Print("gmesh.txt");
    cmesh->Print("cmesh.txt");
    
    //Solve Problem
    SolveProblem(cmesh);
    
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

void SolveProblem(CompMesh *cmesh){

    LinearAnalysis an(cmesh,SolverType::ELDLt);
       
    std::vector<std::string> ScalarNames, VectorNames;
    ScalarNames = {"SigmaX","SigmaY","TauXY"};
    VectorNames = {"Displacement"}; //, "ExactDisplacement"}

    an.Run();
    //anG.PrintGlobalMatrix();
    //anG.PrintSolution();
    //anG.PrintGlobalRhs();

    //VecDouble errors(4);
    //anG.PostProcessError(errors);

    EigenSpMatrix *spMat = dynamic_cast<EigenSpMatrix *>(an.GlobalMatrix());
    if (!spMat) {
        std::cerr << "Error: GlobalMatrix is not of type EigenSpMatrix." << std::endl;
        return;
    }

    VecDouble sol = spMat->Solution();
    VecDouble rhs = spMat->Rhs();
    double strainEnergy = (sol.dot(rhs))/2;
    std::cout << std::fixed << std::setprecision(10) << "Strain Energy: "<< strainEnergy << std::endl;

    VTUGenerator::PrintResults(cmesh,"mefResult",ScalarNames,VectorNames);
}
