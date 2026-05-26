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

void CreateGlobalModel(CompMesh *cmeshG);
void SolveGlobalProblem(CompMesh *cmeshG);
void CreateLocalModel(CompMesh *cmeshG, CompMesh *cmeshL);
void SolveLocalProblem(CompMesh *cmeshL);
void LocalToGlobalCorrespondence(CompMesh *cmeshG, CompMesh *cmeshL);
void LocalToGlobalCorrespondenceBoundary(CompMesh *cmeshG, CompMesh *cmeshL);
void CreateEnrichedModel(CompMesh *cmeshG, CompMesh *cmeshL);
void SolveEnrichedProblem(CompMesh *cmeshG);
bool CheckConvergence(int it);

int overlappingRegion;
int overlappingNHDirichletBoundary;
int overlappingNHNeumannBoundary;
int globalLocalIterations;
double globalLocalTolerance;
double ModElasticity = 1.;
double PoissonRatio = 0.49999;

VecDouble previousSol;
VecDouble currentSol;
//Local index to global index correspondence for elements
std::map<int64_t,int64_t> globalElementCorrespondence;
//Local node/integration point to global node/integration point correspondence 
std::map<int64_t, MatrixDouble> globalNodeCorrespondence;
//Global connects to be enriched to the index of the new connect associated with the enriched degree of freedom
std::map<int64_t,int64_t> enrichedConnects;

int main(int argc, char **args) { 

    overlappingRegion = 2;
    globalLocalIterations = 2;
    globalLocalTolerance = 1e-4;
    overlappingNHDirichletBoundary = 2;
    //overlappingNHNeumannBoundary = 3;

    //GLOBAL MODEL
    std::cout << "Solve Global Problem \n";

    //Create Global Model
    GeoMesh *gmeshG = new GeoMesh();
    GmshTools::Read(*gmeshG,"../chapaLGlobal.msh");
    CompMesh *cmeshG = new CompMesh(gmeshG,ApproxType::EIsoparametric);
    CreateGlobalModel(cmeshG);
    //gmeshG->Print("gmeshGlobal.txt");
    cmeshG->Print("cmeshGlobal.txt");
    
    //Solve Global Problem
    SolveGlobalProblem(cmeshG);

    CompMesh *cmeshL = nullptr;
    GeoMesh * gmeshL = new GeoMesh();

    GlobalLocalEnrichment *globalLocal = new GlobalLocalEnrichment(overlappingRegion, dimension, ModElasticity, PoissonRatio, false, 1.0);
    globalLocal->SetForcingFunction(forcingFunction);

    int it = 0;
    //for(int it = 0; it < globalLocalIterations; it++){
    while(it < globalLocalIterations){
        std::cout << "\nGlobal-Local Iteration " << it + 1 << ": ";

        //LOCAL MODEL
        std::cout << "\nSolve Local Problem \n";
        if(it == 0){
            //Create Local Model
            GmshTools::Read(*gmeshL,"../chapaLLocal.msh");
            cmeshL = new CompMesh(gmeshL,ApproxType::EIsoparametric);
            CreateLocalModel(cmeshG, cmeshL);
            gmeshL->Print("gmeshLocal.txt");
            //cmeshL->Print("cmeshLocal.txt");

            MEFGGlobalLocalTools::LocalToGlobalCorrespondence(cmeshG,cmeshL, globalElementCorrespondence, globalNodeCorrespondence, overlappingRegion, overlappingNHDirichletBoundary);
            MEFGGlobalLocalTools::LocalToGlobalCorrespondenceBoundary(cmeshG,cmeshL, globalElementCorrespondence, globalNodeCorrespondence, overlappingNHNeumannBoundary);
        }
        //Solve Local Problem
        SolveLocalProblem(cmeshL);

        //ENRICHED GLOBAL MODEL
        std::cout << "\nSolve Enriched Global Problem \n";
        if(it == 0) MEFGGlobalLocalTools::CreateEnrichedModel(cmeshG,cmeshL, globalElementCorrespondence, globalNodeCorrespondence, enrichedConnects, overlappingNHNeumannBoundary, globalLocal);

        cmeshG->Print("cmeshEnrichedGlobal.txt");
        SolveEnrichedProblem(cmeshG);

        if(CheckConvergence(it)) break;
        it++;
   };
}   

void CreateGlobalModel(CompMesh *cmeshG){

    Elasticity2D *matelasticityG1 = new Elasticity2D(1, ModElasticity, PoissonRatio, false, 1.0); //domínio global
    Elasticity2D *matelasticityG2 = new Elasticity2D(2, ModElasticity, PoissonRatio, false, 1.0); //domínio local

    //matelasticityG1->SetExactSolution(exactSol);

    cmeshG->InsertMaterial(matelasticityG1);
    cmeshG->InsertMaterial(matelasticityG2);

    //Chapa retangular tracionada
    //enrichedConnects[0]=-1;
    //enrichedConnects[5]=-1;

    //Chapa L 
    enrichedConnects[8]=-1;

    /*enrichedConnects[28]=-1;
    enrichedConnects[29]=-1;
    enrichedConnects[51]=-1;
    enrichedConnects[52]=-1;
    enrichedConnects[62]=-1;
    enrichedConnects[63]=-1;
    enrichedConnects[64]=-1;*/

    /*enrichedConnects[44]=-1;
    enrichedConnects[45]=-1;
    enrichedConnects[113]=-1;
    enrichedConnects[116]=-1;
    enrichedConnects[206]=-1;
    enrichedConnects[215]=-1;
    enrichedConnects[221]=-1;*/

    /*enrichedConnects[43]=-1;
    enrichedConnects[44]=-1;
    enrichedConnects[45]=-1;
    enrichedConnects[46]=-1;
    enrichedConnects[112]=-1;
    enrichedConnects[113]=-1;
    enrichedConnects[115]=-1;
    enrichedConnects[116]=-1;
    enrichedConnects[198]=-1;
    enrichedConnects[200]=-1;
    enrichedConnects[204]=-1;
    enrichedConnects[206]=-1;
    enrichedConnects[207]=-1;
    enrichedConnects[209]=-1;
    enrichedConnects[213]=-1;
    enrichedConnects[215]=-1;
    enrichedConnects[216]=-1;
    enrichedConnects[217]=-1;
    enrichedConnects[220]=-1;
    enrichedConnects[221]=-1;*/

    /*enrichedConnects[76]=-1;
    enrichedConnects[77]=-1;
    enrichedConnects[237]=-1;
    enrichedConnects[244]=-1;
    enrichedConnects[734]=-1;
    enrichedConnects[783]=-1;
    enrichedConnects[809]=-1;*/

    /*enrichedConnects[73]=-1;
    enrichedConnects[74]=-1;
    enrichedConnects[75]=-1;
    enrichedConnects[76]=-1;
    enrichedConnects[77]=-1;
    enrichedConnects[78]=-1;
    enrichedConnects[79]=-1;
    enrichedConnects[80]=-1;
    enrichedConnects[234]=-1;
    enrichedConnects[235]=-1;
    enrichedConnects[236]=-1;
    enrichedConnects[237]=-1;
    enrichedConnects[241]=-1;
    enrichedConnects[242]=-1;
    enrichedConnects[243]=-1;
    enrichedConnects[244]=-1;
    enrichedConnects[686]=-1;
    enrichedConnects[688]=-1;
    enrichedConnects[692]=-1;
    enrichedConnects[694]=-1;
    enrichedConnects[702]=-1;
    enrichedConnects[705]=-1;
    enrichedConnects[719]=-1;
    enrichedConnects[725]=-1;
    enrichedConnects[727]=-1;
    enrichedConnects[728]=-1;
    enrichedConnects[729]=-1;
    enrichedConnects[730]=-1;
    enrichedConnects[731]=-1;
    enrichedConnects[732]=-1;
    enrichedConnects[733]=-1;
    enrichedConnects[734]=-1;
    enrichedConnects[735]=-1;
    enrichedConnects[737]=-1;
    enrichedConnects[741]=-1;
    enrichedConnects[743]=-1;
    enrichedConnects[751]=-1;
    enrichedConnects[754]=-1;
    enrichedConnects[768]=-1;
    enrichedConnects[774]=-1;
    enrichedConnects[776]=-1;
    enrichedConnects[777]=-1;
    enrichedConnects[778]=-1;
    enrichedConnects[779]=-1;
    enrichedConnects[780]=-1;
    enrichedConnects[781]=-1;
    enrichedConnects[782]=-1;
    enrichedConnects[783]=-1;
    enrichedConnects[784]=-1;
    enrichedConnects[785]=-1;
    enrichedConnects[788]=-1;
    enrichedConnects[789]=-1;
    enrichedConnects[796]=-1;
    enrichedConnects[799]=-1;
    enrichedConnects[805]=-1;
    enrichedConnects[806]=-1;
    enrichedConnects[807]=-1;
    enrichedConnects[808]=-1;
    enrichedConnects[809]=-1;
    enrichedConnects[810]=-1;
    enrichedConnects[811]=-1;
    enrichedConnects[812]=-1;
    enrichedConnects[813]=-1;
    enrichedConnects[814]=-1;*/
    
    //BC 
    MatrixDouble val1(2,2);
    val1.setZero();
    VecDouble val2(2);
    val2.setZero();

    //Chapa L teste: apoio fixo no lado BC e carregamento unitário uniforme ao longe de FA
    // L2Projection * matbcG1 = new L2Projection(3,dimension-1,BoundaryConditionType::kDirichlet,val1,val2);
    // L2Projection * matbcG2 = new L2Projection(4,dimension-1,BoundaryConditionType::kNeumann,val1,val2);
    // matbcG2->SetForcingFunction(forcingFunction);

    //Chapa L: apoio horizontal no pontos B e F, apoio vertical no ponto A
    val2[0] = 1.0;
    L2Projection * matbcG1 = new L2Projection(3,dimension-2,BoundaryConditionType::kDirectionalHomogeneousDirichlet,val1,val2);
    val2.setZero();
    val2[1] = 1.0;
    L2Projection * matbcG2 = new L2Projection(4,dimension-2,BoundaryConditionType::kDirectionalHomogeneousDirichlet,val1,val2);
    
    //Chapa L: apoio fixo no ponto D, apoio vertical no ponto A
    //L2Projection * matbcG1 = new L2Projection(3,dimension-2,BoundaryConditionType::kDirichlet,val1,val2);
    //val2[0] = 1.0;
    //L2Projection * matbcG2 = new L2Projection(4,dimension-2,BoundaryConditionType::kDirectionalHomogeneousDirichlet,val1,val2);

    val2.setZero();
    L2Projection * matbcG3 = new L2Projection(5,dimension-1,BoundaryConditionType::kNeumann,val1,val2);
    matbcG3->SetForcingFunction(forcingFunctionAB);
    L2Projection * matbcG4 = new L2Projection(6,dimension-1,BoundaryConditionType::kNeumann,val1,val2);
    matbcG4->SetForcingFunction(forcingFunctionBC);
    L2Projection * matbcG5 = new L2Projection(7,dimension-1,BoundaryConditionType::kNeumann,val1,val2);
    matbcG5->SetForcingFunction(forcingFunctionEF);
    L2Projection * matbcG6 = new L2Projection(8,dimension-1,BoundaryConditionType::kNeumann,val1,val2);
    matbcG6->SetForcingFunction(forcingFunctionFA);

    cmeshG->InsertMaterial(matbcG1);
    cmeshG->InsertMaterial(matbcG2);
    cmeshG->InsertMaterial(matbcG3);
    cmeshG->InsertMaterial(matbcG4);
    cmeshG->InsertMaterial(matbcG5);
    cmeshG->InsertMaterial(matbcG6);
    
    cmeshG->AutoBuild();
}

void SolveGlobalProblem(CompMesh *cmeshG){

    LinearAnalysis anG(cmeshG,SolverType::ELDLt);
       
    std::vector<std::string> ScalarNames, VectorNames;
    ScalarNames = {"SigmaX","SigmaY","TauXY"};
    VectorNames = {"Displacement"}; //, "ExactDisplacement"}

    anG.Run();
    //anG.PrintGlobalMatrix();
    //anG.PrintSolution();
    //anG.PrintGlobalRhs();

    //VecDouble errors(4);
    //anG.PostProcessError(errors);

    EigenSpMatrix *spMat = dynamic_cast<EigenSpMatrix *>(anG.GlobalMatrix());
    if (!spMat) {
        std::cerr << "Error: GlobalMatrix is not of type EigenSpMatrix." << std::endl;
        return;
    }

    VecDouble sol = spMat->Solution();
    VecDouble rhs = spMat->Rhs();
    double strainEnergy = (sol.dot(rhs))/2;
    std::cout << std::fixed << std::setprecision(10) << "Strain Energy: "<< strainEnergy << std::endl;

    VTUGenerator::PrintResults(cmeshG,"globalResult",ScalarNames,VectorNames);
}

void CreateLocalModel(CompMesh *cmeshG, CompMesh * cmeshL){

    Elasticity2D* matelasticityL = new Elasticity2D(1, ModElasticity, PoissonRatio, false, 1.0);
    cmeshL->InsertMaterial(matelasticityL);

    //BC;
    MatrixDouble val1(2,2);
    val1.setZero();
    VecDouble val2(2);
    val2.setZero();
    
    InterpolatedBC * matbcL1 = new InterpolatedBC(2,dimension-1,2,BoundaryConditionType::kDirichlet,&globalElementCorrespondence,&globalNodeCorrespondence,cmeshG);
    //Chapa L: apoio fixo no ponto D, apoio vertical no ponto A
    //L2Projection * matbcL2 = new L2Projection(3,dimension-2,BoundaryConditionType::kDirichlet,val1,val2);
    
    //val2[1] = 1.0;
    //L2Projection * matbcL2 = new L2Projection(3,dimension-2,BoundaryConditionType::kDirichlet,val1,val2);
    //val2.setZero();
    //val2[0] = -2.;
    //L2Projection * matbcL2 = new L2Projection(2,dimension-1,BoundaryConditionType::kDirectionalNonHomogeneousDirichlet,val1,val2);
    //L2Projection * matbcL3 = new L2Projection(4,dimension-1,BoundaryConditionType::kNeumann,val1,val2);
    //matbcL3->SetForcingFunction(forcingFunction);

    cmeshL->InsertMaterial(matbcL1);
    //cmeshL->InsertMaterial(matbcL2);
    //cmeshL->InsertMaterial(matbcL3);
    
    cmeshL->AutoBuild();
}

void SolveLocalProblem(CompMesh *cmeshL){

    LinearAnalysis anL(cmeshL,SolverType::ELDLt);

    anL.Run();
    // anL.PrintGlobalMatrix();
    // anL.PrintGlobalRhs();
    // anL.PrintSolution();

    std::vector<std::string> ScalarNames, VectorNames;
    ScalarNames = {"SigmaX","SigmaY","TauXY"};
    VectorNames = {"Displacement"};

    EigenSpMatrix *spMat = dynamic_cast<EigenSpMatrix *>(anL.GlobalMatrix());
    if (!spMat) {
        std::cerr << "Error: GlobalMatrix is not of type EigenSpMatrix." << std::endl;
        return;
    }

    VecDouble sol = spMat->Solution();
    VecDouble rhs = spMat->Rhs();

    double bignumber = 1.e20;
    for (int64_t i = 0; i < rhs.size(); i++){
        if (rhs[i]>1.e10) rhs[i] = 0.;
        if (rhs[i]<-1.e10) rhs[i] = 0.;
    }

    double strainEnergy = (sol.dot(rhs))/2;
    
    

    std::cout << "Strain Energy: "<< strainEnergy << std::endl;

    VTUGenerator::PrintResults(cmeshL,"localResult",ScalarNames,VectorNames); 
};

void SolveEnrichedProblem(CompMesh *cmeshG){

    LinearAnalysis anE(cmeshG,SolverType::ELU);

    anE.Run();
    // anE.PrintGlobalMatrix();
    // anE.PrintGlobalRhs();
    // anE.PrintSolution();
    
    EigenSpMatrix *spMat = dynamic_cast<EigenSpMatrix *>(anE.GlobalMatrix());
    if (!spMat) {
        std::cerr << "Error: GlobalMatrix is not of type EigenSpMatrix." << std::endl;
        return;
    }

    currentSol = spMat->Solution();

    VecDouble rhs = spMat->Rhs();
    double strainEnergy = (currentSol.dot(rhs))/2;
    std::cout << "Strain Energy: "<< strainEnergy << std::endl;

    std::vector<std::string> ScalarNames, VectorNames;
    // ScalarNames = {"SigmaX","SigmaY","TauXY"};
    VectorNames = {"Displacement"};

    // VTUGenerator::PrintResults(cmeshG,"Enriched",ScalarNames,VectorNames); 
};

bool CheckConvergence(int it){

    if(it == 0){
        previousSol = currentSol;
        return false;
    }else{
        MatrixDouble diff = currentSol - previousSol;

        if(diff.norm() <= globalLocalTolerance){
            return true;
        }else{
            previousSol = currentSol;
            return false;
        }
    }
};