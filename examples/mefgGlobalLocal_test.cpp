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
    force[0] = -1*y;
    // force[0] = -1;
};

void CreateGlobalModel(CompMesh *cmeshG);
void SolveGlobalProblem(CompMesh *cmeshG);
void CreateLocalModel(CompMesh *cmeshG, CompMesh *cmeshL);
void SolveLocalProblem(CompMesh *cmeshL);
void CreateEnrichedModel(CompMesh *cmeshG, CompMesh *cmeshL);
void SolveEnrichedProblem(CompMesh *cmeshG);

int overlappingRegion;
int overlappingNHDirichletBoundary;
int overlappingNHNeumannBoundary;
int globalLocalIterations;
REAL ModElasticity = 1.;
REAL PoissonRatio = 0.0;

//Local index to global index correspondence for elements
std::map<int,int> globalElementCorrespondence;
//Local node/integration point to global node/integration point correspondence 
std::map<int, MatrixDouble> globalNodeCorrespondence;
//Global connects to be enriched to the index of the new connect associated with the enriched degree of freedom
std::map<int,int> enrichedConnects;

int main(int argc, char **args) { 

    overlappingRegion = 2;
    globalLocalIterations = 1;
    overlappingNHDirichletBoundary = 2;
    overlappingNHNeumannBoundary = 4;

    //GLOBAL MODEL
    std::cout << "Solve Global Problem \n";

    //Create Global Model
    GeoMesh *gmeshG = new GeoMesh();
    GmshTools::Read(*gmeshG,"../global.msh");
    CompMesh *cmeshG = new CompMesh(gmeshG,ApproxType::EHierarquic);
    cmeshG->SetDefaultOrder(1);
    CreateGlobalModel(cmeshG);
    gmeshG->Print("gmeshGlobal.txt");
    cmeshG->Print("cmeshGlobal.txt");
    
    //Solve Global Problem
    SolveGlobalProblem(cmeshG);
    
    //Automatic search the connects to be enriched in the global model
    for (auto celEl:cmeshG->ElementVec()){
        auto geoEl = celEl->Reference();
        //Create the enriched element associated with the global element 
        //Now put the rule to find the connects.
        // Ex: if all element node coordinate x < .1, then they will be enriched. 
        bool enriched = true;
        auto geoNodes = geoEl->getGeometricNodes();
        for (auto node:geoNodes){
            if (gmeshG->NodeVec()[node]->getCoordinateValue(0) > 0.1){
                enriched = false;
            }
        }
        if (enriched){
            auto connects = celEl->getConnectivity();
            for (int i = 0; i < connects.size(); i++){
                enrichedConnects[connects[i]->Index()]=-1;
            }
        }
    }    

    CompMesh *cmeshL = nullptr;
    GeoMesh * gmeshL = new GeoMesh();

    //Create Local Model
    GmshTools::Read(*gmeshL,"../local.msh");
    cmeshL = new CompMesh(gmeshL,ApproxType::EHierarquic);
    cmeshL->SetDefaultOrder(1);
    CreateLocalModel(cmeshG, cmeshL);
    gmeshL->Print("gmeshLocal.txt");
    cmeshL->Print("cmeshLocal.txt");

    MEFGGlobalLocalTools::LocalToGlobalCorrespondence(cmeshG,cmeshL, globalElementCorrespondence, globalNodeCorrespondence, overlappingRegion, overlappingNHDirichletBoundary);
    MEFGGlobalLocalTools::LocalToGlobalCorrespondenceBoundary(cmeshG,cmeshL, globalElementCorrespondence, globalNodeCorrespondence, overlappingNHNeumannBoundary);

    // //Print correspondence for debugging
    // std::cout << "Global to local element correspondence: \n";
    // for (auto &elem:globalElementCorrespondence){
    //     std::cout << "Global element index: " << elem.second << " \nLocal element index: \n" << elem.first << "\n";
    // }
    // std::cout << "Global to local node correspondence: \n";
    // for (auto &node:globalNodeCorrespondence){
    //     std::cout << "Global node index: " << node.first << " \nLocal node indexes: \n" << node.second.transpose() << "\n";
    // }
    // std::cout << "Global connects to be enriched: \n";
    // for (auto &connect:enrichedConnects){
    //     std::cout << "Global connect index: " << connect.first << "\n";
    // }

    GlobalLocalEnrichment *globalLocal = new GlobalLocalEnrichment(1,dimension,ModElasticity,PoissonRatio);
    globalLocal->SetForcingFunction(forcingFunction);

    for(int n = 0; n < globalLocalIterations; n++){

        std::cout << "\nGlobal-Local Iteration " << n + 1 << ": ";

        //Solve Local Problem
        SolveLocalProblem(cmeshL);

        //ENRICHED GLOBAL MODEL
        std::cout << "\nSolve Enriched Global Problem \n";
        if(n == 0) MEFGGlobalLocalTools::CreateEnrichedModel(cmeshG,cmeshL, globalElementCorrespondence, globalNodeCorrespondence, enrichedConnects, overlappingNHNeumannBoundary, globalLocal);

        cmeshG->Print("cmeshEnrichedGlobal.txt");
        SolveEnrichedProblem(cmeshG);
   }
}   

void CreateGlobalModel(CompMesh *cmeshG){

    Elasticity2D *matelasticityG1 = new Elasticity2D(1, ModElasticity, PoissonRatio); //região de sobreposição do domínio local no domínio global
    Elasticity2D *matelasticityG2 = new Elasticity2D(overlappingRegion, ModElasticity, PoissonRatio); //domínio global complementar
    cmeshG->InsertMaterial(matelasticityG1);
    cmeshG->InsertMaterial(matelasticityG2);


    //BC 
    MatrixDouble val1(2,2);
    val1.setZero();
    VecDouble val2(2),val3(2);
    val2.setZero();
    val3.setZero();
    L2Projection * matbcG1 = new L2Projection(3,dimension-1,BoundaryConditionType::kDirichlet,val1,val2);
    val3[1] = 1.0;
    L2Projection * matbcG2 = new L2Projection(4,dimension-1,BoundaryConditionType::kDirectionalHomogeneousDirichlet,val1,val3);
    L2Projection * matbcG3 = new L2Projection(5,dimension-1,BoundaryConditionType::kNeumann,val1,val2);
    matbcG3->SetForcingFunction(forcingFunction);

    cmeshG->InsertMaterial(matbcG1);
    cmeshG->InsertMaterial(matbcG2);
    cmeshG->InsertMaterial(matbcG3);
    
    cmeshG->AutoBuild();
}

void SolveGlobalProblem(CompMesh *cmeshG){

    LinearAnalysis anG(cmeshG,SolverType::ELDLt);
       
    std::vector<std::string> ScalarNames, VectorNames;
    ScalarNames = {"SigmaX","SigmaY","TauXY"};
    VectorNames = {"Displacement"};

    anG.Run();
    // anG.PrintSolution();
    anG.PrintGlobalRhs();

    EigenSpMatrix *spMat = dynamic_cast<EigenSpMatrix *>(anG.GlobalMatrix());
    if (!spMat) {
        std::cerr << "Error: GlobalMatrix is not of type EigenSpMatrix." << std::endl;
        return;
    }

    VecDouble sol = spMat->Solution();
    VecDouble rhs = spMat->Rhs();
    REAL strainEnergy = (sol.dot(rhs))/2;
    std::cout << "Strain Energy: "<< strainEnergy << std::endl;

    //anG.PrintGlobalRhs();

    // VTUGenerator::PrintResults(cmeshG,"globalResult",ScalarNames,VectorNames);
}

void CreateLocalModel(CompMesh *cmeshG, CompMesh * cmeshL){

    Elasticity2D* matelasticityL = new Elasticity2D(1, ModElasticity, PoissonRatio);
    cmeshL->InsertMaterial(matelasticityL);

    //BC;
    MatrixDouble val1(2,2);
    val1.setZero();
    VecDouble val2(2);
    val2[1] = 1.0;
    L2Projection * matbcL1 = new L2Projection(3,dimension-1,BoundaryConditionType::kDirectionalHomogeneousDirichlet,val1,val2);
    val2.setZero();
    InterpolatedBC * matbcL2 = new InterpolatedBC(2,dimension-1,2,BoundaryConditionType::kDirichlet,&globalElementCorrespondence,&globalNodeCorrespondence,cmeshG);
    //val2[0] = -2.;
    //L2Projection * matbcL2 = new L2Projection(2,dimension-1,BoundaryConditionType::kDirectionalNonHomogeneousDirichlet,val1,val2);
    L2Projection * matbcL3 = new L2Projection(4,dimension-1,BoundaryConditionType::kNeumann,val1,val2);
    matbcL3->SetForcingFunction(forcingFunction);

    cmeshL->InsertMaterial(matbcL1);
    cmeshL->InsertMaterial(matbcL2);
    cmeshL->InsertMaterial(matbcL3);
    
    cmeshL->AutoBuild();
}

void SolveLocalProblem(CompMesh *cmeshL){

    LinearAnalysis anL(cmeshL,SolverType::ELDLt);

    anL.Run();
    //anL.PrintGlobalMatrix();
    //anL.PrintGlobalRhs();
    // anL.PrintSolution();

    std::vector<std::string> ScalarNames, VectorNames;
    ScalarNames = {"SigmaX","SigmaY","TauXY"};
    VectorNames = {"Displacement"};

    // VTUGenerator::PrintResults(cmeshL,"localResult",ScalarNames,VectorNames); 
};

void SolveEnrichedProblem(CompMesh *cmeshG){

    LinearAnalysis anE(cmeshG,SolverType::ELU);

    anE.Run();
    //anE.PrintGlobalMatrix();
    anE.PrintGlobalRhs();
    // anE.PrintSolution();

    EigenSpMatrix *spMat = dynamic_cast<EigenSpMatrix *>(anE.GlobalMatrix());
    if (!spMat) {
        std::cerr << "Error: GlobalMatrix is not of type EigenSpMatrix." << std::endl;
        return;
    }

    VecDouble currentSol = spMat->Solution();

    VecDouble rhs = spMat->Rhs();
    REAL strainEnergy = (currentSol.dot(rhs))/2;
    std::cout << "Strain Energy: "<< strainEnergy << std::endl;



    std::vector<std::string> ScalarNames, VectorNames;
    // ScalarNames = {"SigmaX","SigmaY","TauXY"};
    VectorNames = {"Displacement"};

    // VTUGenerator::PrintResults(cmeshG,"Enriched",ScalarNames,VectorNames); 
};