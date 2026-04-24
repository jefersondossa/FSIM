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
#include "MixedElasticity.h"
#include "ElementEnriched.h"
#include "Element.h"
#include "GlobalLocalEnrichment.h"
#include "MEFGGlobalLocalTools.h"
#include "MixedGlobalLocalEnrichment.h"

//Comentar erro ao debugar após incluir a weakform GlobalLocalEnrichment e o element ElementEnriched. 

// Defines the problem dimension
const int dimension = 2;
auto forcingFunction = [](const VecDouble &coord, VecDouble &force){
    const auto &x=coord[0];
    const auto &y=coord[1];
    // force[0] = -1*y;
    force[0] = -1.;
};

MixedCompMesh* CreateGlobalModel(GeoMesh *gmesh, std::vector<CompMesh *> &meshvector);
MixedCompMesh* CreateLocalModel(GeoMesh *gmesh, std::vector<CompMesh *> &meshvector, CompMesh *cmeshG);
void SolveGlobalProblem(CompMesh *cmeshG);
void SolveLocalProblem(CompMesh *cmeshL);
void SolveEnrichedProblem(CompMesh *cmeshG);

int overlappingRegion;
int overlappingNHDirichletBoundary;
int overlappingNHNeumannBoundary;
int globalLocalIterations;
double ModElasticity = 1.;
double PoissonRatio = 0.0;

//Local index to global index correspondence for elements
std::map<int64_t,int64_t> globalElementCorrespondence;
//Local node/integration point to global node/integration point correspondence 
std::map<int64_t, MatrixDouble> globalNodeCorrespondence;
//Global connects to be enriched to the index of the new connect associated with the enriched degree of freedom
std::map<int64_t,int64_t> enrichedConnects;

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
    CompMesh *cmeshDispG, *cmeshPressG;
    std::vector<CompMesh *> meshvectorG = {cmeshDispG, cmeshPressG};
    MixedCompMesh *cmeshG = CreateGlobalModel(gmeshG,meshvectorG);

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
    
    GeoMesh * gmeshL = new GeoMesh();
    GmshTools::Read(*gmeshL,"../local.msh");
    CompMesh *cmeshDispL, *cmeshPressL;
    std::vector<CompMesh *> meshvectorL = {cmeshDispL, cmeshPressL};
    MixedCompMesh *cmeshL = CreateLocalModel(gmeshL, meshvectorL,cmeshG);
    gmeshL->Print("gmeshLocal.txt");
    cmeshL->Print("cmeshLocal.txt");

    MEFGGlobalLocalTools::LocalToGlobalCorrespondence(cmeshG,cmeshL, globalElementCorrespondence, globalNodeCorrespondence, overlappingRegion, overlappingNHDirichletBoundary);
    MEFGGlobalLocalTools::LocalToGlobalCorrespondenceBoundary(cmeshG,cmeshL, globalElementCorrespondence, globalNodeCorrespondence, overlappingNHNeumannBoundary);

    //  //Print correspondence for debugging
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

    MixedGlobalLocalEnrichment *globalLocal = new MixedGlobalLocalEnrichment(1,dimension,ModElasticity,PoissonRatio);
    globalLocal->SetForcingFunction(forcingFunction);

    for(int n = 0; n < globalLocalIterations; n++){

        std::cout << "\nGlobal-Local Iteration " << n + 1 << ": ";

        //LOCAL MODEL
        std::cout << "\nSolve Local Problem \n";
        //Solve Local Problem
        SolveLocalProblem(cmeshL);

        //ENRICHED GLOBAL MODEL
        std::cout << "\nSolve Enriched Global Problem \n";
        if(n == 0) MEFGGlobalLocalTools::CreateMixedEnrichedModel(cmeshG,cmeshL, globalElementCorrespondence, globalNodeCorrespondence, enrichedConnects, overlappingNHNeumannBoundary, globalLocal);

        cmeshG->Print("cmeshEnrichedGlobal.txt");
        SolveEnrichedProblem(cmeshG);
        cmeshG->Print("cmeshEnrichedGlobal2.txt");
   }
}   

MixedCompMesh* CreateGlobalModel(GeoMesh *gmesh, std::vector<CompMesh *> &meshvector){

    //Displacement Cmesh
    meshvector[0] = new CompMesh(gmesh, ApproxType::EIsoparametric);
    //meshvector[0]->SetDefaultOrder(2);
    int nstate = 2;
    WeakForm *matelasticityG1 = new WeakForm(1, nstate); //região de sobreposição do domínio local no domínio global
    WeakForm *matelasticityG2 = new WeakForm(2, nstate); //domínio global complementar
    meshvector[0]->InsertMaterial(matelasticityG1);
    meshvector[0]->InsertMaterial(matelasticityG2);
    
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

    meshvector[0]->InsertMaterial(matbcG1);
    meshvector[0]->InsertMaterial(matbcG2);
    meshvector[0]->InsertMaterial(matbcG3);
    
    meshvector[0]->AutoBuild();

    //Pressure Cmesh
    meshvector[1] = new CompMesh(gmesh, ApproxType::EHierarquic);
    meshvector[1]->SetDefaultOrder(0);

    WeakForm* matelasticityG3 = new WeakForm(1, 1);
    WeakForm* matelasticityG4 = new WeakForm(2, 1);

    meshvector[1]->InsertMaterial(matelasticityG3);
    meshvector[1]->InsertMaterial(matelasticityG4);
    
    //BC - All Neumann for pressure
    MatrixDouble val4(1,1);
    VecDouble val5(1);
    val4.setZero();
    val5.setZero();
    L2Projection * matbc1 = new L2Projection(3,1,BoundaryConditionType::kNeumann,val4,val5);
    L2Projection * matbc2 = new L2Projection(4,1,BoundaryConditionType::kNeumann,val4,val5);
    L2Projection * matbc3 = new L2Projection(5,1,BoundaryConditionType::kNeumann,val4,val5);

    meshvector[1]->InsertMaterial(matbc1);
    meshvector[1]->InsertMaterial(matbc2);
    meshvector[1]->InsertMaterial(matbc3);

    meshvector[1]->CreateDisconnectedElements();
    meshvector[1]->AutoBuild();

    //MixedCompMesh

    MixedCompMesh* cmesh = new MixedCompMesh(meshvector);

    MixedElasticity * mat1 = new MixedElasticity(1,2,1.0,0.);
    MixedElasticity * mat2 = new MixedElasticity(2,2,1.0,0.);
    cmesh->InsertMaterial(mat1);
    cmesh->InsertMaterial(mat2);

    cmesh->AutoBuild();

    return cmesh;
}

void SolveGlobalProblem(CompMesh *cmeshG){

    LinearAnalysis anG(cmeshG,SolverType::ELDLt);
       
    std::vector<std::string> ScalarNames, VectorNames;
    ScalarNames = {"SigmaX","SigmaY","TauXY"};
    VectorNames = {"Displacement"};

    anG.Run();
    // anG.PrintGlobalMatrix();
    // anG.PrintSolution();
    // anG.PrintGlobalRhs();

    EigenSpMatrix *spMat = dynamic_cast<EigenSpMatrix *>(anG.GlobalMatrix());
    if (!spMat) {
        std::cerr << "Error: GlobalMatrix is not of type EigenSpMatrix." << std::endl;
        return;
    }

    VecDouble sol = spMat->Solution();
    VecDouble rhs = spMat->Rhs();
    double strainEnergy = (sol.dot(rhs))/2;
    std::cout << "Strain Energy: "<< strainEnergy << std::endl;

    VTUGenerator::PrintResults(cmeshG,"globalResult",ScalarNames,VectorNames);
}

MixedCompMesh* CreateLocalModel(GeoMesh *gmesh, std::vector<CompMesh *> &meshvector, CompMesh *cmeshG){
    //Displacement Cmesh
    meshvector[0] = new CompMesh(gmesh, ApproxType::EIsoparametric);
    //meshvector[0]->SetDefaultOrder(2);

    int nstate = 2;
    WeakForm* matelasticityL = new WeakForm(1, nstate);
    meshvector[0]->InsertMaterial(matelasticityL);

    //BC;
    MatrixDouble val1(2,2);
    val1.setZero();
    VecDouble val2(2);
    // val2[1] = 1.0;
    L2Projection * matbcL1 = new L2Projection(3,dimension-1,BoundaryConditionType::kNeumann,val1,val2);
    val2.setZero();
    InterpolatedBC * matbcL2 = new InterpolatedBC(2,dimension-1,2,BoundaryConditionType::kDirichlet,&globalElementCorrespondence,&globalNodeCorrespondence,cmeshG);
    // val2[0] = -2.;
    // L2Projection * matbcL2 = new L2Projection(2,dimension-1,BoundaryConditionType::kDirichlet,val1,val2);
    L2Projection * matbcL3 = new L2Projection(4,dimension-1,BoundaryConditionType::kNeumann,val1,val2);
    matbcL3->SetForcingFunction(forcingFunction);

    meshvector[0]->InsertMaterial(matbcL1);
    meshvector[0]->InsertMaterial(matbcL2);
    meshvector[0]->InsertMaterial(matbcL3);
    
    meshvector[0]->AutoBuild();

    //Pressure Cmesh
    meshvector[1] = new CompMesh(gmesh, ApproxType::EHierarquic);
    meshvector[1]->SetDefaultOrder(0);

    WeakForm* matelasticityL1 = new WeakForm(1, 1);
    meshvector[1]->InsertMaterial(matelasticityL1);
    
    //BC - All Neumann for pressure
    MatrixDouble val4(1,1);
    VecDouble val5(1);
    val4.setZero();
    val5.setZero();
    L2Projection * matbc1 = new L2Projection(2,1,BoundaryConditionType::kNeumann,val4,val5);
    L2Projection * matbc2 = new L2Projection(3,1,BoundaryConditionType::kNeumann,val4,val5);
    L2Projection * matbc3 = new L2Projection(4,1,BoundaryConditionType::kNeumann,val4,val5);

    meshvector[1]->InsertMaterial(matbc1);
    meshvector[1]->InsertMaterial(matbc2);
    meshvector[1]->InsertMaterial(matbc3);

    meshvector[1]->AutoBuild();

    //MixedCompMesh

    MixedCompMesh* cmesh = new MixedCompMesh(meshvector);

    MixedElasticity * mat1 = new MixedElasticity(1,2,1.0,0.);
    cmesh->InsertMaterial(mat1);

    cmesh->AutoBuild();

    return cmesh;


}

void SolveLocalProblem(CompMesh *cmeshL){

    LinearAnalysis anL(cmeshL,SolverType::ELDLt);

    anL.Run();
    // anL.PrintGlobalMatrix();
    // anL.PrintGlobalRhs();
    anL.PrintSolution();

    EigenSpMatrix *spMat = dynamic_cast<EigenSpMatrix *>(anL.GlobalMatrix());
    if (!spMat) {
        std::cerr << "Error: GlobalMatrix is not of type EigenSpMatrix." << std::endl;
        return;
    }
    MatrixDouble previousSol = spMat->Solution();

    anL.GlobalMatrix()->SolutionNorm();

    std::vector<std::string> ScalarNames, VectorNames;
    ScalarNames = {"SigmaX","SigmaY","TauXY"};
    VectorNames = {"Displacement"};

    VTUGenerator::PrintResults(cmeshL,"localResult",ScalarNames,VectorNames); 
};

void SolveEnrichedProblem(CompMesh *cmeshG){

    LinearAnalysis anE(cmeshG,SolverType::ELU);

    anE.Run();
    anE.PrintGlobalMatrix();
    anE.PrintGlobalRhs();
    anE.PrintSolution();

    EigenSpMatrix *spMat = dynamic_cast<EigenSpMatrix *>(anE.GlobalMatrix());
    if (!spMat) {
        std::cerr << "Error: GlobalMatrix is not of type EigenSpMatrix." << std::endl;
        return;
    }

    VecDouble currentSol = spMat->Solution();

    VecDouble rhs = spMat->Rhs();
    double strainEnergy = (currentSol.dot(rhs))/2;
    std::cout << "Strain Energy: "<< strainEnergy << std::endl;


    std::vector<std::string> ScalarNames, VectorNames;
    // ScalarNames = {"SigmaX","SigmaY","TauXY"};
    VectorNames = {"Displacement"};

    // VTUGenerator::PrintResults(cmeshG,"Enriched",ScalarNames,VectorNames); 
};