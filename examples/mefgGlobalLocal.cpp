#include "GeoMesh.h"
#include "LinearAnalysis.h"
#include "MixedCompMesh.h"
#include "GmshTools.h"
#include "VTUGenerator.h"
#include "L2Projection.h"
#include "Elasticity2D.h"
#include "CompMeshTools.h"

#include "ElementEnriched.h"
#include "Element.h"
#include "GlobalLocalEnrichment.h"

//Comentar erro ao debugar após incluir a weakform GlobalLocalEnrichment e o element ElementEnriched. 

// Defines the problem dimension
const int dimension = 2;
auto forcingFunction = [](const VecDouble &coord, VecDouble &force){
    const auto &x=coord[0];
    const auto &y=coord[1];
    force[0] = -1*y;
};

void CreateGlobalModel(CompMesh *cmeshG);
void SolveGlobalProblem(CompMesh *cmeshG);
void CreateLocalModel(CompMesh *cmeshL);
void SolveLocalProblem(CompMesh *cmeshL);
void LocalToGlobalCorrespondence(CompMesh *cmeshG, CompMesh *cmeshL);

int overlappingRegion;
int overlappingBoundary;
std::map<int64_t,int64_t> globalElementCorrespondence;
std::map<int64_t, MatrixDouble> globalNodeCorrespondence;
std::set<int64_t> enrichedNodes;

int main(int argc, char **args) { 

    //GLOBAL MODEL
    std::cout << "Solve Global Problem \n";

    //Create Global Model
    GeoMesh *gmeshG = new GeoMesh();
    GmshTools::Read(*gmeshG,"../global.msh");
    CompMesh *cmeshG = new CompMesh(gmeshG,ApproxType::EHierarquic);
    CreateGlobalModel(cmeshG);

    //Solve Global Problem
    SolveGlobalProblem(cmeshG);

    //LOCAL MODEL
    std::cout << "\nSolve Local Problem \n";

    //Create Local Model
    GeoMesh * gmeshL = new GeoMesh();
    GmshTools::Read(*gmeshL,"../local.msh");
    CompMesh *cmeshL = new CompMesh(gmeshL,ApproxType::EHierarquic);
    CreateLocalModel(cmeshL);

    //Solve Local Problem
    SolveLocalProblem(cmeshL);

    LocalToGlobalCorrespondence(cmeshG,cmeshL);

    Element *localEl = cmeshL->ElementVec()[6];
    Element *globalEl = cmeshG->ElementVec()[6];
    cmeshG->Print("cmeshG.txt");
    MatrixDouble stiffness(16,16);
    stiffness.setZero();
    VecDouble rhs(16);
    rhs.setZero();

    int64_t nConnects = cmeshG->NConnects();
    int64_t nEnrichedNodes = enrichedNodes.size();
    cmeshG->ConnectVec().resize(nConnects + nEnrichedNodes);
    int order = cmeshG->GetDefaultOrder();
    int count = 0;
    for(auto node:enrichedNodes){
        int64_t prevSeqNum = cmeshG->ConnectVec()[nConnects+count-1]->GetSequenceNumber();
        int64_t seqnum = prevSeqNum + (cmeshG->ConnectVec()[nConnects+count]->GetNShapeFunctions()*dimension);
        Connect* c = new Connect(dimension, 1, order, nConnects+count, seqnum);

    }

    GlobalLocalEnrichment *globalLocal = new GlobalLocalEnrichment(1,dimension,1.0,0.0);
    ElementEnriched *teste = new ElementEnriched(1, localEl, globalEl, cmeshG, globalLocal);
    teste->setCorrespondence(globalElementCorrespondence, globalNodeCorrespondence);
    teste->ComputeElContribution(stiffness,rhs);

}   

void LocalToGlobalCorrespondence(CompMesh *cmeshG, CompMesh *cmeshL){

    Element *localElement = nullptr;
    Element *globalElement = nullptr;

    VecDouble x(dimension);
    VecDouble xsiCorr(2);

    MatrixDouble globalNode;
    int64_t elCorr;

    for (int i=0; i<cmeshL->NElements(); i++){
        if(cmeshL->ElementVec()[i]->Dimension() == cmeshL->Dimension() || cmeshL->ElementVec()[i]->Reference()->Material() == overlappingBoundary){
            localElement = cmeshL->ElementVec()[i];

            int numberOfIntegrationPoints = localElement -> getNumberOfIntegrationPoints();
            localElement -> ComputeIntegPointCoordinates();

            globalNode.resize(numberOfIntegrationPoints, 2);
            globalNode.setZero();

            elCorr = 0;

            for(int index = 0; index < numberOfIntegrationPoints; index++){
                x.setZero();
                xsiCorr.setZero();

                for(int k = 0; k < dimension; k++) x[k] = localElement->getIntegPointCoordinatesValue(index)[k];

                for (elCorr; elCorr<cmeshG->NElements(); elCorr++){
                    if(cmeshG->ElementVec()[elCorr]->Reference()->Material() != overlappingRegion) continue;
                    if(cmeshG->ElementVec()[elCorr]->Dimension() != cmeshG->Dimension()) continue;
                    globalElement = cmeshG->ElementVec()[elCorr];

                    if(CompMeshTools ::searchNodeCorrespondence(x, cmeshG, elCorr, xsiCorr, globalElement)){
                        globalNode(index, 0) = xsiCorr[0]; globalNode(index, 1) = xsiCorr[1];
                        break;
                    }else{
                        elCorr = globalElement->Index();
                    }                 
                }        
            }

            globalNodeCorrespondence[localElement->Index()] = globalNode;
            globalElementCorrespondence[localElement->Index()] = elCorr;
        }
    }
}

void CreateGlobalModel(CompMesh *cmeshG){

    Elasticity2D *matelasticityG1 = new Elasticity2D(1, 1.0, 0.0); //região de sobreposição do domínio local no domínio global
    Elasticity2D *matelasticityG2 = new Elasticity2D(2, 1.0, 0.0); //domínio global complementar
    cmeshG->InsertMaterial(matelasticityG1);
    cmeshG->InsertMaterial(matelasticityG2);

    overlappingRegion = 1;
    enrichedNodes = {0,5};
    
    //BC 
    MatrixDouble val1(2,2);
    val1.setZero();
    VecDouble val2(2);
    val2.setZero();
    L2Projection * matbcG1 = new L2Projection(3,dimension-1,BoundaryConditionType::kDirichlet,val1,val2);
    val2[1] = 1.0;
    L2Projection * matbcG2 = new L2Projection(4,dimension-1,BoundaryConditionType::kDirectionalHomogeneousDirichlet,val1,val2);
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
    //anG.PrintSolution();

    VTUGenerator::PrintResults(cmeshG,"globalResult",ScalarNames,VectorNames);
}

void CreateLocalModel(CompMesh *cmeshL){

    Elasticity2D* matelasticityL = new Elasticity2D(1, 1.0, 0.0);
    cmeshL->InsertMaterial(matelasticityL);

    overlappingBoundary = 2;

    //BC;
    MatrixDouble val1(2,2);
    val1.setZero();
    VecDouble val2(2);
    val2[1] = 1.0;
    L2Projection * matbcL1 = new L2Projection(3,dimension-1,BoundaryConditionType::kDirectionalHomogeneousDirichlet,val1,val2);
    val2.setZero();
    val2[0] = -2.0;
    L2Projection * matbcL2 = new L2Projection(2,dimension-1,BoundaryConditionType::kDirectionalNonHomogeneousDirichlet,val1,val2);
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
    //anL.PrintSolution();

    std::vector<std::string> ScalarNames, VectorNames;
    ScalarNames = {"SigmaX","SigmaY","TauXY"};
    VectorNames = {"Displacement"};

    VTUGenerator::PrintResults(cmeshL,"localResult",ScalarNames,VectorNames); 
};