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

//Comentar erro ao debugar após incluir a weakform GlobalLocalEnrichment e o element ElementEnriched. 

// Defines the problem dimension
const int dimension = 2;
auto forcingFunction = [](const VecDouble &coord, VecDouble &force){
    const auto &x=coord[0];
    const auto &y=coord[1];
    // force[0] = -1*y;
    force[0] = -1;
};

MixedCompMesh* CreateGlobalModel(GeoMesh *gmesh, std::vector<CompMesh *> &meshvector);
MixedCompMesh* CreateLocalModel(GeoMesh *gmesh, std::vector<CompMesh *> &meshvector, CompMesh *cmeshG);
void SolveGlobalProblem(CompMesh *cmeshG);
void SolveLocalProblem(CompMesh *cmeshL);
void LocalToGlobalCorrespondence(CompMesh *cmeshG, CompMesh *cmeshL);
void LocalToGlobalCorrespondenceBoundary(CompMesh *cmeshG, CompMesh *cmeshL);
void CreateEnrichedModel(CompMesh *cmeshG, CompMesh *cmeshL);
void SolveEnrichedProblem(CompMesh *cmeshG);

int overlappingRegion;
int overlappingNHDirichletBoundary;
int overlappingNHNeumannBoundary;
int globalLocalIterations;
//Local index to global index correspondence for elements
std::map<int64_t,int64_t> globalElementCorrespondence;
//Local node/integration point to global node/integration point correspondence 
std::map<int64_t, MatrixDouble> globalNodeCorrespondence;
//Global node to be enriched to the index of the new connect associated with the enriched degree of freedom in the global mesh
std::map<int64_t,int64_t> enrichedNodes;

int main(int argc, char **args) { 

    overlappingRegion = 1;
    globalLocalIterations = 1;
    overlappingNHDirichletBoundary = 2;
    overlappingNHNeumannBoundary = 4;
    enrichedNodes[0]=-1;
    enrichedNodes[5]=-1;

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

    
    GeoMesh * gmeshL = new GeoMesh();
    GmshTools::Read(*gmeshL,"../local.msh");
    CompMesh *cmeshDispL, *cmeshPressL;
    std::vector<CompMesh *> meshvectorL = {cmeshDispL, cmeshPressL};
    MixedCompMesh *cmeshL = CreateLocalModel(gmeshL, meshvectorL,cmeshG);
    gmeshL->Print("gmeshLocal.txt");
    cmeshL->Print("cmeshLocal.txt");

    LocalToGlobalCorrespondence(cmeshG,cmeshL);
    LocalToGlobalCorrespondenceBoundary(cmeshG,cmeshL);

    for(int n = 0; n < globalLocalIterations; n++){

        std::cout << "\nGlobal-Local Iteration " << n + 1 << ": ";

        //LOCAL MODEL
        std::cout << "\nSolve Local Problem \n";
        //Solve Local Problem
        SolveLocalProblem(cmeshL);

        //ENRICHED GLOBAL MODEL
        std::cout << "\nSolve Enriched Global Problem \n";
        if(n == 0) CreateEnrichedModel(cmeshG,cmeshL);

        cmeshG->Print("cmeshEnrichedGlobal.txt");
        SolveEnrichedProblem(cmeshG);
   }
}   

MixedCompMesh* CreateGlobalModel(GeoMesh *gmesh, std::vector<CompMesh *> &meshvector){

    //Displacement Cmesh
    meshvector[0] = new CompMesh(gmesh, ApproxType::EHierarquic);
    meshvector[0]->SetDefaultOrder(2);
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
    // val3[1] = 1.0;
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
    anG.PrintSolution();
    //anG.PrintGlobalRhs();

    // VTUGenerator::PrintResults(cmeshG,"globalResult",ScalarNames,VectorNames);
}

MixedCompMesh* CreateLocalModel(GeoMesh *gmesh, std::vector<CompMesh *> &meshvector, CompMesh *cmeshG){
    //Displacement Cmesh
    meshvector[0] = new CompMesh(gmesh, ApproxType::EHierarquic);
    meshvector[0]->SetDefaultOrder(2);

    int nstate = 2;
    WeakForm* matelasticityL = new WeakForm(1, nstate);
    meshvector[0]->InsertMaterial(matelasticityL);

    //BC;
    MatrixDouble val1(2,2);
    val1.setZero();
    VecDouble val2(2);
    val2[1] = 1.0;
    L2Projection * matbcL1 = new L2Projection(3,dimension-1,BoundaryConditionType::kDirectionalHomogeneousDirichlet,val1,val2);
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
    anL.PrintGlobalMatrix();
    anL.PrintGlobalRhs();
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

    // VTUGenerator::PrintResults(cmeshL,"localResult",ScalarNames,VectorNames); 
};


void LocalToGlobalCorrespondenceBoundary(CompMesh *cmeshG, CompMesh *cmeshL){

    Element *localElement = nullptr;
    Element *globalElement = nullptr;

    VecDouble x(dimension);
    VecDouble xsiCorr(1);

    MatrixDouble globalNode;
    int64_t elCorr;

    for (int i=0; i<cmeshL->NElements(); i++){
        localElement = cmeshL->ElementVec()[i];
        if (!localElement) continue;
        if(cmeshL->ElementVec()[i]->Dimension() == cmeshL->Dimension()) continue;
        if (localElement->Reference()->Material() != overlappingNHNeumannBoundary) continue;

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
                if(cmeshG->ElementVec()[elCorr]->Dimension() == cmeshG->Dimension()) continue;
                globalElement = cmeshG->ElementVec()[elCorr];

                if(CompMeshTools ::searchNodeCorrespondence(x, cmeshG, elCorr, xsiCorr, globalElement)){
                    globalElement->IntegrationData().fA0.resize(1,1);
                    globalElement->IntegrationData().fA0Inv.resize(1,1);
                    globalNode(index, 0) = xsiCorr[0];
                    break;
                }else{
                    globalElement->IntegrationData().fA0.resize(1,1);
                    globalElement->IntegrationData().fA0Inv.resize(1,1);
                    elCorr = globalElement->Index();
                }        
            }        
        }

        globalNodeCorrespondence[localElement->Index()] = globalNode;
        globalElementCorrespondence[localElement->Index()] = elCorr;
    }
};


void LocalToGlobalCorrespondence(CompMesh *cmeshG, CompMesh *cmeshL){

    Element *localElement = nullptr;
    Element *globalElement = nullptr;

    VecDouble x(dimension);
    VecDouble xsiCorr(2);

    MatrixDouble globalNode;
    int64_t elCorr;

    for (int i=0; i<cmeshL->NElements(); i++){
        if(cmeshL->ElementVec()[i]->Dimension() == cmeshL->Dimension() || cmeshL->ElementVec()[i]->Reference()->Material() == overlappingNHDirichletBoundary){
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
};

void CreateEnrichedModel(CompMesh *cmeshG, CompMesh *cmeshL){

    // Create the new connects in the global mesh for the enriched nodes and resize the connect vector of the global mesh accordingly. 
    // The number of new connects is equal to the number of enriched nodes, since we are considering only one degree of freedom per node, 
    // but it can be easily generalized for more degrees of freedom per node.
    int64_t nConnects = cmeshG->NConnects();
    int64_t nEnrichedNodes = enrichedNodes.size();
    cmeshG->ConnectVec().resize(nConnects + nEnrichedNodes);
    int order = cmeshG->GetDefaultOrder();
    int count = 0;
    int64_t SeqNum = cmeshG->NGlobalDOF();
    int nstate = cmeshG->NState();
    for(auto &node:enrichedNodes){;
        Connect* c = new Connect(dimension, 1, order, nConnects+count, SeqNum);
        SeqNum += nstate;
        node.second = nConnects+count;
        cmeshG->ConnectVec()[nConnects+count] = c;
        count++;
    }

    //cmeshG->Print("cmeshGEnriched.txt");
    
    GlobalLocalEnrichment *globalLocal = new GlobalLocalEnrichment(1,dimension,1.0,0.0);
    globalLocal->SetForcingFunction(forcingFunction);

    //Count the number local elements to enrich in the global mesh
    int nelsToEnrich = 0;
    for (auto el:cmeshL->ElementVec()){
        if (el->Dimension() != cmeshL->Dimension() && el->Reference()->Material() != overlappingNHNeumannBoundary) continue;
        nelsToEnrich++;
    }
    int64_t nElementsG = cmeshG->NElements();
    cmeshG->ElementVec().resize(cmeshG->NElements() + nelsToEnrich);
    //Create the enriched elements in the global mesh
    count = 0;
    for (auto localEl:cmeshL->ElementVec()){
        if (localEl->Dimension() != cmeshL->Dimension() && localEl->Reference()->Material() != overlappingNHNeumannBoundary) continue;
        Element* globalEl = cmeshG->ElementVec()[globalElementCorrespondence[localEl->Index()]];
        ElementEnriched *enrichedEl = new ElementEnriched(nElementsG+count, localEl, globalEl, cmeshG, globalLocal);
        enrichedEl->setCorrespondence(&globalElementCorrespondence, &globalNodeCorrespondence);
        
        //Sets which node will have enriched solution
        enrichedEl->SetEnrichmentData(&enrichedNodes);

        //Remove global element weak form, for skipping it when contributing in the global stiffness matrix and rhs.
        if (globalEl->Dimension() == cmeshG->Dimension()){
            globalEl->SetWeakForm(nullptr);
        }

        //Seek how many nodes will be enriched in the global element and construct the proper connectivity for the enriched element.
        VecInt geoNodes = globalEl->Reference()->getGeometricNodes();
        VecInt enrichedConnects = globalEl->getConnectivityIndices();
        for (int i = 0; i < geoNodes.size(); i++){
            if (enrichedNodes.find(geoNodes[i]) != enrichedNodes.end()){
                enrichedConnects.conservativeResize(enrichedConnects.size() + 1); // Increase size by 1
                enrichedConnects(enrichedConnects.size() - 1) = enrichedNodes[geoNodes[i]];          
            }
        }
        enrichedEl->getConnectivity().resize(enrichedConnects.size());
        enrichedEl->setConnectivity(enrichedConnects);
        cmeshG->ElementVec()[nElementsG+count] = enrichedEl;
        count++;
    }
    // Update the problem size
    int64_t fNGlobalDOF  = 0;
    for (int64_t i = 0; i < cmeshG->NConnects(); i++){
        int nstate = cmeshG->NState();
        fNGlobalDOF += cmeshG->ConnectVec()[i]->GetNShapeFunctions() * nstate;
    }
    cmeshG->NGlobalDOF() = fNGlobalDOF;

    //cmeshG->Print("cmeshGEnriched2.txt");
};

void SolveEnrichedProblem(CompMesh *cmeshG){

    LinearAnalysis anE(cmeshG,SolverType::ELU);

    anE.Run();
    //anE.PrintGlobalMatrix();
    //anE.PrintGlobalRhs();
    anE.PrintSolution();


    std::vector<std::string> ScalarNames, VectorNames;
    // ScalarNames = {"SigmaX","SigmaY","TauXY"};
    VectorNames = {"Displacement"};

    // VTUGenerator::PrintResults(cmeshG,"Enriched",ScalarNames,VectorNames); 
};