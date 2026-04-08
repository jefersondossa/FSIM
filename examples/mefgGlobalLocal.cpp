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

//Comentar erro ao debugar após incluir a weakform GlobalLocalEnrichment e o element ElementEnriched. 

// Defines the problem dimension
const int dimension = 2;
auto forcingFunction = [](const VecDouble &coord, VecDouble &force){
    const auto &x=coord[0];
    const auto &y=coord[1];

    force[0] = -1*y;
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
VecDouble previousSol;
VecDouble currentSol;
//Local index to global index correspondence for elements
std::map<int64_t,int64_t> globalElementCorrespondence;
//Local node/integration point to global node/integration point correspondence 
std::map<int64_t, MatrixDouble> globalNodeCorrespondence;
//Global node to be enriched to the index of the new connect associated with the enriched degree of freedom in the global mesh
std::map<int64_t,int64_t> enrichedNodes;

int main(int argc, char **args) { 

    //GLOBAL MODEL
    std::cout << "Solve Global Problem \n";

    //Create Global Model
    GeoMesh *gmeshG = new GeoMesh();
    GmshTools::Read(*gmeshG,"../global.msh");
    CompMesh *cmeshG = new CompMesh(gmeshG,ApproxType::EIsoparametric);
    CreateGlobalModel(cmeshG);
    gmeshG->Print("gmeshGlobal.txt");
    cmeshG->Print("cmeshGlobal.txt");
    
    //Solve Global Problem
    SolveGlobalProblem(cmeshG);

    CompMesh *cmeshL = nullptr;
    GeoMesh * gmeshL = new GeoMesh();

    int it = 0;
    //for(int it = 0; it < globalLocalIterations; it++){
    while(it < globalLocalIterations){
        std::cout << "\nGlobal-Local Iteration " << it + 1 << ": ";

        //LOCAL MODEL
        std::cout << "\nSolve Local Problem \n";
        if(it == 0){
            //Create Local Model
            GmshTools::Read(*gmeshL,"../local.msh");
            cmeshL = new CompMesh(gmeshL,ApproxType::EIsoparametric);
            CreateLocalModel(cmeshG, cmeshL);
            //gmeshL->Print("gmeshLocal.txt");
            //cmeshL->Print("cmeshLocal.txt");

            LocalToGlobalCorrespondence(cmeshG,cmeshL);
            LocalToGlobalCorrespondenceBoundary(cmeshG,cmeshL);
        }
        //Solve Local Problem
        SolveLocalProblem(cmeshL);

        //ENRICHED GLOBAL MODEL
        std::cout << "\nSolve Enriched Global Problem \n";
        if(it == 0) CreateEnrichedModel(cmeshG,cmeshL);

        //cmeshG->Print("cmeshEnrichedGlobal.txt");
        SolveEnrichedProblem(cmeshG);

        if(CheckConvergence(it)) break;
        it++;
   };
}   

void CreateGlobalModel(CompMesh *cmeshG){

    Elasticity2D *matelasticityG1 = new Elasticity2D(1, 1.0, 0.0); //domínio global
    Elasticity2D *matelasticityG2 = new Elasticity2D(2, 1.0, 0.0); //domínio local

    //matelasticityG1->SetExactSolution(exactSol);

    cmeshG->InsertMaterial(matelasticityG1);
    cmeshG->InsertMaterial(matelasticityG2);

    overlappingRegion = 2;
    globalLocalIterations = 1;
    globalLocalTolerance = 1e-4;

    //Chapa retangular tracionada
    enrichedNodes[0]=-1;
    enrichedNodes[6]=-1;
    enrichedNodes[7]=-1;
    enrichedNodes[5]=-1;

    //Chapa L 
    // enrichedNodes[8]=-1;
    // enrichedNodes[28]=-1;
    // enrichedNodes[29]=-1;
    // enrichedNodes[51]=-1;
    // enrichedNodes[52]=-1;
    // enrichedNodes[62]=-1;
    // enrichedNodes[63]=-1;
    // enrichedNodes[64]=-1;
    
    //BC 
    MatrixDouble val1(2,2);
    val1.setZero();
    VecDouble val2(2);
    val2.setZero();

    //Chapa L teste: apoio fixo no lado BC e carregamento unitário uniforme ao longe de FA
    //L2Projection * matbcG1 = new L2Projection(3,dimension-1,BoundaryConditionType::kDirichlet,val1,val2);
    //L2Projection * matbcG2 = new L2Projection(4,dimension-1,BoundaryConditionType::kNeumann,val1,val2);
    //matbcG2->SetForcingFunction(forcingFunction);

    //Chapa L: apoio horizontal no pontos B e F, apoio vertical no ponto A
    // val2[0] = 1.0;
    // L2Projection * matbcG1 = new L2Projection(3,dimension-2,BoundaryConditionType::kDirectionalHomogeneousDirichlet,val1,val2);
    // val2.setZero();
    // val2[1] = 1.0;
    // L2Projection * matbcG2 = new L2Projection(4,dimension-2,BoundaryConditionType::kDirectionalHomogeneousDirichlet,val1,val2);
    
    //Chapa L: apoio fixo no ponto D, apoio vertical no ponto A
    //L2Projection * matbcG1 = new L2Projection(3,dimension-2,BoundaryConditionType::kDirichlet,val1,val2);
    //val2[0] = 1.0;
    //L2Projection * matbcG2 = new L2Projection(4,dimension-2,BoundaryConditionType::kDirectionalHomogeneousDirichlet,val1,val2);
    
    // L2Projection * matbcG3 = new L2Projection(5,dimension-1,BoundaryConditionType::kNeumann,val1,val2);
    // matbcG3->SetForcingFunction(forcingFunctionAB);
    // L2Projection * matbcG4 = new L2Projection(6,dimension-1,BoundaryConditionType::kNeumann,val1,val2);
    // matbcG4->SetForcingFunction(forcingFunctionBC);
    // L2Projection * matbcG5 = new L2Projection(7,dimension-1,BoundaryConditionType::kNeumann,val1,val2);
    // matbcG5->SetForcingFunction(forcingFunctionEF);
    // L2Projection * matbcG6 = new L2Projection(8,dimension-1,BoundaryConditionType::kNeumann,val1,val2);
    // matbcG6->SetForcingFunction(forcingFunctionFA);

    //Chapa retangular tracionada
    L2Projection * matbcG1 = new L2Projection(3,dimension-1,BoundaryConditionType::kDirichlet,val1,val2);
    val2[1] = 1.0;
    L2Projection * matbcG2 = new L2Projection(4,dimension-1,BoundaryConditionType::kDirectionalHomogeneousDirichlet,val1,val2);
    L2Projection * matbcG3 = new L2Projection(5,dimension-1,BoundaryConditionType::kNeumann,val1,val2);
    matbcG3->SetForcingFunction(forcingFunction);

    cmeshG->InsertMaterial(matbcG1);
    cmeshG->InsertMaterial(matbcG2);
    cmeshG->InsertMaterial(matbcG3);
    // cmeshG->InsertMaterial(matbcG4);
    // cmeshG->InsertMaterial(matbcG5);
    // cmeshG->InsertMaterial(matbcG6);
    
    cmeshG->AutoBuild();
}

void SolveGlobalProblem(CompMesh *cmeshG){

    LinearAnalysis anG(cmeshG,SolverType::ELDLt);
       
    std::vector<std::string> ScalarNames, VectorNames;
    ScalarNames = {"SigmaX","SigmaY","TauXY"};
    VectorNames = {"Displacement"}; //, "ExactDisplacement"}

    anG.Run();
    //anG.PrintGlobalMatrix();
    //anG.PrintGlobalRhs();
    anG.PrintSolution();

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
    std::cout << "Strain Energy: "<< strainEnergy << std::endl;

    VTUGenerator::PrintResults(cmeshG,"globalResult",ScalarNames,VectorNames);
}

void CreateLocalModel(CompMesh *cmeshG, CompMesh * cmeshL){

    Elasticity2D* matelasticityL = new Elasticity2D(1, 1.0, 0.0);
    cmeshL->InsertMaterial(matelasticityL);

    overlappingNHDirichletBoundary = 2;
    overlappingNHNeumannBoundary = 4;

    //BC;
    MatrixDouble val1(2,2);
    val1.setZero();
    VecDouble val2(2);
    val2.setZero();
    
    InterpolatedBC * matbcL1 = new InterpolatedBC(2,dimension-1,2,BoundaryConditionType::kDirichlet,&globalElementCorrespondence,&globalNodeCorrespondence,cmeshG);
    //Chapa L: apoio fixo no ponto D, apoio vertical no ponto A
    //L2Projection * matbcL2 = new L2Projection(3,dimension-2,BoundaryConditionType::kDirichlet,val1,val2);
    
    //Chapa reatangular tracionada
    val2[1] = 1.0;
    L2Projection * matbcL2 = new L2Projection(3,dimension-1,BoundaryConditionType::kDirectionalHomogeneousDirichlet,val1,val2);
    val2.setZero();
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
    anL.PrintSolution();

    std::vector<std::string> ScalarNames, VectorNames;
    ScalarNames = {"SigmaX","SigmaY","TauXY"};
    VectorNames = {"Displacement"};

    VTUGenerator::PrintResults(cmeshL,"localResult",ScalarNames,VectorNames); 
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