#include "MEFGGlobalLocalTools.h"
#include "GeoMesh.h"
#include "CompMesh.h"
#include "Elasticity2D.h"
#include "L2Projection.h"
#include "CompMeshTools.h"
#include "ElementEnriched.h"
#include "ElementMixed.h"

void MEFGGlobalLocalTools::LocalToGlobalCorrespondence(CompMesh *cmeshG, CompMesh *cmeshL, 
                                                       std::map<int64_t,int64_t> &globalElementCorrespondence, 
                                                       std::map<int64_t, MatrixDouble> &globalNodeCorrespondence,
                                                       int overlappingRegion, int overlappingNHDirichletBoundary){
    
    Element *localElement = nullptr;
    Element *globalElement = nullptr;
    int dimension = cmeshG->Dimension();

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
                xsiCorr.fill(1e50);

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

void MEFGGlobalLocalTools::LocalToGlobalCorrespondenceBoundary(CompMesh *cmeshG, CompMesh *cmeshL,
                                                           std::map<int64_t,int64_t> &globalElementCorrespondence,
                                                           std::map<int64_t, MatrixDouble> &globalNodeCorrespondence,
                                                           int overlappingNHNeumannBoundary){

    Element *localElement = nullptr;
    Element *globalElement = nullptr;
    int dimension = cmeshG->Dimension();

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


void MEFGGlobalLocalTools::CreateEnrichedModel(CompMesh *cmeshG, CompMesh *cmeshL, 
                                               std::map<int64_t,int64_t> &globalElementCorrespondence, 
                                               std::map<int64_t, MatrixDouble> &globalNodeCorrespondence,
                                               std::map<int64_t,int64_t> &enrichedConnects,
                                               int overlappingNHNeumannBoundary,
                                               GlobalLocalEnrichment *globalLocal){

    int dimension = cmeshG->Dimension();

    // Create the new connects in the global mesh for the enriched nodes and resize the connect vector of the global mesh accordingly. 
    // The number of new connects is equal to the number of enriched nodes, since we are considering only one degree of freedom per node, 
    // but it can be easily generalized for more degrees of freedom per node.
    int64_t nConnects = cmeshG->NConnects();
    int64_t nEnrichedConnects = enrichedConnects.size();
    cmeshG->ConnectVec().resize(nConnects + nEnrichedConnects);
    int count = 0;
    int64_t SeqNum = cmeshG->NGlobalDOF();
    int nstate = cmeshG->NState();
    for(auto &con:enrichedConnects){;
        Connect* originalConnect = cmeshG->ConnectVec()[con.first];
        int nshape = originalConnect->GetNShapeFunctions();
        int order = originalConnect->GetOrder();
        Connect* c = new Connect(dimension, nshape, order, nConnects+count, SeqNum);
        SeqNum += nstate;
        con.second = nConnects+count;
        cmeshG->ConnectVec()[nConnects+count] = c;
        count++;
    }

    //cmeshG->Print("cmeshGEnriched.txt");

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
        enrichedEl->SetEnrichmentData(&enrichedConnects);

        //Remove global element weak form, for skipping it when contributing in the global stiffness matrix and rhs.
        if (globalEl->Dimension() == cmeshG->Dimension()){
            globalEl->SetWeakForm(nullptr);
        }

        //Seek how many connects will be enriched in the global element and construct the proper connectivity for the enriched element.
        auto elConnects = globalEl->getConnectivity();
        VecInt enrichedCon = globalEl->getConnectivityIndices();
        for (int i = 0; i < elConnects.size(); i++){
            if (enrichedConnects.find(elConnects[i]->Index()) != enrichedConnects.end()){
                enrichedCon.conservativeResize(enrichedCon.size() + 1); // Increase size by 1
                enrichedCon(enrichedCon.size() - 1) = enrichedConnects[elConnects[i]->Index()];          
            }
        }
        enrichedEl->getConnectivity().resize(enrichedCon.size());
        enrichedEl->setConnectivity(enrichedCon);
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
    

};


void MEFGGlobalLocalTools::CreateMixedEnrichedModel(CompMesh *cmeshG, CompMesh *cmeshL, 
                                               std::map<int64_t,int64_t> &globalElementCorrespondence, 
                                               std::map<int64_t, MatrixDouble> &globalNodeCorrespondence,
                                               std::map<int64_t,int64_t> &enrichedConnects,
                                               int overlappingNHNeumannBoundary,
                                               MixedGlobalLocalEnrichment *globalLocal){

    MixedCompMesh *mixedCmeshG = dynamic_cast<MixedCompMesh *>(cmeshG);
    MixedCompMesh *mixedCmeshL = dynamic_cast<MixedCompMesh *>(cmeshL);

    // Create the new connects in the global mesh for the enriched nodes and resize the connect vector of the global mesh accordingly. 
    // The number of new connects is equal to the number of enriched nodes, since we are considering only one degree of freedom per node, 
    // but it can be easily generalized for more degrees of freedom per node.
    int64_t nConnects = mixedCmeshG->MeshVector()[0]->NConnects();
    int64_t nEnrichedConnects = enrichedConnects.size();
    mixedCmeshG->MeshVector()[0]->ConnectVec().resize(nConnects + nEnrichedConnects);
    int count = 0;
    int64_t SeqNum = cmeshG->NGlobalDOF();
    for(auto &con:enrichedConnects){;
        Connect* originalConnect = mixedCmeshG->MeshVector()[0]->ConnectVec()[con.first];
        int nshape = originalConnect->GetNShapeFunctions();
        int order = originalConnect->GetOrder();
        int nstate = originalConnect->GetNStateVariables();
        Connect* c = new Connect(nstate, nshape, order, nConnects+count, SeqNum);
        SeqNum += nstate;
        con.second = nConnects+count;
        mixedCmeshG->MeshVector()[0]->ConnectVec()[nConnects+count] = c;
        count++;
    }

    // mixedCmeshG->MeshVector()[0]->Print("cmeshGEnriched.txt");

    //Count the number local elements to enrich in the global mesh
    int nelsToEnrich = 0;
    for (auto el:cmeshL->ElementVec()){
        if (el->Dimension() != cmeshL->Dimension() && el->Reference()->Material() != overlappingNHNeumannBoundary) continue;
        nelsToEnrich++;
    }
    int64_t nElementsG = mixedCmeshG->MeshVector()[0]->NElements();
    mixedCmeshG->MeshVector()[0]->ElementVec().resize(mixedCmeshG->MeshVector()[0]->NElements() + nelsToEnrich);
    cmeshG->ElementVec().resize(cmeshG->NElements() + nelsToEnrich);
    //Create the enriched elements in the global mesh
    count = 0;
    for (auto localEl:cmeshL->ElementVec()){
        if (localEl->Dimension() != cmeshL->Dimension() && localEl->Reference()->Material() != overlappingNHNeumannBoundary) continue;
        Element* globalEl = cmeshG->ElementVec()[globalElementCorrespondence[localEl->Index()]];
        ElementEnriched *enrichedEl = new ElementEnriched(nElementsG+count, localEl, globalEl, mixedCmeshG, globalLocal);
        enrichedEl->setCorrespondence(&globalElementCorrespondence, &globalNodeCorrespondence);
        ElementMixed* globalElMixed = dynamic_cast<ElementMixed*>(globalEl);
        // std::vector<Element *> elMixedEnriched = {enrichedEl, localEl};
        // ElementMixed *enrichedElMixed = new ElementMixed(nElementsG+count, elMixedEnriched, mixedCmeshG, globalLocal);

        //Sets which node will have enriched solution
        enrichedEl->SetEnrichmentData(&enrichedConnects);

        //Remove global element weak form, for skipping it when contributing in the global stiffness matrix and rhs.
        if (globalEl->Dimension() == mixedCmeshG->MeshVector()[0]->Dimension()){
            mixedCmeshG->MeshVector()[0]->ElementVec()[globalElementCorrespondence[localEl->Index()]]->SetWeakForm(nullptr);
            globalEl->SetWeakForm(nullptr);
        }

        //Seek how many connects will be enriched in the global element and construct the proper connectivity for the enriched element.
        auto elConnects = globalEl->getConnectivity();
        std::vector<Connect *> enrichedCon = globalElMixed->SubElements()[0]->getConnectivity();
        for (int i = 0; i < elConnects.size(); i++){
            if (enrichedConnects.find(elConnects[i]->Index()) != enrichedConnects.end()){
                enrichedCon.resize(enrichedCon.size() + 1); // Increase size by 1
                enrichedCon[enrichedCon.size() - 1] = mixedCmeshG->MeshVector()[0]->ConnectVec()[enrichedConnects[elConnects[i]->Index()]];          
            }
        }
        std::vector<Connect *> presscon = globalElMixed->SubElements()[1]->getConnectivity();
        enrichedCon.insert(enrichedCon.end(), presscon.begin(), presscon.end());
        
        enrichedEl->setConnectivity(enrichedCon);
        mixedCmeshG->MeshVector()[0]->ElementVec()[nElementsG+count] = enrichedEl;
        cmeshG->ElementVec()[nElementsG+count] = enrichedEl;
        count++;
    }
    // Update the problem size
    int64_t fNGlobalDOF  = 0;
    for (int imesh = 0; imesh < mixedCmeshG->MeshVector().size(); imesh++){
        for (int64_t i = 0; i < mixedCmeshG->MeshVector()[imesh]->NConnects(); i++){
            Connect* con = mixedCmeshG->MeshVector()[imesh]->ConnectVec()[i];
            int nstate = con->GetNStateVariables();
            fNGlobalDOF += mixedCmeshG->MeshVector()[imesh]->ConnectVec()[i]->GetNShapeFunctions() * nstate;
        }
    }
    cmeshG->NGlobalDOF() = fNGlobalDOF;

}