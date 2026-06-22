#include "MEFGGlobalLocalTools.h"
#include "GeoMesh.h"
#include "CompMesh.h"
#include "Elasticity2D.h"
#include "L2Projection.h"
#include "CompMeshTools.h"
#include "ElementEnriched.h"
#include "ElementMixed.h"

void MEFGGlobalLocalTools::LocalToGlobalCorrespondence(CompMesh *cmeshG, CompMesh *cmeshL, 
                                                       std::map<int,int> &globalElementCorrespondence, 
                                                       std::map<int, MatrixDouble> &globalNodeCorrespondence,
                                                       int overlappingRegion, int overlappingNHDirichletBoundary){
    
    Element *localElement = nullptr;
    Element *globalElement = nullptr;
    int dimension = cmeshG->Dimension();

    VecDouble x(dimension);
    VecDouble xsiCorr(2);

    MatrixDouble globalNode;
    int elCorr;

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
                                                           std::map<int,int> &globalElementCorrespondence,
                                                           std::map<int, MatrixDouble> &globalNodeCorrespondence,
                                                           int overlappingNHNeumannBoundary){

    Element *localElement = nullptr;
    Element *globalElement = nullptr;
    int dimension = cmeshG->Dimension();

    VecDouble x(dimension);
    VecDouble xsiCorr(1);

    MatrixDouble globalNode;
    int elCorr;

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
                                               std::map<int,int> &globalElementCorrespondence, 
                                               std::map<int, MatrixDouble> &globalNodeCorrespondence,
                                               std::map<int,int> &enrichedConnects,
                                               int overlappingNHNeumannBoundary,
                                               GlobalLocalEnrichment *globalLocal){

    int dimension = cmeshG->Dimension();

    // Create the new connects in the global mesh for the enriched nodes and resize the connect vector of the global mesh accordingly. 
    // The number of new connects is equal to the number of enriched nodes, since we are considering only one degree of freedom per node, 
    // but it can be easily generalized for more degrees of freedom per node.
    int nConnects = cmeshG->NConnects();
    int nEnrichedConnects = enrichedConnects.size();
    cmeshG->ConnectVec().resize(nConnects + nEnrichedConnects);
    int count = 0;
    int SeqNum = cmeshG->NGlobalDOF();
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
    int nElementsG = cmeshG->NElements();
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
    int fNGlobalDOF  = 0;
    for (int i = 0; i < cmeshG->NConnects(); i++){
        int nstate = cmeshG->NState();
        fNGlobalDOF += cmeshG->ConnectVec()[i]->GetNShapeFunctions() * nstate;
    }
    cmeshG->NGlobalDOF() = fNGlobalDOF;
    

};

void MEFGGlobalLocalTools::CreateMixedEnrichedModel(
    CompMesh *cmeshG,
    CompMesh *cmeshL,
    std::map<int,int> &globalElementCorrespondence,
    std::map<int, MatrixDouble> &globalNodeCorrespondence,
    std::map<int,int> &enrichedConnects,
    int overlappingNHNeumannBoundary,
    MixedGlobalLocalEnrichment *globalLocal)
{
    if (!cmeshG || !cmeshL || !globalLocal) {
        PanicButton();
        return;
    }

    MixedCompMesh *mixedCmeshG = dynamic_cast<MixedCompMesh *>(cmeshG);
    MixedCompMesh *mixedCmeshL = dynamic_cast<MixedCompMesh *>(cmeshL);

    if (!mixedCmeshG || !mixedCmeshL) {
        PanicButton();
        return;
    }

    if (mixedCmeshG->MeshVector().empty()) {
        PanicButton();
        return;
    }

    CompMesh *dispMeshG = mixedCmeshG->MeshVector()[0];

    if (!dispMeshG) {
        PanicButton();
        return;
    }

    // ------------------------------------------------------------
    // 1. Create enriched connects in the displacement mesh
    // ------------------------------------------------------------

    const int nConnects = dispMeshG->NConnects();
    const int nEnrichedConnects =
        static_cast<int>(enrichedConnects.size());

    dispMeshG->ConnectVec().resize(nConnects + nEnrichedConnects);

    int seqNum = cmeshG->NGlobalDOF();
    int count = 0;

    for (auto &con : enrichedConnects) {

        const int originalConnectIndex = con.first;

        if (originalConnectIndex < 0 || originalConnectIndex >= nConnects) {
            PanicButton();
            return;
        }

        Connect *originalConnect =
            dispMeshG->ConnectVec()[originalConnectIndex];

        if (!originalConnect) {
            PanicButton();
            return;
        }

        const int nshape = originalConnect->GetNShapeFunctions();
        const int order  = originalConnect->GetOrder();
        const int nstate = originalConnect->GetNStateVariables();

        const int newConnectIndex = nConnects + count;

        Connect *newConnect =
            new Connect(nstate, nshape, order, newConnectIndex, seqNum);

        seqNum += nshape * nstate;

        con.second = newConnectIndex;
        dispMeshG->ConnectVec()[newConnectIndex] = newConnect;

        count++;
    }

    // ------------------------------------------------------------
    // 2. Count local elements to enrich
    // ------------------------------------------------------------

    int nelsToEnrich = 0;

    for (auto localEl : cmeshL->ElementVec()) {

        if (!localEl) continue;
        if (!localEl->Reference()) continue;

        const bool isDomainElement =
            localEl->Dimension() == cmeshL->Dimension();

        const bool isNeumannBoundary =
            localEl->Reference()->Material() == overlappingNHNeumannBoundary;

        if (!isDomainElement && !isNeumannBoundary) continue;

        nelsToEnrich++;
    }

    const int nElementsDispG  = dispMeshG->NElements();
    const int nElementsMixedG = cmeshG->NElements();

    dispMeshG->ElementVec().resize(nElementsDispG + nelsToEnrich);
    cmeshG->ElementVec().resize(nElementsMixedG + nelsToEnrich);

    // ------------------------------------------------------------
    // 3. Create enriched elements
    // ------------------------------------------------------------

    count = 0;

    for (auto localEl : cmeshL->ElementVec()) {

        if (!localEl) continue;
        if (!localEl->Reference()) continue;

        const bool isDomainElement =
            localEl->Dimension() == cmeshL->Dimension();

        const bool isNeumannBoundary =
            localEl->Reference()->Material() == overlappingNHNeumannBoundary;

        if (!isDomainElement && !isNeumannBoundary) continue;

        auto itCorrespondence =
            globalElementCorrespondence.find(localEl->Index());

        if (itCorrespondence == globalElementCorrespondence.end()) {
            PanicButton();
            return;
        }

        const int globalElIndex = itCorrespondence->second;

        if (globalElIndex < 0 ||
            globalElIndex >= static_cast<int>(cmeshG->ElementVec().size()))
        {
            PanicButton();
            return;
        }

        Element *globalEl = cmeshG->ElementVec()[globalElIndex];

        if (!globalEl) {
            PanicButton();
            return;
        }

        ElementMixed *globalElMixed =
            dynamic_cast<ElementMixed *>(globalEl);

        if (!globalElMixed) {
            PanicButton();
            return;
        }

        if (globalElMixed->SubElements().size() < 2 ||
            !globalElMixed->SubElements()[0] ||
            !globalElMixed->SubElements()[1])
        {
            PanicButton();
            return;
        }

        const int newDispElementIndex  = nElementsDispG  + count;
        const int newMixedElementIndex = nElementsMixedG + count;

        ElementEnriched *enrichedEl =
            new ElementEnriched(
                newMixedElementIndex,
                localEl,
                globalEl,
                mixedCmeshG,
                globalLocal
            );

        enrichedEl->setCorrespondence(
            &globalElementCorrespondence,
            &globalNodeCorrespondence
        );

        enrichedEl->SetEnrichmentData(&enrichedConnects);

        // --------------------------------------------------------
        // Remove original global element weak form
        // --------------------------------------------------------

        if (globalEl->Dimension() == dispMeshG->Dimension()) {

            if (globalElIndex >= 0 &&
                globalElIndex < static_cast<int>(dispMeshG->ElementVec().size()) &&
                dispMeshG->ElementVec()[globalElIndex])
            {
                dispMeshG->ElementVec()[globalElIndex]->SetWeakForm(nullptr);
            }

            globalEl->SetWeakForm(nullptr);
        }

        // --------------------------------------------------------
        // Connectivity order:
        //
        // [standard displacement connects]
        // [enriched displacement connects]
        // [pressure connects]
        // --------------------------------------------------------

        std::vector<Connect *> enrichedConnectivity =
            globalElMixed->SubElements()[0]->getConnectivity();

        auto globalDispConnects =
            globalElMixed->SubElements()[0]->getConnectivity();

        for (auto connect : globalDispConnects) {

            if (!connect) continue;

            auto itEnriched =
                enrichedConnects.find(connect->Index());

            if (itEnriched == enrichedConnects.end()) continue;

            const int enrichedConnectIndex = itEnriched->second;

            if (enrichedConnectIndex < 0 ||
                enrichedConnectIndex >=
                    static_cast<int>(dispMeshG->ConnectVec().size()))
            {
                PanicButton();
                return;
            }

            Connect *enrichedConnect =
                dispMeshG->ConnectVec()[enrichedConnectIndex];

            if (!enrichedConnect) {
                PanicButton();
                return;
            }

            enrichedConnectivity.push_back(enrichedConnect);
        }

        std::vector<Connect *> pressureConnectivity =
            globalElMixed->SubElements()[1]->getConnectivity();

        enrichedConnectivity.insert(
            enrichedConnectivity.end(),
            pressureConnectivity.begin(),
            pressureConnectivity.end()
        );

        enrichedEl->setConnectivity(enrichedConnectivity);

        if (newDispElementIndex < 0 ||
            newDispElementIndex >=
                static_cast<int>(dispMeshG->ElementVec().size()))
        {
            PanicButton();
            return;
        }

        if (newMixedElementIndex < 0 ||
            newMixedElementIndex >=
                static_cast<int>(cmeshG->ElementVec().size()))
        {
            PanicButton();
            return;
        }

        dispMeshG->ElementVec()[newDispElementIndex] = enrichedEl;
        cmeshG->ElementVec()[newMixedElementIndex] = enrichedEl;

        count++;
    }

    // ------------------------------------------------------------
    // 4. Update total number of global DOFs
    // ------------------------------------------------------------

    int nGlobalDOF = 0;

    for (int imesh = 0;
         imesh < static_cast<int>(mixedCmeshG->MeshVector().size());
         imesh++)
    {
        CompMesh *mesh = mixedCmeshG->MeshVector()[imesh];

        if (!mesh) continue;

        for (int ic = 0; ic < mesh->NConnects(); ic++) {

            Connect *con = mesh->ConnectVec()[ic];

            if (!con) {
                PanicButton();
                return;
            }

            const int nshape = con->GetNShapeFunctions();
            const int nstate = con->GetNStateVariables();

            nGlobalDOF += nshape * nstate;
        }
    }

    cmeshG->NGlobalDOF() = nGlobalDOF;
}