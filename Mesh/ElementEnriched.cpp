#include "ElementEnriched.h"
#include "Elasticity2D.h"
#include "ElasticityPositional2D.h"
#include "ElasticTruss.h"
#include "PositionalTruss.h"
#include "GlobalLocalEnrichment.h"
#include "MixedGlobalLocalEnrichment.h"
#include "ElementMixed.h"

ElementEnriched::ElementEnriched(int index, Element* localEl, Element* globalEl, CompMesh* cmesh, WeakForm *wf) : Element(){
    this->fIndex = index;
    fMesh = cmesh;
    fLocalElement = localEl;
    fGlobalElement = globalEl;
    this->fWeakForm = wf;
    ElementMixed *mixed = dynamic_cast<ElementMixed*>(globalEl);
    if (mixed){
        this->nLocDOF = mixed->SubElements()[0]->NLocDOF()*2 + mixed->SubElements()[1]->NLocDOF();
    } else {
        this->nLocDOF = globalEl->NLocDOF() * 2;
    }
    this->fReference = localEl->Reference();
    fLocalElement->IntegrationData().fNeedsSol = true;
    fLocalElement->IntegrationData().fNeedsDSol = true;
    fGlobalElement->IntegrationData().fNeedsSol = true;
    fGlobalElement->IntegrationData().fNeedsDSol = true;
};

Element *ElementEnriched::Clone() const {
    return new ElementEnriched(*this);
};

void ElementEnriched::ComputeElContribution(MatrixDouble &jacobianNRMatrix, VecDouble &rhsVector){

    if (!this->fWeakForm) return;

    GlobalLocalEnrichment *globalLocal = dynamic_cast<GlobalLocalEnrichment*>(this->fWeakForm);
    MixedGlobalLocalEnrichment *globalLocalMixed = dynamic_cast<MixedGlobalLocalEnrichment*>(this->fWeakForm);
    
    if (!globalLocal && !globalLocalMixed){
        PanicButton();
        return;
    }
    
    int DIM = fLocalElement->Dimension();
    int index = 0;
    //int elGlobalIndex = globalElementCorrespondence->at(fLocalElement->Index());
    MatrixDouble elGlobalXsi = globalNodeCorrespondence->at(fLocalElement->Index());
    // auto intrule = fLeftElement->GetIntRule();
    ElementMixed *mixedGlobal = dynamic_cast<ElementMixed*>(fGlobalElement);
    ElementMixed *mixedLocal = dynamic_cast<ElementMixed*>(fLocalElement);

    // fGlobalElement->IntegrationData().fJacA1 = 0;

    for(int it = 0; it < fLocalElement->getNumberOfIntegrationPoints(); it++){

        // //Defines the integration points adimentional coordinates 
        for (int k = 0; k < DIM; k++){
            REAL coord = fLocalElement->IntegPointCoordinate(index,k);
            fLocalElement->IntegrationData().fAdimCoord[k] = coord;
            fGlobalElement->IntegrationData().fAdimCoord[k] = elGlobalXsi(index,k); 
        }         

        // //Returns the quadrature integration weight
        fLocalElement->IntegrationData().fWeight = fLocalElement->IntegPointWeight(index);
        fGlobalElement->IntegrationData().fWeight = fLocalElement->IntegPointWeight(index);

        //Computes the jacobian matrix
        fLocalElement->Reference()->ComputeJacobian(fLocalElement->IntegrationData());
        fGlobalElement->Reference()->ComputeJacobian(fGlobalElement->IntegrationData());

        //Computes spatial derivatives
        fLocalElement->ComputeSpatialDerivatives();
        std::vector<IntPointData *> data(2);
        if (mixedGlobal){
            for (int i = 0; i < mixedGlobal->SubElements().size(); i++){
                mixedGlobal->SubElements()[i]->IntegrationData().fAdimCoord = fGlobalElement->IntegrationData().fAdimCoord;
                mixedGlobal->SubElements()[i]->IntegrationData().fWeight = fGlobalElement->IntegrationData().fWeight;
                mixedGlobal->SubElements()[i]->IntegrationData().fA0 = fGlobalElement->IntegrationData().fA0;
                mixedGlobal->SubElements()[i]->IntegrationData().fA0Inv = fGlobalElement->IntegrationData().fA0Inv;
                mixedGlobal->SubElements()[i]->IntegrationData().fAxes0 = fGlobalElement->IntegrationData().fAxes0;
                mixedGlobal->SubElements()[i]->IntegrationData().fJacA0 = fGlobalElement->IntegrationData().fJacA0;
                mixedGlobal->SubElements()[i]->IntegrationData().fX = fGlobalElement->IntegrationData().fX;

                //Computes spatial derivatives
                mixedGlobal->SubElements()[i]->ComputeSpatialDerivatives();
                mixedGlobal->SubElements()[i]->interpolateSolution();
                mixedGlobal->SubElements()[i]->interpolateSolDerivatives();
                data[i] = &mixedGlobal->SubElements()[i]->IntegrationData();
            }
            for (int i = 0; i < mixedLocal->SubElements().size(); i++){
                mixedLocal->SubElements()[i]->IntegrationData().fAdimCoord = fLocalElement->IntegrationData().fAdimCoord;
                mixedLocal->SubElements()[i]->IntegrationData().fWeight = fLocalElement->IntegrationData().fWeight;
                mixedLocal->SubElements()[i]->IntegrationData().fA0 = fLocalElement->IntegrationData().fA0;
                mixedLocal->SubElements()[i]->IntegrationData().fA0Inv = fLocalElement->IntegrationData().fA0Inv;
                mixedLocal->SubElements()[i]->IntegrationData().fAxes0 = fLocalElement->IntegrationData().fAxes0;
                mixedLocal->SubElements()[i]->IntegrationData().fJacA0 = fLocalElement->IntegrationData().fJacA0;
                mixedLocal->SubElements()[i]->IntegrationData().fX = fLocalElement->IntegrationData().fX;

                //Computes spatial derivatives
                mixedLocal->SubElements()[i]->ComputeSpatialDerivatives();
                mixedLocal->SubElements()[i]->interpolateSolution();
                mixedLocal->SubElements()[i]->interpolateSolDerivatives();
            }
            
        } else {
            fGlobalElement->ComputeSpatialDerivatives();
            if (fGlobalElement->IntegrationData().fNeedsSol) fGlobalElement->interpolateSolution();
            if (fGlobalElement->IntegrationData().fNeedsDSol) fGlobalElement->interpolateSolDerivatives();
        }
        

        if (fLocalElement->IntegrationData().fNeedsSol) fLocalElement->interpolateSolution();
        if (fLocalElement->IntegrationData().fNeedsDSol) fLocalElement->interpolateSolDerivatives();
        
        //Computes the element diffusion/viscosity matrix
        if (globalLocal) globalLocal->ComputeStiffness(index, fLocalElement->IntegrationData(), fGlobalElement->IntegrationData(), jacobianNRMatrix);
        if (globalLocalMixed) globalLocalMixed->ComputeStiffness(index, fLocalElement->IntegrationData(), data, jacobianNRMatrix);

        //Computes the RHS vector
        if (globalLocal) globalLocal->ComputeResidual(index, fLocalElement->IntegrationData(), fGlobalElement->IntegrationData(), rhsVector); 
        if (globalLocalMixed) globalLocalMixed->ComputeResidual(index, fLocalElement->IntegrationData(), data, rhsVector); 

        index++;        
    };

    if (globalLocal) AccountForEnrichment(jacobianNRMatrix, rhsVector);
    if (globalLocalMixed) AccountForEnrichmentMixed(jacobianNRMatrix, rhsVector);

    if (globalLocal){
        std::cout << "Global index: " << fGlobalElement->Index() << ", Local index: " << fLocalElement->Index() << '\n';
        std::cout << "Área global: " << fGlobalElement->IntegrationData().fJacA0 << '\n';
        std::cout << "Área global/local: " << fGlobalElement->IntegrationData().fJacA1 << '\n';
    }
    if (globalLocalMixed){
        std::cout << "Global index: " << fGlobalElement->Index() << ", Local index: " << fLocalElement->Index() << '\n';
        std::cout << "Área global: " << mixedGlobal->SubElements()[0]->IntegrationData().fJacA0 << '\n';
        std::cout << "Área global/local: " << mixedGlobal->SubElements()[0]->IntegrationData().fJacA1 << '\n';
    }
   
    //PrintMathematica(jacobianNRMatrix, "Stiffness");
    // std::cout << "Stiffness \n" << jacobianNRMatrix << '\n';
    // std::cout << "Rhs \n" << rhsVector << '\n';

    return;
};


void ElementEnriched::AccountForEnrichment(MatrixDouble &Stiffness, VecDouble &Rhs){

    //Resize the stiffnes matrix and vector to account only the enriched DOFs
    int nDOFcorrect=0;
    int nglobalDOF = fGlobalElement->NLocDOF();
    for (int i = 0; i < fConnect.size(); i++){
        Connect *c = fConnect[i];
        nDOFcorrect += c->GetNStateVariables() * c->GetNShapeFunctions();
    }

    MatrixDouble StiffnessCorrect(nDOFcorrect, nDOFcorrect);
    VecDouble RhsCorrect(nDOFcorrect);
    StiffnessCorrect.setZero(); RhsCorrect.setZero();
    //Account the standard DOFs in the stiffness matrix and rhs vector
    for (int i = 0; i < nglobalDOF; i++){
        RhsCorrect[i] = Rhs[i];
        for (int j = 0; j < nglobalDOF; j++){
            StiffnessCorrect(i,j) = Stiffness(i,j);   
        }
    }

    //Account the enriched DOFs in the stiffness matrix and rhs vector
    std::set<int> locIndexes;
    std::map<int,int> locIndToEnrichIndex;
    auto globalConnects = fGlobalElement->getConnectivity();
    for (size_t i = 0; i < globalConnects.size(); i++){
        if (!(*connectEnrichment)[globalConnects[i]->Index()])continue;
        int enrichConnectIndex = (*connectEnrichment)[globalConnects[i]->Index()];
        Connect *c;
        int locConnectIndex = -1;
        for(int ic = 0; ic < fConnect.size(); ic++) {
            if(fConnect[ic]->Index() == enrichConnectIndex) {
                c = fConnect[ic];
                locConnectIndex = ic;
                locIndToEnrichIndex[locConnectIndex-fGlobalElement->NSides()] = i;
                break;
            }
        }
        locConnectIndex -= fGlobalElement->NSides();
        locIndexes.insert(locConnectIndex);
        int nshapei = c->GetNShapeFunctions();
        int nstatei = c->GetNStateVariables();
        for (int j = 0; j < nshapei*nstatei; j++){
            RhsCorrect[nglobalDOF + locConnectIndex*nstatei + j] = Rhs[nglobalDOF + i*nstatei + j];
            //Global DOFs x Enriched DOFs
            for (int k = 0; k < nglobalDOF; k++){
                StiffnessCorrect(nglobalDOF + locConnectIndex*nstatei + j,k) = Stiffness(nglobalDOF + i*nstatei + j, k);   
                StiffnessCorrect(k,nglobalDOF + locConnectIndex*nstatei + j) = Stiffness(k,nglobalDOF + i*nstatei + j);   
            }
            for (int k = 0; k < nshapei*nstatei; k++){
                //Enriched DOFs x Enriched DOFs - Diagonal
                StiffnessCorrect(nglobalDOF + locConnectIndex*nstatei + j,nglobalDOF + locConnectIndex*nstatei + k) = Stiffness(nglobalDOF + i*nstatei + j,nglobalDOF + i*nstatei + k);   
                //Enriched DOFs x Enriched DOFs - Off-diagonal
                for (auto locconnect:locIndexes){
                    if (locconnect == locConnectIndex) continue;
                    StiffnessCorrect(nglobalDOF + locConnectIndex*nstatei + j,nglobalDOF + locconnect*nstatei + k) = Stiffness(nglobalDOF + i*nstatei + j,nglobalDOF + locIndToEnrichIndex[locconnect]*nstatei + k);
                    StiffnessCorrect(nglobalDOF + locconnect*nstatei + k,nglobalDOF + locConnectIndex*nstatei + j) = Stiffness(nglobalDOF + locIndToEnrichIndex[locconnect]*nstatei + k,nglobalDOF + i*nstatei + j);
                }
            }
        }
    }
    // std::cout << "Stiffness before \n" << jacobianNRMatrix << '\n';
    // std::cout << "Rhs \n" << rhsVector << '\n';

    
    Stiffness = StiffnessCorrect;
    Rhs = RhsCorrect;

}


void ElementEnriched::AccountForEnrichmentMixed(MatrixDouble &Stiffness,
                                                VecDouble &Rhs)
{
    if (this->Dimension() != this->Mesh()->Dimension()) return;

    ElementMixed *mixed = dynamic_cast<ElementMixed*>(fGlobalElement);
    if (!mixed) {
        PanicButton();
        return;
    }

    if (Stiffness.rows() != Stiffness.cols()) {
        PanicButton();
        return;
    }

    if (Rhs.size() != Stiffness.rows()) {
        PanicButton();
        return;
    }

    auto globalConnects = fGlobalElement->getConnectivity();

    const int nGlobalConnects = static_cast<int>(globalConnects.size());
    const int nPressureConnects =
        static_cast<int>(mixed->SubElements()[1]->getConnectivity().size());

    if (nPressureConnects > nGlobalConnects) {
        PanicButton();
        return;
    }

    const int nDispGlobalConnects = nGlobalConnects - nPressureConnects;

    const int nglobalDOF = mixed->SubElements()[0]->NLocDOF();

    // ------------------------------------------------------------
    // Original assumed layout:
    //
    // [ standard displacement DOFs ]
    // [ candidate enriched displacement DOFs ]
    // [ pressure DOFs ]
    // ------------------------------------------------------------

    std::vector<int> keptOldDofs;

    // ------------------------------------------------------------
    // 1. Keep standard displacement DOFs
    // ------------------------------------------------------------
    for (int idof = 0; idof < nglobalDOF; idof++) {
        keptOldDofs.push_back(idof);
    }

    // ------------------------------------------------------------
    // 2. Build accumulated offsets for candidate enriched DOFs
    // ------------------------------------------------------------
    std::vector<int> oldEnrichedStart(nDispGlobalConnects + 1, 0);
    oldEnrichedStart[0] = nglobalDOF;

    for (int ic = 0; ic < nDispGlobalConnects; ic++) {

        Connect *c = globalConnects[ic];

        const int ndof =
            c->GetNShapeFunctions() *
            c->GetNStateVariables();

        oldEnrichedStart[ic + 1] = oldEnrichedStart[ic] + ndof;
    }

    const int oldPressureOffset = oldEnrichedStart[nDispGlobalConnects];

    // ------------------------------------------------------------
    // 3. Identify existing enriched connects in this element
    // ------------------------------------------------------------
    struct EnrichedInfo {
        int localConnectIndex;
        int globalDispConnectIndex;
        int oldStart;
        int ndof;
    };

    std::vector<EnrichedInfo> enrichedInfos;

    for (int ig = 0; ig < nDispGlobalConnects; ig++) {

        const int globalConnectIndex = globalConnects[ig]->Index();

        auto it = connectEnrichment->find(globalConnectIndex);
        if (it == connectEnrichment->end()) continue;

        const int enrichedConnectIndex = it->second;

        Connect *localConnect = nullptr;
        int localConnectIndex = -1;

        for (int iloc = 0; iloc < static_cast<int>(fConnect.size()); iloc++) {
            if (fConnect[iloc]->Index() == enrichedConnectIndex) {
                localConnect = fConnect[iloc];
                localConnectIndex = iloc;
                break;
            }
        }

        if (!localConnect || localConnectIndex < 0) {
            PanicButton();
            return;
        }

        const int ndof =
            localConnect->GetNShapeFunctions() *
            localConnect->GetNStateVariables();

        const int oldStart = oldEnrichedStart[ig];

        enrichedInfos.push_back({
            localConnectIndex,
            ig,
            oldStart,
            ndof
        });
    }

    std::sort(
        enrichedInfos.begin(),
        enrichedInfos.end(),
        [](const EnrichedInfo &a, const EnrichedInfo &b) {
            return a.localConnectIndex < b.localConnectIndex;
        }
    );

    // ------------------------------------------------------------
    // 4. Keep enriched DOFs that actually exist in this element
    // ------------------------------------------------------------
    for (const auto &info : enrichedInfos) {

        for (int idof = 0; idof < info.ndof; idof++) {

            const int oldDof = info.oldStart + idof;

            if (oldDof < 0 || oldDof >= Stiffness.rows()) {
                PanicButton();
                return;
            }

            keptOldDofs.push_back(oldDof);
        }
    }

    // ------------------------------------------------------------
    // 5. Keep pressure DOFs
    // ------------------------------------------------------------
    const int pressureStartConnect =
        static_cast<int>(fConnect.size()) - nPressureConnects;

    if (pressureStartConnect < 0) {
        PanicButton();
        return;
    }

    int pressureOffset = 0;

    for (int ic = pressureStartConnect;
         ic < static_cast<int>(fConnect.size());
         ic++)
    {
        Connect *c = fConnect[ic];

        const int ndof =
            c->GetNShapeFunctions() *
            c->GetNStateVariables();

        for (int idof = 0; idof < ndof; idof++) {

            const int oldDof =
                oldPressureOffset +
                pressureOffset +
                idof;

            if (oldDof < 0 || oldDof >= Stiffness.rows()) {
                PanicButton();
                return;
            }

            keptOldDofs.push_back(oldDof);
        }

        pressureOffset += ndof;
    }

    // ------------------------------------------------------------
    // 6. Build old-to-new DOF map
    // ------------------------------------------------------------
    std::vector<int> oldToNew(Stiffness.rows(), -1);

    for (int newDof = 0; newDof < static_cast<int>(keptOldDofs.size()); newDof++) {

        const int oldDof = keptOldDofs[newDof];

        if (oldDof < 0 || oldDof >= Stiffness.rows()) {
            PanicButton();
            return;
        }

        if (oldToNew[oldDof] != -1) {
            // Same old DOF inserted twice
            PanicButton();
            return;
        }

        oldToNew[oldDof] = newDof;
    }

    const int nCorrectDOF = static_cast<int>(keptOldDofs.size());

    MatrixDouble StiffnessCorrect(nCorrectDOF, nCorrectDOF);
    VecDouble RhsCorrect(nCorrectDOF);

    StiffnessCorrect.setZero();
    RhsCorrect.setZero();

    // ------------------------------------------------------------
    // 7. Copy RHS and stiffness matrix using the map
    // ------------------------------------------------------------
    for (int iold = 0; iold < Stiffness.rows(); iold++) {

        const int inew = oldToNew[iold];

        if (inew < 0) continue;

        RhsCorrect[inew] = Rhs[iold];

        for (int jold = 0; jold < Stiffness.cols(); jold++) {

            const int jnew = oldToNew[jold];

            if (jnew < 0) continue;

            StiffnessCorrect(inew, jnew) = Stiffness(iold, jold);
        }
    }

    Stiffness = StiffnessCorrect;
    Rhs = RhsCorrect;
}