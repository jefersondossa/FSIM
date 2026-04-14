#include "ElementEnriched.h"
#include "Elasticity2D.h"
#include "ElasticityPositional2D.h"
#include "ElasticTruss.h"
#include "PositionalTruss.h"
#include "GlobalLocalEnrichment.h"
#include "MixedGlobalLocalEnrichment.h"
#include "ElementMixed.h"

ElementEnriched::ElementEnriched(int64_t index, Element* localEl, Element* globalEl, CompMesh* cmesh, WeakForm *wf) : Element(){
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
    //int64_t elGlobalIndex = globalElementCorrespondence->at(fLocalElement->Index());
    MatrixDouble elGlobalXsi = globalNodeCorrespondence->at(fLocalElement->Index());
    // auto intrule = fLeftElement->GetIntRule();
    ElementMixed *mixed = dynamic_cast<ElementMixed*>(fGlobalElement);


    for(int it = 0; it < fLocalElement->getNumberOfIntegrationPoints(); it++){

        // //Defines the integration points adimentional coordinates 
        for (int k = 0; k < DIM; k++){
            double coord = fLocalElement->IntegPointCoordinate(index,k);
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
        if (mixed){
            for (int i = 0; i < mixed->SubElements().size(); i++){
            mixed->SubElements()[i]->IntegrationData().fA0 = fGlobalElement->IntegrationData().fA0;
            mixed->SubElements()[i]->IntegrationData().fA0Inv = fGlobalElement->IntegrationData().fA0Inv;
            mixed->SubElements()[i]->IntegrationData().fAxes0 = fGlobalElement->IntegrationData().fAxes0;
            mixed->SubElements()[i]->IntegrationData().fJacA0 = fGlobalElement->IntegrationData().fJacA0;
            mixed->SubElements()[i]->IntegrationData().fX = fGlobalElement->IntegrationData().fX;

            //Computes spatial derivatives
            mixed->SubElements()[i]->ComputeSpatialDerivatives();
            mixed->SubElements()[i]->interpolateSolution();
            mixed->SubElements()[i]->interpolateSolDerivatives();
            data[i] = &mixed->SubElements()[i]->IntegrationData();
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
        RhsCorrect[i] = rhsVector[i];
        for (int j = 0; j < nglobalDOF; j++){
            StiffnessCorrect(i,j) = jacobianNRMatrix(i,j);   
        }
    }

    //Account the enriched DOFs in the stiffness matrix and rhs vector
    std::set<int> locIndexes;
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
                break;
            }
        }
        locConnectIndex -= fGlobalElement->NSides();
        locIndexes.insert(locConnectIndex);
        int nshapei = c->GetNShapeFunctions();
        int nstatei = c->GetNStateVariables();
        for (int j = 0; j < nshapei*nstatei; j++){
            RhsCorrect[nglobalDOF + locConnectIndex*nstatei + j] = rhsVector[nglobalDOF + i*nstatei + j];
            //Global DOFs x Enriched DOFs
            for (int k = 0; k < nglobalDOF; k++){
                StiffnessCorrect(nglobalDOF + locConnectIndex*nstatei + j,k) = jacobianNRMatrix(nglobalDOF + i*nstatei + j, k);   
                StiffnessCorrect(k,nglobalDOF + locConnectIndex*nstatei + j) = jacobianNRMatrix(k,nglobalDOF + i*nstatei + j);   
            }
            for (int k = 0; k < nshapei*nstatei; k++){
                //Enriched DOFs x Enriched DOFs - Diagonal
                StiffnessCorrect(nglobalDOF + locConnectIndex*nstatei + j,nglobalDOF + locConnectIndex*nstatei + k) = jacobianNRMatrix(nglobalDOF + i*nstatei + j,nglobalDOF + i*nstatei + k);   
                //Enriched DOFs x Enriched DOFs - Off-diagonal
                for (auto locconnect:locIndexes){
                    if (locconnect == locConnectIndex) continue;
                    StiffnessCorrect(nglobalDOF + locConnectIndex*nstatei + j,nglobalDOF + locconnect*nstatei + k) = jacobianNRMatrix(nglobalDOF + i*nstatei + j,nglobalDOF + locconnect*nstatei + k);
                    StiffnessCorrect(nglobalDOF + locconnect*nstatei + k,nglobalDOF + locConnectIndex*nstatei + j) = jacobianNRMatrix(nglobalDOF + locconnect*nstatei + k,nglobalDOF + i*nstatei + j);
                }
            }
        }
    }
    // std::cout << "Stiffness before \n" << jacobianNRMatrix << '\n';
    // std::cout << "Rhs \n" << rhsVector << '\n';

    
    jacobianNRMatrix = StiffnessCorrect;
    rhsVector = RhsCorrect;

    // std::cout << "Stiffness \n" << jacobianNRMatrix << '\n';
    // std::cout << "Rhs \n" << rhsVector << '\n';

    return;
};