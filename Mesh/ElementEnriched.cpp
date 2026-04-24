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
    ElementMixed *mixedGlobal = dynamic_cast<ElementMixed*>(fGlobalElement);
    ElementMixed *mixedLocal = dynamic_cast<ElementMixed*>(fLocalElement);


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
                    StiffnessCorrect(nglobalDOF + locConnectIndex*nstatei + j,nglobalDOF + locconnect*nstatei + k) = Stiffness(nglobalDOF + i*nstatei + j,nglobalDOF + locconnect*nstatei + k);
                    StiffnessCorrect(nglobalDOF + locconnect*nstatei + k,nglobalDOF + locConnectIndex*nstatei + j) = Stiffness(nglobalDOF + locconnect*nstatei + k,nglobalDOF + i*nstatei + j);
                }
            }
        }
    }
    // std::cout << "Stiffness before \n" << jacobianNRMatrix << '\n';
    // std::cout << "Rhs \n" << rhsVector << '\n';

    
    Stiffness = StiffnessCorrect;
    Rhs = RhsCorrect;

}

void ElementEnriched::AccountForEnrichmentMixed(MatrixDouble &Stiffness, VecDouble &Rhs){

    if (this->Dimension()!= this->Mesh()->Dimension()) return;

    // std::cout << "Stiffness before \n" << Stiffness << '\n';

    //Resize the stiffnes matrix and vector to account only the enriched DOFs
    int nDOFcorrect=0;
    // int nglobalDOF = fGlobalElement->NLocDOF();
    ElementMixed *mixed = dynamic_cast<ElementMixed*>(fGlobalElement);
    if (!mixed){
        PanicButton();
        return;
    }
    int nglobalDOF = mixed->SubElements()[0]->NLocDOF();

    // for (int i = 0; i < mixed->SubElements().size(); i++){
        // auto connect = mixed->SubElements()[i]->getConnectivity();
        for (int j = 0; j < fConnect.size(); j++){
            Connect * c = fConnect[j];
            nDOFcorrect += c->GetNStateVariables() * c->GetNShapeFunctions();
        }
    // }
    
    

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

    int nstateEnriched = 0;
    //Account the enriched DOFs in the stiffness matrix and rhs vector
    std::set<int> locIndexes;
    auto globalConnects = fGlobalElement->getConnectivity();
    int nEnrichedConnects = this->getConnectivity().size() - globalConnects.size();
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
        nstateEnriched = nstatei;
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
                    // tem alguma coisa errada aqui, revisar depois
                    StiffnessCorrect(nglobalDOF + locConnectIndex*nstatei + j,nglobalDOF + locconnect*nstatei + k) = Stiffness(nglobalDOF + i*nstatei + j,nglobalDOF + locconnect*nstatei + k);

                    // StiffnessCorrect(nglobalDOF + locConnectIndex*nstatei + j,nglobalDOF + locconnect*nstatei + k) = Stiffness(nglobalDOF + i*nstatei + j,nglobalDOF + locconnect*nstatei + k);
                    StiffnessCorrect(nglobalDOF + locconnect*nstatei + k,nglobalDOF + locConnectIndex*nstatei + j) = Stiffness(nglobalDOF + i*nstatei + j,nglobalDOF + locconnect*nstatei + k);
                }
            }
        }
    }

    //Account pressure DOFs in the stiffness matrix and rhs vector
    int nPressureConnects = mixed->SubElements()[1]->getConnectivity().size();
    int pressureCount = this->getConnectivity().size() - nPressureConnects;
    for (int i = pressureCount; i < this->getConnectivity().size(); i++){
        Connect * c = fConnect[i];
        int nshapei = c->GetNShapeFunctions();
        int nstatei = c->GetNStateVariables();
        for (int j = 0; j < nshapei*nstatei; j++){
            RhsCorrect[nglobalDOF + locIndexes.size()*nstatei + (i-pressureCount)*nstatei + j] = Rhs[nglobalDOF + (pressureCount + i - pressureCount)*nstatei + j];
            //Global DOFs x Pressure DOFs
            for (int k = 0; k < nglobalDOF; k++){
                StiffnessCorrect(nglobalDOF + locIndexes.size()*nstateEnriched + (i-pressureCount)*nstatei + j,k) = Stiffness(nglobalDOF*nstateEnriched + (i-pressureCount)*nstatei + j, k);   
                StiffnessCorrect(k,nglobalDOF + locIndexes.size()*nstateEnriched + (i-pressureCount)*nstatei + j) = Stiffness(k, nglobalDOF*nstateEnriched + (i-pressureCount)*nstatei + j);   
            }
            //Enriched DOFs x Pressure DOFs
            // for (auto locconnect:locIndexes){
            //     for (int k = 0; k < nshapei*nstatei; k++){
            //         StiffnessCorrect(nglobalDOF + locconnect*nstatei + k, nglobalDOF + locIndexes.size()*nstatei + (i-pressureCount)*nstatei + j) = Stiffness(nglobalDOF + locconnect*nstatei + k, nglobalDOF*nstateEnriched + (i - pressureCount)*nstatei + j);
            //         StiffnessCorrect(nglobalDOF + locIndexes.size()*nstatei + (i-pressureCount)*nstatei + j, nglobalDOF + locconnect*nstatei + k) = Stiffness(nglobalDOF*nstateEnriched + (i - pressureCount)*nstatei + j, nglobalDOF + locconnect*nstatei + k);
            //     }
            // }
            //Pressure DOFs x Pressure DOFs
            for (int k = 0; k < nshapei*nstatei; k++){
                StiffnessCorrect(nglobalDOF + locIndexes.size()*nstateEnriched + (i-pressureCount)*nstatei + j, nglobalDOF + locIndexes.size()*nstateEnriched + (i-pressureCount)*nstatei + k) = Stiffness(nglobalDOF*nstateEnriched + (i - pressureCount)*nstatei + j, nglobalDOF*nstateEnriched + (i - pressureCount)*nstatei + k);
            }
        }
        
    }
    

    // std::cout << "Stiffness before \n" << jacobianNRMatrix << '\n';
    // std::cout << "Rhs \n" << rhsVector << '\n';

    // Stiffness.setZero(); Rhs.setZero();
    Stiffness = StiffnessCorrect;
    Rhs = RhsCorrect;


}