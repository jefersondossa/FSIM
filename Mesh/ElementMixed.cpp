#include "ElementMixed.h"

ElementMixed::ElementMixed(int index, std::vector<Element *> elvector, MixedCompMesh* mesh, WeakForm *wf){
    int DIM = elvector[0]->Dimension();
    fSubElements = elvector;
    this->fIndex = index;
    this->fReference = elvector[0]->Reference();
    this->fMesh = mesh;
    this->fWeakForm = wf;

    fIntegData.fAdimCoord.resize(DIM);

    fIntegData.fA0Inv.resize(DIM,DIM);
    fIntegData.fA0Inv.setZero();
    fIntegData.fAxes0.resize(3,DIM);
    fIntegData.fAxes0.setZero();
    fIntegData.fA0.resize(DIM,DIM);
    fIntegData.fA0.setZero();
    fIntegData.fX.resize(3);
    fIntegData.fX.setZero();
    
    this->nLocDOF = 0;

    for (int i = 0; i < fSubElements.size(); i++){
        int locdof = fSubElements[i]->NLocDOF();
        this->nLocDOF += locdof;
    }
    
};

Element *ElementMixed::Clone() const {
    return new ElementMixed(*this);
};

void ElementMixed::ComputeElContribution(MatrixDouble &jacobianNRMatrix, VecDouble &rhsVector){

    if (!this->fWeakForm && !fSubElements[0]->GetWeakForm()) return;

    int DIM = this->fReference->Dimension();
    int nstate = 0;
    if (fSubElements[0]->GetWeakForm()){
        nstate = fSubElements[0]->GetWeakForm()->NState();
    } else {
        nstate = this->fWeakForm->NState();
    }

    for (int i = 0; i < fSubElements.size(); i++){
        auto &integdata = fSubElements[i]->IntegrationData();
        integdata.fNeedsDSol = true;
        if (fSubElements[i]->GetWeakForm()){
            nstate = fSubElements[i]->GetWeakForm()->NState();
        }
        integdata.fDSolDx.resize(nstate, DIM);
        integdata.fNeedsSol = true;
        integdata.fSol.resize(nstate);
    }
    
    int index = 0;
    
    int nsub = fSubElements.size();
    for (int i = 0; i < nsub; i++){
        fSubElements[i]->IntegrationData().fElementIndex = this->fIndex;
    }
    
    //The integration is always performed with basis in the first sub element.
    for(int it = 0; it < fSubElements[0]->getNumberOfIntegrationPoints(); it++){

        //Defines the integration points adimentional coordinates
        for (int k = 0; k < DIM; k++) {
            REAL coord = fSubElements[0]->IntegPointCoordinate(index,k);
            this->fIntegData.fAdimCoord[k] = coord;
            for (int i = 0; i < nsub; i++){
                fSubElements[i]->IntegrationData().fAdimCoord[k] = coord;
            }
        }

        //Returns the quadrature integration weight
        REAL weight = fSubElements[0]->IntegPointWeight(index);
        this->fIntegData.fWeight = weight;
        for (int i = 0; i < nsub; i++){
            fSubElements[i]->IntegrationData().fWeight = weight;
        }
            
        //Computes the jacobian matrix
        this->fReference->ComputeJacobian(this->fIntegData);
        for (int i = 0; i < nsub; i++){
            fSubElements[i]->IntegrationData().fA0 = this->fIntegData.fA0;
            fSubElements[i]->IntegrationData().fA0Inv = this->fIntegData.fA0Inv;
            fSubElements[i]->IntegrationData().fAxes0 = this->fIntegData.fAxes0;
            fSubElements[i]->IntegrationData().fJacA0 = this->fIntegData.fJacA0;
            fSubElements[i]->IntegrationData().fX = this->fIntegData.fX;

            //Computes spatial derivatives
            fSubElements[i]->ComputeSpatialDerivatives();
            fSubElements[i]->interpolateSolution();
            fSubElements[i]->interpolateSolDerivatives();
        }


        //Computes the element stiffness matrix and residual vector
        if (!this->fWeakForm){
            for (int i = 0; i < nsub; i++){
                fSubElements[i]->GetWeakForm()->ComputeStiffness(index, fSubElements[i]->IntegrationData(), jacobianNRMatrix);
                fSubElements[i]->GetWeakForm()->ComputeResidual(index, fSubElements[i]->IntegrationData(), rhsVector);
            }
            
        } else {
            std::vector<IntPointData *> data(nsub);
            for (int i = 0; i < nsub; i++){
                data[i] = &fSubElements[i]->IntegrationData();
            }
            this->fWeakForm->ComputeStiffness(index, data, jacobianNRMatrix);
            this->fWeakForm->ComputeResidual(index, data, rhsVector); 
        }
        
        index++;        
    };  

    // std::cout << "Stiffness \n" << jacobianNRMatrix << std::endl;
    // std::cout << "Rhs \n" << rhsVector << std::endl;

    return;
};


//------------------------------------------------------------------------------
//-----------------------TRANSIENT NAVIER-STOKES PROBEM-------------------------
//------------------------------------------------------------------------------
void ElementMixed::ComputeElContribution(std::vector<MatrixDouble> &jacobianNRMatrix, std::vector<VecDouble> &rhsVector){
    PanicButton();
    if (!this->fWeakForm) return;

    int DIM = this->fReference->Dimension();
    
    this->fIntegData.fA0Inv.resize(DIM,DIM);
    this->fIntegData.fAdimCoord.resize(DIM);

    int index = 0;

    // for(int it = 0; it < this->fIntRule.NPoints(); it++){

    //     //Defines the integration points adimentional coordinates
    //     for (int k = 0; k < DIM; k++) this->fIntegData.fAdimCoord[k] = this->fIntRule.PointList(index,k);

    //     //Returns the quadrature integration weight
    //     this->fIntegData.fWeight = this->fIntRule.WeightList(index);

    //     //Computes the jacobian matrix
    //     this->fReference->ComputeJacobian(this->fIntegData);

    //     //Computes spatial derivatives
    //     this->ComputeSpatialDerivatives();
        
    //     //Computes the element diffusion/viscosity matrix
    //     this->fWeakForm->ComputeStiffness(index, this->fIntegData, jacobianNRMatrix);

    //     if (this->fIntegData.fNeedsSol) this->interpolateSolution();
    //     if (this->fIntegData.fNeedsDSol) this->interpolateSolDerivatives();
    //     if (this->fIntegData.fNeedsTimeDerivatives) this->interpolateSolDTimeDerivatives();
        
    //     //Computes the RHS vector
    //     this->fWeakForm->ComputeResidual(index, this->fIntegData, rhsVector); 

    //     index++;        
    // };  
    // // std::cout << "\nStiffness Element " << this->Index() << "\n" << jacobianNRMatrix[1];
    // // std::cout << "\nrhsVector Element " << this->Index() << "\n" << rhsVector[1];

    return;
};

void ElementMixed::ComputeError(VecDouble &errors){
    if (!this->fWeakForm) return;

    int DIM = this->fReference->Dimension();
    
    this->fIntegData.fA0Inv.resize(DIM,DIM);
    this->fIntegData.fAdimCoord.resize(DIM);

    int index = 0;
    int nsub = fSubElements.size();

    for(int it = 0; it < fSubElements[0]->getNumberOfIntegrationPoints(); it++){

        //Defines the integration points adimentional coordinates
        for (int k = 0; k < DIM; k++) this->fIntegData.fAdimCoord[k] = fSubElements[0]->IntegPointCoordinate(index,k);

        //Returns the quadrature integration weight
        this->fIntegData.fWeight = fSubElements[0]->IntegPointWeight(index);

        //Computes the jacobian matrix
        //Computes the jacobian matrix
        this->fReference->ComputeJacobian(this->fIntegData);
        for (int i = 0; i < nsub; i++){
            fSubElements[i]->IntegrationData().fA0 = this->fIntegData.fA0;
            fSubElements[i]->IntegrationData().fA0Inv = this->fIntegData.fA0Inv;
            fSubElements[i]->IntegrationData().fAxes0 = this->fIntegData.fAxes0;
            fSubElements[i]->IntegrationData().fJacA0 = this->fIntegData.fJacA0;
            fSubElements[i]->IntegrationData().fX = this->fIntegData.fX;

            //Computes spatial derivatives
            fSubElements[i]->ComputeSpatialDerivatives();
            fSubElements[i]->interpolateSolution();
            fSubElements[i]->interpolateSolDerivatives();
        }
        
        //Computes the element error
        std::vector<IntPointData *> data(fSubElements.size());
        for (int i = 0; i < fSubElements.size(); i++){
            data[i] = &fSubElements[i]->IntegrationData();
        }
        this->fWeakForm->ComputeError(data, errors); 

        index++;        
    };
    // std::cout << "\nErrors Element " << this->Index() << "\n" << errors;
    return;
};