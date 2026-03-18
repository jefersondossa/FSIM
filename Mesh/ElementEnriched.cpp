#include "ElementEnriched.h"
#include "Elasticity2D.h"
#include "ElasticityPositional2D.h"
#include "ElasticTruss.h"
#include "PositionalTruss.h"
#include "GlobalLocalEnrichment.h"


ElementEnriched::ElementEnriched(int64_t index, Element* localEl, Element* globalEl, CompMesh* cmesh, WeakForm *wf) : Element(){
    this->fIndex = index;
    fLocalElement = localEl;
    fGlobalElement = globalEl;
    this->fWeakForm = wf;
};

Element *ElementEnriched::Clone() const {
    return new ElementEnriched(*this);
};

void ElementEnriched::ComputeElContribution(MatrixDouble &jacobianNRMatrix, VecDouble &rhsVector){

    if (!this->fWeakForm) return;

    GlobalLocalEnrichment *globalLocal = dynamic_cast<GlobalLocalEnrichment*>(this->fWeakForm);
    if (!globalLocal) {
        PanicButton();
        return;
    }

    int DIM = fLocalElement->Dimension();
    int index = 0;
    //int64_t elGlobalIndex = globalElementCorrespondence->at(fLocalElement->Index());
    MatrixDouble elGlobalXsi = globalNodeCorrespondence[fLocalElement->Index()];
    // auto intrule = fLeftElement->GetIntRule();
    
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
        fGlobalElement->ComputeSpatialDerivatives();

        if (fLocalElement->IntegrationData().fNeedsSol) fLocalElement->interpolateSolution();
        if (fLocalElement->IntegrationData().fNeedsDSol) fLocalElement->interpolateSolDerivatives();
        if (fGlobalElement->IntegrationData().fNeedsSol) fGlobalElement->interpolateSolution();
        if (fGlobalElement->IntegrationData().fNeedsDSol) fGlobalElement->interpolateSolDerivatives();

        //Computes the element diffusion/viscosity matrix
        globalLocal->ComputeStiffness(index, fLocalElement->IntegrationData(), fGlobalElement->IntegrationData(), jacobianNRMatrix);

        //Computes the RHS vector
        globalLocal->ComputeResidual(index, fLocalElement->IntegrationData(), fGlobalElement->IntegrationData(), rhsVector); 

        index++;        
    };  

    // std::cout << "Stiffness \n" << jacobianNRMatrix << '\n';
    // std::cout << "Rhs \n" << rhsVector << '\n';

    return;
};