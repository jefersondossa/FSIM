#include "InterpolatedBC.h"
#include "ElementMixed.h"

InterpolatedBC::InterpolatedBC(int matid, int dim, int nstate, 
        BoundaryConditionType bctype, 
        std::map<int64_t,int64_t> *globalElementCorrespondence,
        std::map<int64_t, MatrixDouble> *globalNodeCorrespondence,
        CompMesh *globalmesh) : WeakForm() {
    this->fDimension = dim;
    this->fMatId = matid;
    fNState = nstate;
    BCType = bctype;
    fGlobalElementCorrespondence = globalElementCorrespondence;
    fGlobalNodeCorrespondence = globalNodeCorrespondence;
    fGlobalMesh = globalmesh;
};

void InterpolatedBC::ComputeStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness){
    if (!data.fNeedsSol){
        data.fNeedsSol = true;
        data.fSol.resize(fNState);
    }
    double WJ = data.fWeight * data.fJacA0 * WeakForm::fBigNumber;
    double nphi = data.fPhi.size();
    
    switch (BCType)
    {
    case BoundaryConditionType::kDirichlet: // Dirichlet
        for (int i = 0; i < nphi; i++){
            for (int j = 0; j < nphi; j++){
                for (int istate = 0; istate < fNState; istate++){
                    Stiffness(fNState*i+istate,fNState*j+istate) +=  WJ * data.fPhi[i] * data.fPhi[j];
                }
            }
        }
        break;
    case BoundaryConditionType::kNeumann:
        break;
    
    default:
        std::cout << "BC Type not implemented \n" ;
        PanicButton();
        break;
    }
   
}

void InterpolatedBC::ComputeResidual(int &index, IntPointData &data, VecDouble &Rhs){

    double WJ = data.fWeight * data.fJacA0;
    double nphi = data.fPhi.size();

    int64_t myIndex = data.fElementIndex;
    Element *el = fGlobalMesh->ElementVec()[(*fGlobalElementCorrespondence)[myIndex]];
    int nshape = 0;
    VecDouble Sol;

    ElementMixed *elMixed = dynamic_cast<ElementMixed *>(el);
    if (elMixed){
        elMixed->SubElements()[0]->IntegrationData().fAdimCoord[0] = (*fGlobalNodeCorrespondence)[myIndex](index,0); 
        elMixed->SubElements()[0]->IntegrationData().fAdimCoord[1] = (*fGlobalNodeCorrespondence)[myIndex](index,1); 
        
        nshape = elMixed->SubElements()[0]->NShapeFunctions();
        elMixed->SubElements()[0]->IntegrationData().fPhi.resize(nshape);
        elMixed->SubElements()[0]->IntegrationData().fDPhi.resize(elMixed->SubElements()[0]->Dimension(),nshape);    
        elMixed->SubElements()[0]->Reference()->ComputeJacobian(elMixed->SubElements()[0]->IntegrationData());
        elMixed->SubElements()[0]->ComputeSpatialDerivatives();
        elMixed->SubElements()[0]->interpolateSolution();
        Sol = elMixed->SubElements()[0]->IntegrationData().fSol;
    }else{
        el->IntegrationData().fAdimCoord[0] = (*fGlobalNodeCorrespondence)[myIndex](index,0); 
        el->IntegrationData().fAdimCoord[1] = (*fGlobalNodeCorrespondence)[myIndex](index,1); 
        
        nshape = el->NShapeFunctions();
        el->IntegrationData().fPhi.resize(nshape);
        el->IntegrationData().fDPhi.resize(el->Dimension(),nshape);    
        el->Reference()->ComputeJacobian(el->IntegrationData());
        el->ComputeSpatialDerivatives();
        el->interpolateSolution();
        Sol = el->IntegrationData().fSol;
    }


    switch (BCType)
    {
    case BoundaryConditionType::kDirichlet: // Dirichlet in all state variables
        for (int i = 0; i < nphi; i++){
            for (int istate = 0; istate < fNState; istate++){
                Rhs(fNState*i+istate) +=  WeakForm::fBigNumber * WJ * data.fPhi[i] * Sol[istate];
            }
        }
        break;
    case BoundaryConditionType::kNeumann: // Neumann in all state variables
    {
        for (int i = 0; i < nphi; i++){
            for (int istate = 0; istate < fNState; istate++){
                Rhs(fNState*i+istate) +=  WJ * data.fPhi[i] * (Sol[istate]);
            }
        }
    }
        break;
    default:
        std::cout << "BC Type not implemented \n" ;
        PanicButton();
        break;
    }
};

void InterpolatedBC::ComputeError(IntPointData &data, VecDouble &errors){

}



int InterpolatedBC::VariableIndex(const std::string &name) const{
    
    if(!strcmp("Solution",name.c_str()))        return 1;
    if(!strcmp("DerivativeX",name.c_str()))      return 2;
    if(!strcmp("DerivativeY",name.c_str()))      return 3;
    if(!strcmp("DerivativeZ",name.c_str()))      return 4;
    if(!strcmp("Material",name.c_str()))         return 100;
    
    std::cout << "Post Process variable not implemented \n";
    PanicButton();
    return -1;
};

int InterpolatedBC::NSolutionVariables(int var) const{
    switch (var)
    {
    case 1:
    case 2:
    case 3:
    case 4:
    case 100:
        return 1;

    default:
        PanicButton();
        return -1;
    }
};

void InterpolatedBC::Solution(IntPointData &data, int var, VecDouble &Sol) {

    if (var == 100){
        Sol[0] = fMatId;
        return;
    };

    //Solution
    if (var == 1){
        for (int i = 0; i < fNState; i++){
            Sol[i] = data.fSol[i];
        }
        return;
    };

    //Derivative X
    if (var == 2){
        for (int i = 0; i < fNState; i++){
            Sol[i] = data.fDSolDx(0,i);
        }
        return;
    };
    //Derivative Y
    if (var == 3){
        for (int i = 0; i < fNState; i++){
            Sol[i] = data.fDSolDx(1,i);
        }
        return;
    };
    //Derivative Z
    if (var == 4){
        for (int i = 0; i < fNState; i++){
            Sol[i] = data.fDSolDx(2,i);
        }
        return;
    };

}; 