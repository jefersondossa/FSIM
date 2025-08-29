#include "CompMesh.h"
#include "Connect.h"
#include "TransientWeakForm.h"

GraphMesh* CompMesh::GetGraphMesh(){
    if (!fGraphMesh){
        fGraphMesh = new GraphMesh(this);
    }
    return fGraphMesh;
}

void CompMesh::InsertMaterial(WeakForm *wf){
    fMaterialVector[wf->Id()] = wf;
    if (fNState == 0){
        fNState = wf->NState();            
    } else if (fNState != wf->NState()){
        std::cout << "In the current vertion, the code don't support materials in the \
                      same mesh with different number of state variables. Please check it. \n";
        PanicButton(); 
    }
}

void CompMesh::Integrate(std::set<int> &matIds, std::vector<std::string> &varNames, std::map<std::string,VecDouble> &result){

    for (int64_t iel = 0; iel < NElements(); iel++){
        int elMatid = fElementVector[iel]->GetWeakForm()->Id();
        const bool is_in = matIds.find(elMatid) != matIds.end();
        if (!is_in) continue;
        
        fElementVector[iel]->Integrate(varNames,result);
    }
    

}

void CompMesh::SetSolution(VecDouble &sol){
    if (sol.size() != fNState){
        std::cout << "The solution vector size is different from the number of state variables. Please check it. \n";
        PanicButton();
    }
    for (int64_t inode = 0; inode < NNodes(); inode++){
        if (fNodeVector[inode]->HasBC()) continue;

        for (int istate = 0; istate < fNState; istate++){
            fConnectVector[inode]->SetSolution(istate,sol[istate]);
            fConnectVector[inode]->SetPreviousSolution(istate,sol[istate]);
        }
    }
}

void CompMesh::BuildMesh(){
    BuildConnects();
    
    //If there is any transient material, allocate the time derivatives for all connects
    for (int i = 0; i < fMaterialVector.size(); i++)
    {
        TransientWeakForm *transientmaterial = dynamic_cast<TransientWeakForm * > (fMaterialVector[i]);
        if (transientmaterial){
            for (int j = 0; j < fConnectVector.size(); j++){
                fConnectVector[j]->AllocateTimeDerivatives();
            }
            break;
        }
    }
    
}

void CompMesh::BuildConnects(){
    int nconnects = 0;

    if (fApproxType == ApproxType::EIsogeometric){
        std::cout << "The Isogeometric Analysis is not implemented yet. Please choose another approximation type. \n";
        PanicButton();
    } else if (fApproxType == ApproxType::EHierarquic){
        std::cout << "The Hierarquic basis functions are not implemented yet. Please choose another approximation type. \n";
        PanicButton();
    } else if (fApproxType == ApproxType::EIsoparametric){
        int nconnects = NNodes();
        fConnectVector.resize(nconnects);
        for (int64_t i = 0; i < NNodes(); i++){
            fConnectVector[i] = new Connect(fNState,1,fOrder,i);
        }
    } else {
        std::cout << "Unknown approximation type. Please check it. \n";
        PanicButton();
    }
    
}