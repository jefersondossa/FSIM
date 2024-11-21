#include "CompMesh.h"

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
            fNodeVector[inode]->SetSolution(istate,sol[istate]);
            fNodeVector[inode]->SetPreviousSolution(istate,sol[istate]);
        }
    }
}
