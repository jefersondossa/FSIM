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
        PanicButton(); //We don't know how solve a problem with two materials with different state variables in the same mesh
    }
}