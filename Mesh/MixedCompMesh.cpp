#include "MixedCompMesh.h"

MixedCompMesh::MixedCompMesh(std::vector<CompMesh* > meshvector){
    fMeshVector = meshvector;
    fNSpaces = meshvector.size();
}


void MixedCompMesh::AutoBuild(){
 
   CreateMixedConnects();
}

void MixedCompMesh::CreateMixedConnects(){
    
    int64_t nconnects = 0;
    for (int i = 0; i < fNSpaces; i++){
        nconnects += fMeshVector[i]->NConnects();
    }
    SetNumConnects(nconnects);

    int64_t count = 0;
    for (int i = 0; i < fNSpaces; i++){
        for (int j = 0; j < fMeshVector[i]->NConnects(); j++){
            Connect *c = fMeshVector[i]->ConnectVec()[j];
            c->SetSequenceNumber(count);
            ConnectVec()[count] = c;
            count++;
        }
    }

}