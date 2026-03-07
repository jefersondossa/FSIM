#include "MixedCompMesh.h"
#include "ElementMixed.h"

MixedCompMesh::MixedCompMesh(std::vector<CompMesh* > meshvector){
    fMeshVector = meshvector;
    fNSpaces = meshvector.size();
}


void MixedCompMesh::AutoBuild(){
 
   CreateMixedConnects();
}

void MixedCompMesh::CreateMixedConnects(){
    
    int64_t nEl0 = fMeshVector[0]->NElements();
    int64_t nEl1 = fMeshVector[1]->NElements();    

#ifdef DEBUG_BUILD
    if (nEl0 != nEl1) {
        std::cout << "Error: The number of elements in the meshes must be the same." << std::endl;
        PanicButton();
    }
#endif
    
this->SetNumElements(nEl0);

    for (int64_t iel = 0; iel < nEl0; iel++){
        auto el0 = fMeshVector[0]->ElementVec()[iel];
        auto el1 = fMeshVector[1]->ElementVec()[iel];

#ifdef DEBUG_BUILD
        if (el0->Reference()->Index() != el1->Reference()->Index()) {
            std::cout << "Error: The reference element of the subelements must be the same." << std::endl;
            PanicButton();
        }
#endif  
        std::vector<Element *> elvector = {el0, el1};

        ElementMixed *fMixedEl = new ElementMixed(iel, elvector, this, this->Material(0)); 
    }

    // int64_t nconnects = 0;
    // for (int i = 0; i < fNSpaces; i++){
    //     nconnects += fMeshVector[i]->NConnects();
    // }
    // SetNumConnects(nconnects);

    // int64_t count = 0;
    // for (int i = 0; i < fNSpaces; i++){
    //     for (int j = 0; j < fMeshVector[i]->NConnects(); j++){
    //         Connect *c = fMeshVector[i]->ConnectVec()[j];
    //         c->SetSequenceNumber(count);
    //         ConnectVec()[count] = c;
    //         count++;
    //     }
    // }

}