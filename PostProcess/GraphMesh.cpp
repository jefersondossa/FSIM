#include "GraphMesh.h"


// GraphMesh::GraphMesh(CompMesh *cmesh){
//     int64_t nElements = cmesh->NElements();

//     fElementConnects.resize(nElements);
//     fElementTypes.resize(nElements);

//     int64_t nodecount = 0;

//     for (int64_t iel = 0; iel < nElements; iel++){
//         fElementTypes[iel] = cmesh->ElementVec()[iel]->PrintType();
//         fElementConnects[iel] = cmesh->ElementVec()[iel]->getConnectivity();
//         for (int iconnect = 0; iconnect < fElementConnects[iel].size(); iconnect++){
//             fNodes.push_back(cmesh->NodeVec()[fElementConnects[iel][iconnect]]->getCoordinates());
//             fGNodeToMNode[nodecount] = fElementConnects[iel][iconnect];
//             fMNodeToGNode[fElementConnects[iel][iconnect]] = nodecount;
//             fElementConnects[iel][iconnect] = nodecount;
//             nodecount++;
//         }
//     }
    
    

// }



GraphMesh::GraphMesh(CompMesh *cmesh){
    int64_t nElements = cmesh->NElements();

    fElementConnects.resize(nElements);
    fElementTypes.resize(nElements);

    int64_t nodecount = 0;

    for (int64_t iel = 0; iel < nElements; iel++){
        auto cel = cmesh->ElementVec()[iel];
        int nelnodes = cel->NElNodes();
        fElementTypes[iel] = cel->PrintType();
        fElementConnects[iel].resize(nelnodes);
        for (int inode = 0; inode < nelnodes; inode++){
            fNodes.push_back(cmesh->NodeVec()[cel->getGeometricNodes()[inode]]->getCoordinates());
            fGNodeToMNode[nodecount] = cel->getGeometricNodes()[inode];
            fMNodeToGNode[cel->getGeometricNodes()[inode]] = nodecount;
            fElementConnects[iel][inode] = nodecount;
            nodecount++;
        }
    }
    
    

}