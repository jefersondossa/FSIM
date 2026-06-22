#include "GraphMesh.h"


// GraphMesh::GraphMesh(CompMesh *cmesh){
//     int nElements = cmesh->NElements();

//     fElementConnects.resize(nElements);
//     fElementTypes.resize(nElements);

//     int nodecount = 0;

//     for (int iel = 0; iel < nElements; iel++){
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
    int nElements =0;

    for (int iel = 0; iel < cmesh->NElements(); iel++){
        auto cel = cmesh->ElementVec()[iel];
        if (cel->Dimension() != cmesh->Dimension()) continue;
        if (cel->GetWeakForm() == nullptr) continue;
        nElements++;
    }
    

    fElementConnects.resize(nElements);
    fElementTypes.resize(nElements);

    int nodecount = 0;
    int elcount = 0;
    for (int iel = 0; iel < cmesh->NElements(); iel++){
        auto cel = cmesh->ElementVec()[iel];
        if (cel->Dimension() != cmesh->Dimension()) continue;
        if (cel->GetWeakForm() == nullptr) continue;
        fGElementToMElement[elcount] = iel;
        int nelnodes = cel->Reference()->NGeometricNodes();
        fElementTypes[elcount] = cel->Reference()->PrintType();
        fElementConnects[elcount].resize(nelnodes);
        for (int inode = 0; inode < nelnodes; inode++){
            fNodes.push_back(cmesh->Reference()->NodeVec()[cel->Reference()->getGeometricNodes()[inode]]->getCoordinates());
            fGNodeToMNode[nodecount] = cel->Reference()->getGeometricNodes()[inode];
            fMNodeToGNode[cel->Reference()->getGeometricNodes()[inode]] = nodecount;
            fElementConnects[elcount][inode] = nodecount;
            nodecount++;
        }
        elcount++;
    }
    
    

}



GraphMesh::GraphMesh(CompMesh *cmesh, std::set<int> &elsToPrint){
    int nElements = elsToPrint.size();

    fElementConnects.resize(nElements);
    fElementTypes.resize(nElements);

    int nodecount = 0;
    int elcount = 0;
    for (auto iel : elsToPrint){
        auto cel = cmesh->ElementVec()[iel];
        fGElementToMElement[elcount] = iel;
        int nelnodes = cel->Reference()->NGeometricNodes();
        fElementTypes[elcount] = cel->Reference()->PrintType();
        fElementConnects[elcount].resize(nelnodes);
        for (int inode = 0; inode < nelnodes; inode++){
            fNodes.push_back(cmesh->Reference()->NodeVec()[cel->Reference()->getGeometricNodes()[inode]]->getCoordinates());
            fGNodeToMNode[nodecount] = cel->Reference()->getGeometricNodes()[inode];
            fMNodeToGNode[cel->Reference()->getGeometricNodes()[inode]] = nodecount;
            fElementConnects[elcount][inode] = nodecount;
            nodecount++;
        }
        elcount++;
    }
    
    

}