#include "MixedCompMesh.h"
#include "ElementMixed.h"

MixedCompMesh::MixedCompMesh(std::vector<CompMesh* > meshvector){
    fMeshVector = meshvector;
    fNSpaces = meshvector.size();
    this->fReference = fMeshVector[0]->Reference();

    this->fNGlobalDOF = 0;
    for (int i = 0; i < fNSpaces; i++){
        fNGlobalDOF += fMeshVector[i]->NGlobalDOF();
    }
    fApproxType = fMeshVector[0]->GetApproxType();
}


void MixedCompMesh::AutoBuild(){
 
   CreateMixedConnects();

   part_elem= new int[NElements()]();
}

void MixedCompMesh::CreateMixedConnects(){

    int nElements = fMeshVector[0]->NElements(); 
    int nEl0 = fMeshVector[0]->NElements();
    int nEl1 = fMeshVector[1]->NElements();

#ifdef DEBUG_BUILD
    if (nEl0 != nEl1) {
        std::cout << "nEl0: " << nEl0 << " nEl1: " << nEl1 << std::endl;
        std::cout << "Error: The number of elements in the meshes must be the same." << std::endl;
        PanicButton();
    }
#endif
    
    this->SetNumElements(nEl0);

    //maps the space to the connect
    std::map<int,int> fixedConnects;

    for (int iel = 0; iel < nEl0; iel++){
        std::vector<Element*> elvector(fNSpaces);

        for (int ispaces = 0; ispaces < fNSpaces; ispaces++){
            elvector[ispaces] = fMeshVector[ispaces]->ElementVec()[iel];
        }

#ifdef DEBUG_BUILD
        if (!elvector[0]) continue;
        if (elvector[0]->Reference()->Index() != elvector[0]->Reference()->Index()) {
            std::cout << "Error: The reference element of the subelements must be the same." << std::endl;
            PanicButton();
        }
#endif  
        int matid = elvector[0]->Reference()->Material();
        ElementMixed *fMixedEl = new ElementMixed(iel, elvector, this, this->Material(matid));
        int connsize = 0;
        for (int ispaces = 0; ispaces < fNSpaces; ispaces++){
            connsize += elvector[ispaces]->getConnectivity().size();
        }
        std::vector<Connect *> mixedconnect(connsize);
        int count = 0;
        for (int ispaces = 0; ispaces < fNSpaces; ispaces++){
            for (int i = 0; i < elvector[ispaces]->getConnectivity().size(); i++){
                mixedconnect[count++] = elvector[ispaces]->getConnectivity()[i];
            }
        }
        fMixedEl->setConnectivity(mixedconnect);

        this->ElementVec()[iel] = fMixedEl;
    }

    //fix the sequence number
    int seqnum = 0;
    for (int ispace = 1; ispace < fNSpaces; ispace++){
        seqnum += fMeshVector[ispace-1]->NGlobalDOF();
        for (int iconnect = 0; iconnect < fMeshVector[ispace]->NConnects(); iconnect++){
            int current_seqnum = fMeshVector[ispace]->ConnectVec()[iconnect]->GetSequenceNumber();
            if (current_seqnum < 0) continue;
            fMeshVector[ispace]->ConnectVec()[iconnect]->SetSequenceNumber(current_seqnum+seqnum);
        }
    }


}

void MixedCompMesh::Print(std::string filename){
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    file << "Computational Mesh Information\n";
    file << "===============================\n";
    file << "Number of Elements: " << NElements() << "\n";
    file << "Number of Spaces: " << fNSpaces << "\n\n\n";
    for (int ispace = 0; ispace < fNSpaces; ispace++){
        file << "===============================\n";
        file << "INFORMATION OF SPACE " << ispace << "\n";
        file << "Polynomial Order: " << fMeshVector[ispace]->GetDefaultOrder() << "\n";
        file << "--------------------------------\n";

        file << "Connects Information:\n";
        for (int i = 0; i < fMeshVector[ispace]->NConnects(); i++)
        {
            Connect *c = fMeshVector[ispace]->ConnectVec()[i];
            file << "Connect " << i << ": ";
            file << "NState = " << c->GetNStateVariables() << ", ";
            file << "NShapeFunctions = " << c->GetNShapeFunctions() << ", ";
            file << "Order = " << c->GetOrder() << ", ";
            file << "SeqNum = " << c->GetSequenceNumber() << ", ";
            file << "Solution = [";
            VecDouble &sol = c->Solution();
            for (size_t j = 0; j < sol.size(); j++) {
                file << sol[j];
                if (j < sol.size() - 1) file << ", ";
            }
            file << "]\n";
        }

        file << "--------------------------------\n";
        file << "Elements Information:\n";
        for (int i = 0; i < fMeshVector[ispace]->NElements(); i++)
        {
            Element *el = fMeshVector[ispace]->ElementVec()[i];
            if (!el) continue;
            file << "Element " << el->Index() << ": ";
            if (el->GetWeakForm()){
                file << "Material ID = " << el->GetWeakForm()->Id() << ", ";
            } else {
                file << "Material ID = None, ";
            }
            file << "Geometric Nodes = [";
            VecInt &geoNodes = el->Reference()->getGeometricNodes();
            for (size_t j = 0; j < geoNodes.size(); j++) {
                file << geoNodes[j];
                if (j < geoNodes.size() - 1) file << ", ";
            }
            file << "], Connects = [";
            VecInt connects = el->getConnectivityIndices();
            for (size_t j = 0; j < connects.size(); j++) {
                file << connects[j];
                if (j < connects.size() - 1) file << ", ";
            }
            file << "]\n";
        }  
    }
    

    file.close();
}
