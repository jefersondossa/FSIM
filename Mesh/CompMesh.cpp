#include "CompMesh.h"
#include "Connect.h"
#include "TransientWeakForm.h"
#include "HierarquicalOneD.h"

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

    switch (fApproxType)
    {
    case ApproxType::EHierarquic:
        BuildHierarquicConnects();
        break;
    case ApproxType::EIsoparametric:
        nconnects = NNodes();
        fConnectVector.resize(nconnects);
        for (int64_t i = 0; i < NNodes(); i++){
            fConnectVector[i] = new Connect(fNState,1,fOrder,i);
        }
        for (int64_t iel=0; iel < NElements(); iel++){
            Element *el = fElementVector[iel];
            VecInt &geoNodes = el->getGeometricNodes();
            el->setConnectivity(geoNodes);
        }
        break;
    default:
        std::cout << "Unknown approximation type. Please check it. \n";
        PanicButton();
        break;
    }

    
}

void CompMesh::BuildHierarquicConnects(){
    int nconnects = 0;
    fConnectVector.reserve(NNodes()+NElements());
    std::map<int,int> node_to_connect;
    std::map<std::vector<int>,int> edge_to_connect;
    std::map<std::vector<int>,int> face_to_connect;
    std::map<std::vector<int>,int> volume_to_connect;

    int edgecount = 0;
    int facecount = 0;
    int volumecount = 0;
    for (auto iel = 0; iel < NElements(); iel++){
        Element *el = fElementVector[iel];
        VecInt &geoNodes = el->getGeometricNodes();
        int ncorner = el->NCornerNodes();
        int nedges = el->NEdges();
        int nfaces = el->NFaces();
        int nvolumes = el->NVolumes();
        int nsides = el->NSides();
        VecInt connect(nsides);
        connect.setZero();

        //Start the connectivity with the corner nodes
        for (int i = 0; i < ncorner; i++){
            if (node_to_connect.find(geoNodes[i]) == node_to_connect.end()){
                fConnectVector.push_back(new Connect(fNState,1,1,nconnects));
                node_to_connect[geoNodes[i]] = nconnects;
                connect[i] = node_to_connect[geoNodes[i]];
                nconnects++;
            } else {
                connect[i] = node_to_connect[geoNodes[i]];
            }
            
        }

        //Create connects for the edges
        for (int i = ncorner; i < ncorner + nedges; i++){
            std::vector<int> sideNodes(2);
            sideNodes[0] = geoNodes[0];
            sideNodes[1] = geoNodes[1];

            //Varrer todos os elementos vizinhos e verificar se os dois nos participam
            //do vetor geo nodes. Se sim, então há um connect comum para os dois.
            // Caso contrário, criar um novo connect de aresta.
            int numNeig=el->getNumberOfNeighborElements(); 
            for (int ineig = 0; ineig < numNeig; ineig++){
                auto neig = el->getNeighborElement(ineig);
                if (!fElementVector[neig]) continue;
                if (fElementVector[neig]->Dimension() == 0) continue;

                VecInt &neigNodes = fElementVector[neig]->getGeometricNodes();
                std::set<int> neighNodesVec(neigNodes.data(),neigNodes.data()+neigNodes.size());
                if (neighNodesVec.find(sideNodes[0])!= neighNodesVec.end() &&
                    neighNodesVec.find(sideNodes[1])!= neighNodesVec.end()){
                    //achou um elemento vizinho que tem os dois nos da aresta
                    //verificar se ja existe um connect para essa aresta
                    if (edge_to_connect.find(sideNodes) == edge_to_connect.end()){
                        int nshape = HierarquicalOneD::NShapeFunctions(2,fOrder);
                        fConnectVector.push_back(new Connect(fNState,nshape,fOrder,nconnects));
                        edge_to_connect[sideNodes] = nconnects;
                        connect[i] = nconnects;
                        nconnects++;
                    } else {
                        //ja existe um connect para essa aresta
                        connect[i] = edge_to_connect[sideNodes];
                    }
                }
            }
        }

        //Create connects for the faces
        for (int i=ncorner+nedges; i < ncorner+nedges+nfaces; i++){
            PanicButton();
        }
        el->setConnectivity(connect);

    }
    


}

void CompMesh::Print(std::string filename){
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    // Create a map for ApproxType to string conversion
    static const std::map<ApproxType, std::string> approxTypeNames = {
        {ApproxType::EIsoparametric, "Isoparametric"},
        {ApproxType::EHierarquic, "Hierarquic"},
        {ApproxType::EIsogeometric, "Isogeometric"}
    };

    file << "Computational Mesh Information\n";
    file << "===============================\n";
    file << "Number of Nodes: " << NNodes() << "\n";
    file << "Number of Elements: " << NElements() << "\n";
    file << "Number of Connects: " << NConnects() << "\n";
    file << "Approximation Type: " << approxTypeNames.at(fApproxType) << "\n";
    file << "Polynomial Order: " << fOrder << "\n";
    file << "--------------------------------\n";


    file << "Connects Information:\n";
    for (int64_t i = 0; i < NConnects(); i++)
    {
        Connect *c = fConnectVector[i];
        file << "Connect " << i << ": ";
        file << "NState = " << c->GetNStateVariables() << ", ";
        file << "NShapeFunctions = " << c->GetNShapeFunctions() << ", ";
        file << "Order = " << c->GetOrder() << ", ";
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
    for (int64_t i = 0; i < NElements(); i++)
    {
        Element *el = fElementVector[i];
        file << "Element " << i << ": ";
        file << "Material ID = " << el->GetWeakForm()->Id() << ", ";
        file << "Geometric Nodes = [";
        VecInt &geoNodes = el->getGeometricNodes();
        for (size_t j = 0; j < geoNodes.size(); j++) {
            file << geoNodes[j];
            if (j < geoNodes.size() - 1) file << ", ";
        }
        file << "], Connects = [";
        VecInt &connects = el->getConnectivity();
        for (size_t j = 0; j < connects.size(); j++) {
            file << connects[j];
            if (j < connects.size() - 1) file << ", ";
        }
        file << "]\n";
    }  

    

    file.close();
}
