#include "CompMesh.h"
#include "Connect.h"
#include "TransientWeakForm.h"
#include "PlasticityModel.h"
#include "HierarquicalOneD.h"
#include "ElementT.h"
#include "ElementWithMem.h"
#include "ElementTransient.h"

CompMesh::CompMesh(GeoMesh *gmesh, ApproxType approxType){
    fReference = gmesh;
    fApproxType = approxType;
}

void CompMesh::AutoBuild(){
    // Create the elements based on the geometric mesh
    BuildElements();

    //Build the computational mesh
    BuildConnectivity();

    part_elem= new int[NElements()]();
}


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
    for (int64_t iconnect = 0; iconnect < NConnects() ; iconnect++){
        for (int istate = 0; istate < fNState; istate++){
            fConnectVector[iconnect]->SetSolution(istate,sol[istate]);
            fConnectVector[iconnect]->SetPreviousSolution(istate,sol[istate]);
        }
    }
}

void CompMesh::BuildConnectivity(){

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
    int64_t seqnum = 0;

    switch (fApproxType)
    {
    case ApproxType::EHierarquic:
        BuildHierarquicConnects();
        break;
    case ApproxType::EIsoparametric:
        nconnects = fReference->NNodes();
        fConnectVector.resize(nconnects);
        for (int64_t i = 0; i < nconnects; i++){
            fConnectVector[i] = new Connect(fNState,1,fOrder,i,seqnum);
            seqnum += fNState;
        }
        for (int64_t iel=0; iel < NElements(); iel++){
            Element *el = fElementVector[iel];
            VecInt &geoNodes = el->Reference()->getGeometricNodes();
            el->setConnectivity(geoNodes);
        }
        fNGlobalDOF = nconnects * fNState;
        break;
    default:
        std::cout << "Unknown approximation type. Please check it. \n";
        PanicButton();
        break;
    }

    
}

void CompMesh::BuildHierarquicConnects(){

    // TODO: Decidir como os connects serão criados. Uma opção é criar os connects
    // de nós primeiro e depois os demais. Isso faz com que a banda da matriz aumente,
    // mas possibilita que não sejam realizadas mudanças na implementação da impressão dos
    // resultados. Outra opção é sempre criar uma graph mesh, que teria conectividade
    // própria e não necessariamente coincidente com a do comp mesh. Isso pode ser mais trabalhoso
    // para implementar, mas manteria a banda da matriz menor.

    int nconnects = 0;
    fConnectVector.reserve(NElements());
    std::map<int,int> node_to_connect;
    std::map<std::vector<int>,int> edge_to_connect;
    std::map<std::vector<int>,int> face_to_connect;
    std::map<std::vector<int>,int> volume_to_connect;

    int edgecount = 0;
    int facecount = 0;
    int volumecount = 0;
    int64_t seqnum = 0;

    for (auto iel = NElements()-1; iel >= 0; iel--){
        Element *el = fElementVector[iel];
        VecInt &geoNodes = el->Reference()->getGeometricNodes();
        int ncorner = el->Reference()->NCornerNodes();
        int nedges = el->Reference()->NEdges();
        int nfaces = el->Reference()->NFaces();
        int nvolumes = el->Reference()->NVolumes();
        int nsides = el->NSides();
        VecInt connect(nsides);
        connect.setZero(); 

        //Start the connectivity with the corner nodes
        for (int i = 0; i < ncorner; i++){
            if (node_to_connect.find(geoNodes[i]) == node_to_connect.end()){
                fConnectVector.push_back(new Connect(fNState,1,1,nconnects,seqnum));
                seqnum += fNState; // Nodes have only 1 associated shape function
                node_to_connect[geoNodes[i]] = nconnects;
                connect[i] = node_to_connect[geoNodes[i]];
                nconnects++;
            } else {
                connect[i] = node_to_connect[geoNodes[i]];
            }
        }

        //Create connects for the edges
        for (int i = ncorner; i < ncorner + nedges; i++){
            std::vector<int> sideNodes(2), sideNodesSorted(2);
            for (int j = 0; j < 2; j++){
                sideNodes[j] = geoNodes[el->SideNodeLocIndex(i,j)];
            }
            sideNodesSorted[0] = sideNodes[1];
            sideNodesSorted[1] = sideNodes[0];

            //Varrer todos os elementos vizinhos e verificar se os dois nos participam
            //do vetor geo nodes. Se sim, então há um connect comum para os dois.
            // Caso contrário, criar um novo connect de aresta.
            int numNeig=el->Reference()->getNumberOfNeighborElements(); 
            for (int ineig = 0; ineig < numNeig; ineig++){
                auto neig = el->Reference()->getNeighborElement(ineig);
                if (!fElementVector[neig]) continue;
                if (fElementVector[neig]->Dimension() == 0) continue;

                VecInt &neigNodes = fElementVector[neig]->Reference()->getGeometricNodes();
                std::set<int> neighNodesVec(neigNodes.data(),neigNodes.data()+neigNodes.size());
                if (neighNodesVec.find(sideNodes[0])!= neighNodesVec.end() &&
                    neighNodesVec.find(sideNodes[1])!= neighNodesVec.end()){
                    //achou um elemento vizinho que tem os dois nos da aresta
                    //verificar se ja existe um connect para essa aresta
                    if (edge_to_connect.find(sideNodes) == edge_to_connect.end() &&
                        edge_to_connect.find(sideNodesSorted) == edge_to_connect.end()){
                        int nshape = HierarquicalOneD::NShapeFunctions(2,fOrder);
                        int64_t ef_seqnum = nshape == 0 ? -1 : seqnum;
                        fConnectVector.push_back(new Connect(fNState,nshape,fOrder,nconnects,ef_seqnum));
                        seqnum += nshape * fNState;
                        edge_to_connect[sideNodes] = nconnects;
                        edge_to_connect[sideNodesSorted] = nconnects;
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
            int nsidenodes = el->NSideNodes(i);
            std::vector<int> sideNodes(nsidenodes);
            for (int j = 0; j < nsidenodes; j++){
                sideNodes[j] = geoNodes[j];
            }

            //Varrer todos os elementos vizinhos e verificar se os dois nos participam
            //do vetor geo nodes. Se sim, então há um connect comum para os dois.
            // Caso contrário, criar um novo connect de aresta.
            int numNeig=el->Reference()->getNumberOfNeighborElements(); 
            for (int ineig = 0; ineig < numNeig; ineig++){
                auto neig = el->Reference()->getNeighborElement(ineig);
                if (!fElementVector[neig]) continue;
                if (fElementVector[neig]->Dimension() == 0) continue;

                VecInt &neigNodes = fElementVector[neig]->Reference()->getGeometricNodes();
                std::set<int> neighNodesVec(neigNodes.data(),neigNodes.data()+neigNodes.size());

                bool allNodesFound = true;
                for (int j = 0; j < nsidenodes; j++){
                    if (neighNodesVec.find(sideNodes[j]) == neighNodesVec.end()){
                        allNodesFound = false;
                        break;
                    }
                }

                if (allNodesFound){
                    //achou um elemento vizinho que tem os dois nos da aresta
                    //verificar se ja existe um connect para essa aresta
                    if (face_to_connect.find(sideNodes) == face_to_connect.end()){


                        int nshape = HierarquicalOneD::NShapeFunctions(2,fOrder);
                        
                        
                        int64_t ef_seqnum = nshape == 0 ? -1 : seqnum;
                        fConnectVector.push_back(new Connect(fNState,nshape,fOrder,nconnects,ef_seqnum));
                        seqnum += nshape * fNState;
                        face_to_connect[sideNodes] = nconnects;
                        connect[i] = nconnects;
                        nconnects++;
                    } else {
                        //ja existe um connect para essa face
                        connect[i] = face_to_connect[sideNodes];
                    }
                }
            }
        }


        el->setConnectivity(connect);

    }
    
    fNGlobalDOF = 0;
    for (int64_t i = 0; i < fConnectVector.size(); i++){
        fNGlobalDOF += fConnectVector[i]->GetNShapeFunctions() * fNState;
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
    for (int64_t i = 0; i < NElements(); i++)
    {
        Element *el = fElementVector[i];
        file << "Element " << i << ": ";
        file << "Material ID = " << el->GetWeakForm()->Id() << ", ";
        file << "Geometric Nodes = [";
        VecInt &geoNodes = el->Reference()->getGeometricNodes();
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



#include "ShapeHexahedron.h"
#include "ShapeOneDLin.h"
#include "ShapeOneDQua.h"
#include "ShapeOneDCub.h"
#include "ShapeQuadrilateralLin.h"
#include "ShapeQuadrilateralQua.h"
#include "ShapePoint.h"
#include "ShapeTetrahedronLin.h"
#include "ShapeTetrahedronQua.h"
#include "ShapeTetrahedronCub.h"
#include "ShapeTriangleLin.h"
#include "ShapeTriangleQua.h"
#include "ShapeTriangleCub.h"
#include "HierarquicalOneD.h"
#include "HierarquicalQuad.h"
#include "HierarquicalTriangle.h"

void CompMesh::BuildElements(){
    fElementVector.resize(fReference->NElements());
    for (int64_t iel = 0; iel < NElements(); iel++){
        int elType = fReference->ElementVec()[iel]->PrintType();
        int matid = fReference->ElementVec()[iel]->Material();

        PlasticityModel *plasticmaterial = dynamic_cast<PlasticityModel * > (fMaterialVector[matid]);
        TransientWeakForm *transientmaterial = dynamic_cast<TransientWeakForm * > (fMaterialVector[matid]);

        switch (elType)
        {
        case 3: 
        {//Linear Line 
            if (plasticmaterial){
                switch (fApproxType){
                case ApproxType::EIsoparametric:
                    fElementVector[iel] = new ElementWithMem<ShapeOneDLin>(iel,fReference->ElementVec()[iel],this,Material(fReference->ElementVec()[iel]->Material()));
                    break;
                case ApproxType::EHierarquic:
                    fElementVector[iel] = new ElementWithMem<HierarquicalOneD>(iel,fReference->ElementVec()[iel],this,Material(fReference->ElementVec()[iel]->Material()));
                
                default:
                    break;
                }
            } else if (transientmaterial){
                switch (fApproxType){
                case ApproxType::EIsoparametric:
                    fElementVector[iel] = new ElementTransient<ShapeOneDLin>(iel,fReference->ElementVec()[iel],this,Material(fReference->ElementVec()[iel]->Material()));
                    break;
                case ApproxType::EHierarquic:
                    fElementVector[iel] = new ElementTransient<HierarquicalOneD>(iel,fReference->ElementVec()[iel],this,Material(fReference->ElementVec()[iel]->Material()));
                
                default:
                    break;
                }
            } else {
                switch (fApproxType){
                case ApproxType::EIsoparametric:
                    fElementVector[iel] = new ElementT<ShapeOneDLin>(iel,fReference->ElementVec()[iel],this,Material(fReference->ElementVec()[iel]->Material()));
                    break;
                case ApproxType::EHierarquic:
                    fElementVector[iel] = new ElementT<HierarquicalOneD>(iel,fReference->ElementVec()[iel],this,Material(fReference->ElementVec()[iel]->Material()));
                
                default:
                    break;
                }
            }
            break;
        }

        case 5: 
        {//Linear Triangle
            if (plasticmaterial){
                switch (fApproxType){
                case ApproxType::EIsoparametric:
                    fElementVector[iel] = new ElementWithMem<ShapeTriangleLin>(iel,fReference->ElementVec()[iel],this,Material(fReference->ElementVec()[iel]->Material()));
                    break;
                case ApproxType::EHierarquic:
                    fElementVector[iel] = new ElementWithMem<HierarquicalTriangle>(iel,fReference->ElementVec()[iel],this,Material(fReference->ElementVec()[iel]->Material()));
                
                default:
                    break;
                }
            } else if (transientmaterial){
                switch (fApproxType){
                case ApproxType::EIsoparametric:
                    fElementVector[iel] = new ElementTransient<ShapeTriangleLin>(iel,fReference->ElementVec()[iel],this,Material(fReference->ElementVec()[iel]->Material()));
                    break;
                case ApproxType::EHierarquic:
                    fElementVector[iel] = new ElementTransient<HierarquicalTriangle>(iel,fReference->ElementVec()[iel],this,Material(fReference->ElementVec()[iel]->Material()));
                
                default:
                    break;
                }
            } else {
                switch (fApproxType){
                case ApproxType::EIsoparametric:
                    fElementVector[iel] = new ElementT<ShapeTriangleLin>(iel,fReference->ElementVec()[iel],this,Material(fReference->ElementVec()[iel]->Material()));
                    break;
                case ApproxType::EHierarquic:
                    fElementVector[iel] = new ElementT<HierarquicalTriangle>(iel,fReference->ElementVec()[iel],this,Material(fReference->ElementVec()[iel]->Material()));
                
                default:
                    break;
                }
            }
            break;
        }

        case 9:
        {// Linear Quadrilateral
            if (plasticmaterial){
                switch (fApproxType){
                case ApproxType::EIsoparametric:
                    fElementVector[iel] = new ElementWithMem<ShapeQuadrilateralLin>(iel,fReference->ElementVec()[iel],this,Material(fReference->ElementVec()[iel]->Material()));
                    break;
                case ApproxType::EHierarquic:
                    fElementVector[iel] = new ElementWithMem<HierarquicalQuad>(iel,fReference->ElementVec()[iel],this,Material(fReference->ElementVec()[iel]->Material()));
                
                default:
                    break;
                }
            } else if (transientmaterial){
                switch (fApproxType){
                case ApproxType::EIsoparametric:
                    fElementVector[iel] = new ElementTransient<ShapeQuadrilateralLin>(iel,fReference->ElementVec()[iel],this,Material(fReference->ElementVec()[iel]->Material()));
                    break;
                case ApproxType::EHierarquic:
                    fElementVector[iel] = new ElementTransient<HierarquicalQuad>(iel,fReference->ElementVec()[iel],this,Material(fReference->ElementVec()[iel]->Material()));
                
                default:
                    break;
                }
            } else {
                switch (fApproxType){
                case ApproxType::EIsoparametric:
                    fElementVector[iel] = new ElementT<ShapeQuadrilateralLin>(iel,fReference->ElementVec()[iel],this,Material(fReference->ElementVec()[iel]->Material()));
                    break;
                case ApproxType::EHierarquic:
                    fElementVector[iel] = new ElementT<HierarquicalQuad>(iel,fReference->ElementVec()[iel],this,Material(fReference->ElementVec()[iel]->Material()));
                
                default:
                    break;
                }
            }
            break;
        }

        case 10:
        {//Linear tetrahedron
            if (plasticmaterial){
                switch (fApproxType){
                case ApproxType::EIsoparametric:
                    fElementVector[iel] = new ElementWithMem<ShapeTetrahedronLin>(iel,fReference->ElementVec()[iel],this,Material(fReference->ElementVec()[iel]->Material()));
                    break;
                case ApproxType::EHierarquic:
                    PanicButton(); //Hierarquic tetrahedron not implemented yet
                default:
                    break;
                }
            } else if (transientmaterial){
                switch (fApproxType){
                case ApproxType::EIsoparametric:
                    fElementVector[iel] = new ElementTransient<ShapeTetrahedronLin>(iel,fReference->ElementVec()[iel],this,Material(fReference->ElementVec()[iel]->Material()));
                    break;
                case ApproxType::EHierarquic:
                    PanicButton(); //Hierarquic tetrahedron not implemented yet
                
                default:
                    break;
                }
            } else {
                switch (fApproxType){
                case ApproxType::EIsoparametric:
                    fElementVector[iel] = new ElementT<ShapeTetrahedronLin>(iel,fReference->ElementVec()[iel],this,Material(fReference->ElementVec()[iel]->Material()));
                    break;
                case ApproxType::EHierarquic:
                    PanicButton(); //Hierarquic tetrahedron not implemented yet
                
                default:
                    break;
                }
            }
            break;
        }
            
        case 12:
        {// Linear Hexahedron
            if (plasticmaterial){
                switch (fApproxType){
                case ApproxType::EIsoparametric:
                    fElementVector[iel] = new ElementWithMem<ShapeHexahedron>(iel,fReference->ElementVec()[iel],this,Material(fReference->ElementVec()[iel]->Material()));
                    break;
                case ApproxType::EHierarquic:
                    PanicButton(); //Hierarquic hexahedron not implemented yet
                default:
                    break;
                }
            } else if (transientmaterial){
                switch (fApproxType){
                case ApproxType::EIsoparametric:
                    fElementVector[iel] = new ElementTransient<ShapeHexahedron>(iel,fReference->ElementVec()[iel],this,Material(fReference->ElementVec()[iel]->Material()));
                    break;
                case ApproxType::EHierarquic:
                    PanicButton(); //Hierarquic hexahedron not implemented yet
                
                default:
                    break;
                }
            } else {
                switch (fApproxType){
                case ApproxType::EIsoparametric:
                    fElementVector[iel] = new ElementT<ShapeHexahedron>(iel,fReference->ElementVec()[iel],this,Material(fReference->ElementVec()[iel]->Material()));
                    break;
                case ApproxType::EHierarquic:
                    PanicButton(); //Hierarquic hexahedron not implemented yet
                
                default:
                    break;
                }
            }
            break;
        }

        case 21:
        {// Quadratic Line
            if (plasticmaterial){
                switch (fApproxType){
                case ApproxType::EIsoparametric:
                    fElementVector[iel] = new ElementWithMem<ShapeOneDQua>(iel,fReference->ElementVec()[iel],this,Material(fReference->ElementVec()[iel]->Material()));
                    break;
                case ApproxType::EHierarquic:
                    fElementVector[iel] = new ElementWithMem<HierarquicalOneD>(iel,fReference->ElementVec()[iel],this,Material(fReference->ElementVec()[iel]->Material()));
                
                default:
                    break;
                }
            } else if (transientmaterial){
                switch (fApproxType){
                case ApproxType::EIsoparametric:
                    fElementVector[iel] = new ElementTransient<ShapeOneDQua>(iel,fReference->ElementVec()[iel],this,Material(fReference->ElementVec()[iel]->Material()));
                    break;
                case ApproxType::EHierarquic:
                    fElementVector[iel] = new ElementTransient<HierarquicalOneD>(iel,fReference->ElementVec()[iel],this,Material(fReference->ElementVec()[iel]->Material()));
                
                default:
                    break;
                }
            } else {
                switch (fApproxType){
                case ApproxType::EIsoparametric:
                    fElementVector[iel] = new ElementT<ShapeOneDQua>(iel,fReference->ElementVec()[iel],this,Material(fReference->ElementVec()[iel]->Material()));
                    break;
                case ApproxType::EHierarquic:
                    fElementVector[iel] = new ElementT<HierarquicalOneD>(iel,fReference->ElementVec()[iel],this,Material(fReference->ElementVec()[iel]->Material()));
                
                default:
                    break;
                }
            }
            break;
        }

        case 22:
        {// Quadratic Triangle
            if (plasticmaterial){
                switch (fApproxType){
                case ApproxType::EIsoparametric:
                    fElementVector[iel] = new ElementWithMem<ShapeTriangleQua>(iel,fReference->ElementVec()[iel],this,Material(fReference->ElementVec()[iel]->Material()));
                    break;
                case ApproxType::EHierarquic:
                    fElementVector[iel] = new ElementWithMem<HierarquicalTriangle>(iel,fReference->ElementVec()[iel],this,Material(fReference->ElementVec()[iel]->Material()));
                
                default:
                    break;
                }
            } else if (transientmaterial){
                switch (fApproxType){
                case ApproxType::EIsoparametric:
                    fElementVector[iel] = new ElementTransient<ShapeTriangleQua>(iel,fReference->ElementVec()[iel],this,Material(fReference->ElementVec()[iel]->Material()));
                    break;
                case ApproxType::EHierarquic:
                    fElementVector[iel] = new ElementTransient<HierarquicalTriangle>(iel,fReference->ElementVec()[iel],this,Material(fReference->ElementVec()[iel]->Material()));
                
                default:
                    break;
                }
            } else {
                switch (fApproxType){
                case ApproxType::EIsoparametric:
                    fElementVector[iel] = new ElementT<ShapeTriangleQua>(iel,fReference->ElementVec()[iel],this,Material(fReference->ElementVec()[iel]->Material()));
                    break;
                case ApproxType::EHierarquic:
                    fElementVector[iel] = new ElementT<HierarquicalTriangle>(iel,fReference->ElementVec()[iel],this,Material(fReference->ElementVec()[iel]->Material()));
                
                default:
                    break;
                }
            }
            break;
        }

        case 28:
        {// Quadratic Quad
            if (plasticmaterial){
                switch (fApproxType){
                case ApproxType::EIsoparametric:
                    fElementVector[iel] = new ElementWithMem<ShapeQuadrilateralQua>(iel,fReference->ElementVec()[iel],this,Material(fReference->ElementVec()[iel]->Material()));
                    break;
                case ApproxType::EHierarquic:
                    fElementVector[iel] = new ElementWithMem<HierarquicalQuad>(iel,fReference->ElementVec()[iel],this,Material(fReference->ElementVec()[iel]->Material()));
                
                default:
                    break;
                }
            } else if (transientmaterial){
                switch (fApproxType){
                case ApproxType::EIsoparametric:
                    fElementVector[iel] = new ElementTransient<ShapeQuadrilateralQua>(iel,fReference->ElementVec()[iel],this,Material(fReference->ElementVec()[iel]->Material()));
                    break;
                case ApproxType::EHierarquic:
                    fElementVector[iel] = new ElementTransient<HierarquicalQuad>(iel,fReference->ElementVec()[iel],this,Material(fReference->ElementVec()[iel]->Material()));
                
                default:
                    break;
                }
            } else {
                switch (fApproxType){
                case ApproxType::EIsoparametric:
                    fElementVector[iel] = new ElementT<ShapeQuadrilateralQua>(iel,fReference->ElementVec()[iel],this,Material(fReference->ElementVec()[iel]->Material()));
                    break;
                case ApproxType::EHierarquic:
                    fElementVector[iel] = new ElementT<HierarquicalQuad>(iel,fReference->ElementVec()[iel],this,Material(fReference->ElementVec()[iel]->Material()));
                
                default:
                    break;
                }
            }
            break;
        }

        case 1:
        {//Point
            if (plasticmaterial){
                fElementVector[iel] = new ElementWithMem<ShapePoint>(iel,fReference->ElementVec()[iel],this,Material(fReference->ElementVec()[iel]->Material()));
            } else if (transientmaterial){
                fElementVector[iel] = new ElementTransient<ShapePoint>(iel,fReference->ElementVec()[iel],this,Material(fReference->ElementVec()[iel]->Material()));
            } else {
                fElementVector[iel] = new ElementT<ShapePoint>(iel,fReference->ElementVec()[iel],this,Material(fReference->ElementVec()[iel]->Material()));
            }
            break;
        }

        case 69:
        {// Cubic Triangle
            if (plasticmaterial){
                switch (fApproxType){
                case ApproxType::EIsoparametric:
                    fElementVector[iel] = new ElementWithMem<ShapeTriangleCub>(iel,fReference->ElementVec()[iel],this,Material(fReference->ElementVec()[iel]->Material()));
                    break;
                case ApproxType::EHierarquic:
                    fElementVector[iel] = new ElementWithMem<HierarquicalTriangle>(iel,fReference->ElementVec()[iel],this,Material(fReference->ElementVec()[iel]->Material()));
                
                default:
                    break;
                }
            } else if (transientmaterial){
                switch (fApproxType){
                case ApproxType::EIsoparametric:
                    fElementVector[iel] = new ElementTransient<ShapeTriangleCub>(iel,fReference->ElementVec()[iel],this,Material(fReference->ElementVec()[iel]->Material()));
                    break;
                case ApproxType::EHierarquic:
                    fElementVector[iel] = new ElementTransient<HierarquicalTriangle>(iel,fReference->ElementVec()[iel],this,Material(fReference->ElementVec()[iel]->Material()));
                
                default:
                    break;
                }
            } else {
                switch (fApproxType){
                case ApproxType::EIsoparametric:
                    fElementVector[iel] = new ElementT<ShapeTriangleCub>(iel,fReference->ElementVec()[iel],this,Material(fReference->ElementVec()[iel]->Material()));
                    break;
                case ApproxType::EHierarquic:
                    fElementVector[iel] = new ElementT<HierarquicalTriangle>(iel,fReference->ElementVec()[iel],this,Material(fReference->ElementVec()[iel]->Material()));
                
                default:
                    break;
                }
            }
            break;
        }

        case 35:
        {// Cubic Line
            if (plasticmaterial){
                switch (fApproxType){
                case ApproxType::EIsoparametric:
                    fElementVector[iel] = new ElementWithMem<ShapeOneDCub>(iel,fReference->ElementVec()[iel],this,Material(fReference->ElementVec()[iel]->Material()));
                    break;
                case ApproxType::EHierarquic:
                    fElementVector[iel] = new ElementWithMem<HierarquicalOneD>(iel,fReference->ElementVec()[iel],this,Material(fReference->ElementVec()[iel]->Material()));
                
                default:
                    break;
                }
            } else if (transientmaterial){
                switch (fApproxType){
                case ApproxType::EIsoparametric:
                    fElementVector[iel] = new ElementTransient<ShapeOneDCub>(iel,fReference->ElementVec()[iel],this,Material(fReference->ElementVec()[iel]->Material()));
                    break;
                case ApproxType::EHierarquic:
                    fElementVector[iel] = new ElementTransient<HierarquicalOneD>(iel,fReference->ElementVec()[iel],this,Material(fReference->ElementVec()[iel]->Material()));
                
                default:
                    break;
                }
            } else {
                switch (fApproxType){
                case ApproxType::EIsoparametric:
                    fElementVector[iel] = new ElementT<ShapeOneDCub>(iel,fReference->ElementVec()[iel],this,Material(fReference->ElementVec()[iel]->Material()));
                    break;
                case ApproxType::EHierarquic:
                    fElementVector[iel] = new ElementT<HierarquicalOneD>(iel,fReference->ElementVec()[iel],this,Material(fReference->ElementVec()[iel]->Material()));
                
                default:
                    break;
                }
            }
            break;
        }

        default:
            std::cout << "Element type " << elType << " not implemented yet. Please check it. \n";
            PanicButton();
            break;
        }


    }
}
