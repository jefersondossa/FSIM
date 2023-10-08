#include "GmshTools.h"
#include "fluidDomain.h"
#include "ElPoisson.h"
#include "ElElasticity2D.h"
#include "ElStokes.h"
#include "ElNavierStokes.h"

#include<cstdlib>
#include<fstream>
#include<iostream>
#include <vector>

void GmshTools::MeshReading(Geometry* &geometry_, const std::string& mshfile,CompMesh *cmesh){

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //+++++++++++++++++++++++++++++OPPENING FILES+++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    //opening the .msh file
    std::ifstream file(mshfile);
    std::string line;
    std::getline(file, line); std::getline(file, line); std::getline(file, line); std::getline(file, line);
  
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++READIN MESH+++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //+++++++++++++++++++++++++++PHYSICAL ENTITIES++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    int number_physical_entities;
    file >> number_physical_entities;
    std::getline(file, line);
    std::unordered_map<int, std::string> physicalEntities;
    physicalEntities.reserve(number_physical_entities);

    for (int i = 0; i < number_physical_entities; i++)
    {
        std::getline(file, line);
        std::vector<std::string> tokens = SplitLine(line, " ");
        int index;
        std::istringstream(tokens[1]) >> index;
        physicalEntities[index] = tokens[2].substr(1, tokens[2].size() - 2);
    }
    std::getline(file, line); std::getline(file, line);

    ReadNodes(file,cmesh);
    ReadElements(geometry_,file,physicalEntities,cmesh);

    
    RenumberConnectivity(cmesh);

    BoundaryConstrains(cmesh);

    BoundarySides(cmesh);
    // // if (rank == 0) std::cout << "Number of elements " << number_elements << " " 
    //                          // << numElem << " " << numBoundElems << std::endl;
    // std::cout << std::endl << "Element Connectivity" << std::endl;        
    
    // for (int jel = 0; jel < cmesh->NElements(); jel++){
    //     VecInt connec = cmesh->ElementVec()[jel] -> getConnectivity();       
    //     for (int i=0; i < cmesh->NElNodes(); i++){
    //         std::cout << connec[i] << " ";
    //     };
    //     std::cout << std::endl;
    // };

}


std::vector<std::string> GmshTools::SplitLine(std::string str, std::string delim)
{
	std::istringstream is(str);
	std::vector<std::string> values;
	std::string token;
	while (getline(is, token, ' '))
		values.push_back(token);
	return values;
}


void GmshTools::ReadNodes(std::ifstream &file, CompMesh *cmesh){
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //+++++++++++++++++++++++++++++++++NODES++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    int64_t numNodes;
    std::string line;
    file >> numNodes;
    cmesh->NodeVec().resize(numNodes);
    std::getline(file, line);
    int64_t index = 0;

    std::cout << "Number of Nodes " << " " << numNodes << std::endl;
    for (int i = 0; i < numNodes; i++){
        VecDouble x(3);
        std::getline(file, line);
        std::vector<std::string> tokens = SplitLine(line, " ");
        
        for (int j = 0; j < 3; j++) std::istringstream(tokens[j+1]) >> x[j];

        int nstate = cmesh->NState();        

        Node *node = new Node(x,index,nstate);
        cmesh->NodeVec()[index] = node;
        index++;
    }
    std::getline(file, line); std::getline(file, line);
    

    int DIM = cmesh->Dimension();
    if (cmesh->getProblemParameters().ProbType() == ProblemType::ENavierStokes || cmesh->getProblemParameters().ProbType() == ProblemType::EStokes){
        cmesh->NGlobalDOF() = (DIM+1) * cmesh->NNodes();
    } else if (cmesh->getProblemParameters().ProbType() == ProblemType::EPoisson) {
        cmesh->NGlobalDOF() = cmesh->NNodes();
    } else if (cmesh->getProblemParameters().ProbType() == ProblemType::EElastic){
        cmesh->NGlobalDOF() = cmesh->NNodes() * DIM;
    } else {
        PanicButton();
    }
    

    return;
}


void GmshTools::ReadElements(Geometry* &geometry_, std::ifstream &file, std::unordered_map<int, std::string> &physicalEntities, CompMesh *cmesh){
    //defyning the maps that are used to store the elements information
    std::unordered_map<int, std::string> gmshElement = { {1, "line"}, {2, "triangle"}, {3, "quadrilateral"}, {8, "line3"}, {9, "triangle6"}, {10, "quadrilateral9"}, {15, "vertex"}, {16, "quadrilateral8"}, {20, "triangle9"}, {21, "triangle10"}, {26, "line4"}, {36, "quadrilateral16"}, {39, "quadrilateral12"} };
    std::unordered_map<std::string, int> numNodes2 = { {"vertex", 1}, {"line", 2}, {"triangle", 3}, {"quadrilateral", 4}, {"line3", 3}, {"triangle6", 6}, {"quadrilateral8", 8}, {"quadrilateral9", 9}, {"line4", 4}, {"triangle", 9}, {"triangle10", 10}, {"quadrilateral12", 12}, {"quadrilateral16", 16}};
    std::unordered_map<std::string, std::string> supportedElements = { {"triangle", "T3"}, {"triangle6", "T6"}, {"triangle10", "T10"}, {"quadrilateral", "Q4"}, {"quadrilateral8", "Q8"}, {"quadrilateral9", "Q9"}, {"quadrilateral12", "Q12"}, {"quadrilateral16", "Q16"}, {"tetrahedron4", "TET4"}, {"tetrahedron10", "TET10"}, {"tetrahedron20", "TET20"} };
    std::unordered_map<Line*, std::vector< std::vector<int> >> lineElements;
    std::string line;
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++ELEMENTS++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    int number_elements;
    file >> number_elements;
    cmesh->ElementVec().reserve(number_elements);
    cmesh->BoundaryVec().reserve(number_elements/10);
    int index = 0;
    std::getline(file, line);
    int cont = 0;

    int numBoundElems = 0;
    // numElem = 0;
    int numFSIInterfaces = 0;

    std::vector<BoundaryCondition*> dirichlet, neumann, glue, FSinterface;
    dirichlet = geometry_->getBoundaryCondition("DIRICHLET"); 
    neumann = geometry_->getBoundaryCondition("NEUMANN"); 
    glue = geometry_->getBoundaryCondition("GLUE");
    FSinterface = geometry_->getBoundaryCondition("FSINTERFACE");

    numFSIInterfaces = FSinterface.size();

    int nElNodes = cmesh->NElNodes();
    int nBdNodes = cmesh->NBdNodes();
    int DEG = cmesh->GetDefaultOrder();
    int DIM = cmesh->Dimension();
    for (int i = 0; i < number_elements; i++)
    {
        std::getline(file, line);
        std::vector<std::string> tokens = SplitLine(line, " ");
        std::vector<int> values(tokens.size(), 0);
        for (size_t j = 0; j < tokens.size(); j++)
            std::istringstream(tokens[j]) >> values[j];
        std::string elementType = gmshElement[values[1]];
        int number_nodes_per_element = numNodes2[elementType];
        std::vector<int> elementNodes;
        elementNodes.reserve(number_nodes_per_element);

        for (size_t j = 5 ; j < values.size(); j++)
            elementNodes.push_back(values[j]-1);
 
        std::string name = physicalEntities[values[3]];
        //Adding domain elements
        if (name[0] == 'v'){
            // if(rank == 0){
                Volume* object = geometry_ -> getVolume(name);
                VecInt connect(nElNodes);

                if (DEG == 2){
                    connect[2] = elementNodes[2];
                    connect[1] = elementNodes[1];
                    connect[3] = elementNodes[3];
                    connect[0] = elementNodes[0];
                    connect[5] = elementNodes[5];
                    connect[9] = elementNodes[8];
                    connect[6] = elementNodes[6];
                    connect[7] = elementNodes[7];
                    connect[4] = elementNodes[4];
                    connect[8] = elementNodes[9];
                } else {
                    for (int k = 0; k < nElNodes; k++) connect[k] = elementNodes[k];
                }
                switch (cmesh->getProblemParameters().ProbType())
                {
                case EPoisson:
                    {
                        ElPoisson *el = new ElPoisson(index++,connect,cmesh);
                        cmesh->ElementVec().push_back(el);
                    }
                    break;
                case EElastic:
                    {
                        ElElasticity2D *el = new ElElasticity2D(index++,connect,cmesh);
                        cmesh->ElementVec().push_back(el);
                    }
                    break;
                case EStokes:
                    {
                        ElStokes *el = new ElStokes(index++,connect,cmesh);
                        cmesh->ElementVec().push_back(el);
                    }
                    break;
                case ENavierStokes:
                    {
                        ElNavierStokes *el = new ElNavierStokes(index++,connect,cmesh);
                        cmesh->ElementVec().push_back(el);
                    }
                    break;
                
                default:
                    PanicButton();
                    break;
                }

                for (int k = 0; k < nElNodes; k++){
                    cmesh->NodeVec()[connect[k]] -> pushInverseIncidence(index);
                };
            // }
        }
        else if (name[0] == 's') {
            if (DIM == 3){
                VecInt connectB(nBdNodes);

                for (int i = 0; i < nBdNodes; i++) connectB[i] = elementNodes[i];

                int ibound;

                std::string::size_type sz;   // alias of size_t
                ibound = std::stoi (&name[1],nullptr,10);

                VecInt constrain(3);
                VecDouble value(3);

                for (int i = 0; i < dirichlet.size(); i++){
                    if (name == dirichlet[i] -> getLineName()){
                        if ((dirichlet[i] -> getComponentX()).size() == 0){
                            constrain[0] = 0; value[0] = 0;
                        }else{
                            std::vector<double> c = dirichlet[i] -> getComponentX();
                            constrain[0] = 1;
                            value[0] = c[0];
                        }
                        if ((dirichlet[i] -> getComponentY()).size() == 0){
                            constrain[1] = 0; value[1] = 0;
                        }else{
                            std::vector<double> c = dirichlet[i] -> getComponentY();
                            constrain[1] = 1;
                            value[1] = c[0];
                        }
                        if ((dirichlet[i] -> getComponentZ()).size() == 0){
                            constrain[2] = 0; value[2] = 0;
                        }else{
                            std::vector<double> c = dirichlet[i] -> getComponentZ();
                            constrain[2] = 1;
                            value[2] = c[0];
                        }
                    }
                }
                for (int i = 0; i < neumann.size(); i++){
                    if (name == neumann[i] -> getLineName()){
                        if ((neumann[i] -> getComponentX()).size() == 0){
                            constrain[0] = 0; value[0] = 0;
                        }else{
                            std::vector<double> c = neumann[i] -> getComponentX();
                            constrain[0] = 0;
                            value[0] = c[0];
                        }
                        if ((neumann[i] -> getComponentY()).size() == 0){
                            constrain[1] = 0; value[1] = 0;
                        }else{
                            std::vector<double> c = neumann[i] -> getComponentY();
                            constrain[1] = 0;
                            value[1] = c[0];
                        }
                        if ((neumann[i] -> getComponentZ()).size() == 0){
                            constrain[2] = 0; value[2] = 0;
                        }else{
                            std::vector<double> c = neumann[i] -> getComponentZ();
                            constrain[2] = 0;
                            value[2] = c[0];
                        }
                    }
                }  
                for (int i = 0; i < glue.size(); i++){
                    if (name == glue[i] -> getLineName()){
                        if ((glue[i] -> getComponentX()).size() == 0){
                            constrain[0] = 2; value[0] = 0;
                        }else{
                            std::vector<double> c = glue[i] -> getComponentX();
                            constrain[0] = 2;
                            value[0] = c[0];
                        }
                        if ((glue[i] -> getComponentY()).size() == 0){
                            constrain[1] = 2; value[1] = 0;
                        }else{
                            std::vector<double> c = glue[i] -> getComponentY();
                            constrain[1] = 2;
                            value[1] = c[0];
                        }
                        if ((glue[i] -> getComponentZ()).size() == 0){
                            constrain[2] = 2; value[2] = 0;
                        }else{
                            std::vector<double> c = glue[i] -> getComponentZ();
                            constrain[2] = 2;
                            value[2] = c[0];
                        }
                    }
                }              
                for (int i = 0; i < FSinterface.size(); i++){
                    if (name == FSinterface[i] -> getLineName()){
                        if ((FSinterface[i] -> getComponentX()).size() == 0){
                            constrain[0] = 3; value[0] = 0;
                        }else{
                            std::vector<double> c = FSinterface[i] -> getComponentX();
                            constrain[0] = 3;
                            value[0] = c[0];
                        }
                        if ((FSinterface[i] -> getComponentY()).size() == 0){
                            constrain[1] = 3; value[1] = 0;
                        }else{
                            std::vector<double> c = FSinterface[i] -> getComponentY();
                            constrain[1] = 3;
                            value[1] = c[0];
                        }
                        if ((FSinterface[i] -> getComponentZ()).size() == 0){
                            constrain[2] = 3; value[2] = 0;
                        }else{
                            std::vector<double> c = FSinterface[i] -> getComponentZ();
                            constrain[2] = 3;
                            value[2] = c[0];
                        }
                    }
                }  
                Boundary * bound = new Boundary(connectB, numBoundElems++, constrain, value, ibound, cmesh);
                // std::cout << "asdasd " << rank << " " << ibound << std::endl;
                cmesh->BoundaryVec().push_back(bound);
            } else {
                // if(rank == 0){
                    Surface* object = geometry_ -> getSurface(name);
                    VecInt connect(nElNodes);
                    for (int j = 0 ; j < nElNodes; j++) connect[j] = elementNodes[j];

                    switch (cmesh->getProblemParameters().ProbType())
                    {
                    case EPoisson:
                        {
                            ElPoisson *el = new ElPoisson(index++,connect,cmesh);
                            cmesh->ElementVec().push_back(el);
                        }
                        break;
                    case EElastic:
                        {
                            ElElasticity2D *el = new ElElasticity2D(index++,connect,cmesh);
                            cmesh->ElementVec().push_back(el);
                        }
                        break;
                    case EStokes:
                        {
                            ElStokes *el = new ElStokes(index++,connect,cmesh);
                            cmesh->ElementVec().push_back(el);
                        }
                        break;
                    case ENavierStokes:
                        {
                            ElNavierStokes *el = new ElNavierStokes(index++,connect,cmesh);
                            cmesh->ElementVec().push_back(el);
                        }
                        break;
                    
                    default:
                        PanicButton();
                        break;
                    }

                    for (int k = 0; k < nElNodes; k++){
                        cmesh->NodeVec()[connect[k]] -> pushInverseIncidence(index);
                    };
                // }
            }
        } else if ((name[0] == 'l') && (DIM == 2)) {
            VecInt connectB(nBdNodes);

            for (int i = 0; i < nBdNodes; i++) connectB[i] = elementNodes[i];
            
            int ibound;

            std::string::size_type sz;   // alias of size_t
            ibound = std::stoi (&name[1],nullptr,10);

            VecInt constrain(3);
            VecDouble value(3);

            for (int i = 0; i < dirichlet.size(); i++){
                if (name == dirichlet[i] -> getLineName()){
                    if ((dirichlet[i] -> getComponentX()).size() == 0){
                        constrain[0] = 0; value[0] = 0;
                    }else{
                        std::vector<double> c = dirichlet[i] -> getComponentX();
                        constrain[0] = 1;
                        value[0] = c[0];
                    }
                    if ((dirichlet[i] -> getComponentY()).size() == 0){
                        constrain[1] = 0; value[1] = 0;
                    }else{
                        std::vector<double> c = dirichlet[i] -> getComponentY();
                        constrain[1] = 1;
                        value[1] = c[0];
                    }
                }
            }
            for (int i = 0; i < neumann.size(); i++){
                if (name == neumann[i] -> getLineName()){
                    if ((neumann[i] -> getComponentX()).size() == 0){
                        constrain[0] = 0; value[0] = 0;
                    }else{
                        std::vector<double> c = neumann[i] -> getComponentX();
                        constrain[0] = 0;
                        value[0] = c[0];
                    }
                    if ((neumann[i] -> getComponentY()).size() == 0){
                        constrain[1] = 0; value[1] = 0;
                    }else{
                        std::vector<double> c = neumann[i] -> getComponentY();
                        constrain[1] = 0;
                        value[1] = c[0];
                    }
                }
            }  
            for (int i = 0; i < glue.size(); i++){
                if (name == glue[i] -> getLineName()){
                    if ((glue[i] -> getComponentX()).size() == 0){
                        constrain[0] = 2; value[0] = 0;
                    }else{
                        std::vector<double> c = glue[i] -> getComponentX();
                        constrain[0] = 2;
                        value[0] = c[0];
                    }
                    if ((glue[i] -> getComponentY()).size() == 0){
                        constrain[1] = 2; value[1] = 0;
                    }else{
                        std::vector<double> c = glue[i] -> getComponentY();
                        constrain[1] = 2;
                        value[1] = c[0];
                    }//std::cout <<"aqui " << std::endl;
                }
            }              
            for (int i = 0; i < FSinterface.size(); i++){
                if (name == FSinterface[i] -> getLineName()){
                    if ((FSinterface[i] -> getComponentX()).size() == 0){
                        constrain[0] = 3; value[0] = 0;
                    }else{
                        std::vector<double> c = FSinterface[i] -> getComponentX();
                        constrain[0] = 3;
                        value[0] = c[0];
                    }
                    if ((FSinterface[i] -> getComponentY()).size() == 0){
                        constrain[1] = 3; value[1] = 0;
                    }else{
                        std::vector<double> c = FSinterface[i] -> getComponentY();
                        constrain[1] = 3;
                        value[1] = c[0];
                    }
                }
            }
            Boundary * bound = new Boundary(connectB, numBoundElems++, constrain, value, ibound, cmesh);
            // std::cout << "asdasd " << rank << " " << ibound << std::endl;
            cmesh->BoundaryVec().push_back(bound);
        }   
    }
    cmesh->part_elem= new int[cmesh->NElements()]();
}


void GmshTools::RenumberConnectivity(CompMesh *cmesh){
    // Renumber nodes - start
    std::vector<int > neighborNodes;

    int* xadj;
    int numNd = cmesh->NNodes();
    std::vector<int> adjncy;
    xadj = new int[numNd+1]();
    adjncy.reserve(10*numNd);

    for (int iNode = 0; iNode < numNd; iNode++){
        neighborNodes.reserve(cmesh->NodeVec()[iNode] -> getNumberOfElements()*3);
        neighborNodes.push_back(iNode);

        for (int j = 0; j < cmesh->NodeVec()[iNode] -> getNumberOfElements(); j++){
            int elem = cmesh->NodeVec()[iNode] -> getInverseIncidenceElement(j);
            VecInt connec = cmesh->ElementVec()[elem-1] -> getConnectivity();

            // std::cout << "COMM " << connec[0] << " " << connec[4] << std::endl;
            bool flag = false;
            for (int i = 0; i < cmesh->NElNodes(); i++){
                for (int iNeig = 0; iNeig < neighborNodes.size(); iNeig++){
                    if (connec[i] == neighborNodes[iNeig]){
                        flag = true;
                        break;
                    }
                }
                if (flag == false) neighborNodes.push_back(connec[i]);
                flag = false;
            }
        }
        //Save nodal adjacency for domain partitioning
        xadj[iNode+1] = xadj[iNode] + neighborNodes.size() - 1;
        for (int i = 1; i < neighborNodes.size(); i++){
            adjncy.push_back(neighborNodes[i]);
        }
        neighborNodes.clear();
        neighborNodes.shrink_to_fit();
    }

    //Save a second adjacency vector in idx_t format
    int *adjncy2;
    int adj_size = adjncy.size();
    adjncy2 = new int[adj_size];
    for (int i = 0; i < adjncy.size(); i++) adjncy2[i] = adjncy[i];
    
    adjncy.clear(); adjncy.shrink_to_fit();

    int* perm;
    int* iperm;
    perm = new int[numNd];
    iperm = new int[numNd];

    // Call METIS for node renumbering

    // MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
    //                         numNodes, numNodes,
    //                         1,NULL,1,NULL,&A); 
    
    // // MatGetOrdering(A, MATORDERINGMETISND, IS *rperm, IS *cperm)

    // MatDestroy(&A);

    METIS_NodeND(&numNd, xadj, adjncy2, NULL, NULL, perm, iperm);

    //Reorder nodes
    for (int i = 0, j; i < cmesh->NNodes(); ++i) {
        for (j = iperm[i]; j < i; j = iperm[j]);
        if (j == i) while (j = iperm[j],j != i) std::swap(cmesh->NodeVec()[i],cmesh->NodeVec()[j]);
    }

    // Update connectivity
    for (int i = 0; i < cmesh->NElements(); i++){
        VecInt connect = cmesh->ElementVec()[i] -> getConnectivity();

        //Reorder connectivity
        for (int k = 0; k < cmesh->NElNodes(); k++) connect[k] = iperm[connect[k]];
        cmesh->ElementVec()[i] -> setConnectivity(connect);
    }
    // Update boundary connectivity
    for (int ibound = 0; ibound < cmesh->NBoundElements(); ibound++){
        VecInt connectB = cmesh->BoundaryVec()[ibound] -> getBoundaryConnectivity();

        for (int k = 0; k < cmesh->NBdNodes(); k++) connectB[k] = iperm[connectB[k]];
        cmesh->BoundaryVec()[ibound] -> setBoundaryConnectivity(connectB);
    }
    
    for (int i = 0; i < cmesh->NNodes(); i++) cmesh->NodeVec()[i] -> clearInverseIncidence();

    for (int i = 0; i < cmesh->NElements(); i++){
        VecInt connect = cmesh->ElementVec()[i] -> getConnectivity();

        for (int k = 0; k < cmesh->NElNodes(); k++) cmesh->NodeVec()[connect[k]] -> pushInverseIncidence(i);
    }

    PetscLogDouble bytes = 0;
    PetscMemoryGetCurrentUsage(&bytes);
    PetscPrintf(PETSC_COMM_WORLD,"Memory used-1 %g M\n",bytes/(1024*1024));
    
    for (int i = 0; i < cmesh->NNodes(); i++){
        for (int j = 0; j < cmesh->NodeVec()[i] -> getNumberOfElements(); j++){
            int elJ = cmesh->NodeVec()[i] -> getInverseIncidenceElement(j);
            for (int k = 0; k <cmesh-> NodeVec()[i] -> getNumberOfElements(); k++)
                cmesh->ElementVec()[elJ] -> pushNeighborElement(cmesh->NodeVec()[i] -> getInverseIncidenceElement(k)); 
        }
    }

    PetscMemoryGetCurrentUsage(&bytes);
    PetscPrintf(PETSC_COMM_WORLD,"Memory used00 %g M\n",bytes/(1024*1024));

    for (int i = 0; i < cmesh->NElements(); i++) cmesh->ElementVec()[i] -> sortEraseNeighborElements();

    PetscMemoryGetCurrentUsage(&bytes);
    PetscPrintf(PETSC_COMM_WORLD,"Memory used11 %g M\n",bytes/(1024*1024));        

    // for (int i = 0; i < NElements(); i++){
    //     if (rank == 0) std::cout << "Neighbor " << i << " ";
    //     for (int j = 0; j < ElementVec()[i] -> getNumberOfNeighborElements(); j++){
    //         std::cout << ElementVec()[i] -> getNeighborElement(j) << " ";
    //     }
    //     std::cout << std::endl;
    // }

    delete [] perm;
    delete [] iperm;
    delete [] adjncy2;
    delete [] xadj;
    // Renumber nodes - end

}

void GmshTools::BoundaryConstrains(CompMesh * cmesh){
    // if (rank == 0) std::cout << "8/9 Setting boundary conditions..." << std::endl;

    for (int ibound = 0; ibound < cmesh->NBoundElements(); ibound++){
        
        VecInt connectB = cmesh->BoundaryVec()[ibound] -> getBoundaryConnectivity();

        for (int j = 0; j < cmesh->NBdNodes(); j++){
            int nstate = cmesh->NodeVec()[connectB[j]]->GetNStateVariables();
            for (int istate = 0; istate < nstate; istate++){
                if ((cmesh->BoundaryVec()[ibound] -> getConstrain(istate) == 1) || (cmesh->BoundaryVec()[ibound] -> getConstrain(istate) == 3)){
                    // for (int j = 0; j < nBdNodes; j++) 
                    //     NodeVec()[connectB[j]] -> setConstrains(k,BoundaryVec()[ibound] -> getConstrain(k),
                    //                                         BoundaryVec()[ibound] -> getConstrainValue(k));
            
                    VecDouble exactSol(nstate);
                    MatrixDouble gradExactSol(cmesh->Dimension(),nstate);
                    auto exact = cmesh->getProblemParameters().getExactSolution();
                    if (exact){
                        VecDouble x = cmesh->NodeVec()[connectB[j]]->getCoordinates();
                        exact(x,exactSol,gradExactSol);
                        
                        cmesh->NodeVec()[connectB[j]] -> SetBoundaryCondition(istate,cmesh->BoundaryVec()[ibound] -> getConstrain(istate),
                                                                exactSol[istate]);
                    } else {
                        cmesh->NodeVec()[connectB[j]] -> SetBoundaryCondition(istate,cmesh->BoundaryVec()[ibound] -> getConstrain(istate),
                                                                              cmesh->BoundaryVec()[ibound] -> getConstrainValue(istate));
                    }
                }
            };
        }
    };
}

void GmshTools::BoundarySides(CompMesh * cmesh){
    // if (rank == 0) std::cout << "9/9 Setting boundary sides..." << std::endl;
    int nBdNodes = cmesh->NBdNodes();
    int DIM = cmesh->Dimension();
    //Sets fluid elements and sides on interface boundaries
    for (int i=0; i<cmesh->NBoundElements(); i++){

        int group = cmesh->BoundaryVec()[i] -> getBoundaryGroup();

       if ((cmesh->BoundaryVec()[i] -> getConstrain(0) > 0) || (cmesh->BoundaryVec()[i] -> getConstrain(1) > 0)) {
            
            VecInt connectB = cmesh->BoundaryVec()[i] -> getBoundaryConnectivity();

            for (int j=0; j<cmesh->NElements(); j++){
                VecInt connect = cmesh->ElementVec()[j] -> getConnectivity();

                int flag = 0;
            
                int side[nBdNodes];
                for (int k=0; k<cmesh->NElNodes(); k++){
                    for (int l = 0; l<nBdNodes; l++){
                        if (connectB[l] == connect[k]){
                            side[flag] = k;
                            flag++;
                        }
                    };
                };
                if (flag == nBdNodes){
                    cmesh->BoundaryVec()[i] -> setElement(j);
                    //Sets element index and side
                    // for (int k=0; k<nBdNodes; k++) std::cout << "BD NODES " << i << " " << j << " " << k << " " << side[k] << std::endl;
                    
                    for (int k=0; k<DIM+1; k++){
                        // std::cout << "DDDDD " << k << std::endl;

                        int* end = side + nBdNodes;
                        int* foo = std::find(side, end, k);

                        if ((foo == end) && (cmesh->ElementVec()[j] -> getElemSideInBoundary() < 0)){
                            int aux = cmesh->BoundaryVec()[i] -> getBoundaryGroup();
                            cmesh->BoundaryVec()[i] -> setBoundaryGroup(aux);
                            cmesh->BoundaryVec()[i] -> setElementSide(k);
                            cmesh->BoundaryVec()[i] -> setElement(j);
                            cmesh->ElementVec()[j] -> setElemSideInBoundary(k);
                            if (cmesh->BoundaryVec()[i] -> getConstrain(0) == 3) cmesh->ElementVec()[j] -> setFSIInterface();
                        }
                    }
                };
            };
        }
    };

}