#include "CompMeshTools.h"
#include "fstream"

#ifdef HAS_METIS
#include <metis.h>
#endif

void CompMeshTools::DomainDecompositionMETIS(CompMesh *cmesh){

    PanicButton();//Need refactor
    // std::string mirror2;
    // mirror2 = "domain_decomposition.txt";
    // std::ofstream mirrorData(mirror2.c_str());
    
    // int size;

    // MPI_Comm_size(PETSC_COMM_WORLD, &size);

    // idx_t objval;
    // idx_t numEl = cmesh->NElements();
    // idx_t numNd = cmesh->NNodes();
    // idx_t ssize = size;
    // idx_t one = 1;
    // int nElNodes = cmesh->NElNodes();
    // idx_t elem_start[numEl+1], elem_connec[nElNodes*numEl];
    // cmesh->part_elem = new int[numEl];
    // cmesh->part_nodes = new int[numNd];


    // for (int i = 0; i < numEl+1; i++){
    //     elem_start[i]=nElNodes*i;
    // };
    // for (int jel = 0; jel < numEl; jel++){
    //     auto connec=cmesh->ElementVec()[jel]->getConnectivity();        
        
    //     for (int i=0; i<nElNodes; i++){
    //     elem_connec[nElNodes*jel+i] = connec[i];
    //     };
    // };

    // //Performs the domain decomposition
    // if (size == 1){
    //     for (int i = 0; i < numNd; i++) cmesh->part_nodes[i] = 0;
    //     for (int i = 0; i < numEl; i++) cmesh->part_elem[i] = 0;
    // } else {
    //     METIS_PartMeshDual(&numEl, &numNd, elem_start, elem_connec, \
    //                             NULL, NULL, &one, &ssize, NULL, NULL,    \
    //                             &objval, cmesh->part_elem, cmesh->part_nodes);
    // }
    
    // mirrorData << std::endl \
    //            << "FLUID MESH DOMAIN DECOMPOSITION - ELEMENTS" << std::endl;
    // for(int i = 0; i < numEl; i++){
    //     mirrorData << "process = " << cmesh->part_elem[i] \
    //                << ", element = " << i << std::endl;
    // };

    // mirrorData << std::endl \
    //            << "FLUID MESH DOMAIN DECOMPOSITION - NODES" << std::endl;
    // for(int i = 0; i < numNd; i++){
    //     mirrorData << "process = " << cmesh->part_nodes[i] \
    //                << ", node = " << i << std::endl;
    // };


    
}