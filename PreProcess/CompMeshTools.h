#ifndef COMPMESHTOOLS_H
#define COMPMESHTOOLS_H

#include "Element.h"
#include "CompMesh.h"
class CompMesh;


namespace CompMeshTools{
    /// @brief Performs domain decomposition, used for MPI parallelism
    /// @param cmesh Computational mesh
    void DomainDecompositionMETIS(CompMesh* cmesh); 
}


#endif