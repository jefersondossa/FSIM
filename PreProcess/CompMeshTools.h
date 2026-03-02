#ifndef COMPMESHTOOLS_H
#define COMPMESHTOOLS_H

#include "Element.h"
#include "CompMesh.h"
class CompMesh;


namespace CompMeshTools{
    /// @brief Performs domain decomposition, used for MPI parallelism
    /// @param cmesh Computational mesh
    void DomainDecompositionMETIS(CompMesh* cmesh); 

    /// @brief Searchs the node correspondence in a geometric mesh
    /// @param x point
    /// @param cmesh the geometric mesh
    /// @param elCorr correpondence element
    /// @param xsiCorr correspondence Xsi
    /// @param elSearch the first guess element
    void searchNodeCorrespondence(VecDouble &x, CompMesh *cmesh, int64_t &elCorr, VecDouble &xsiCorr, int64_t elSearch = 0);
}


#endif