#ifndef COMPMESHTOOLS_H
#define COMPMESHTOOLS_H

#include "Element.h"
#include "CompMesh.h"

class CompMesh;
class Element;

namespace CompMeshTools{
    /// @brief Performs domain decomposition, used for MPI parallelism
    /// @param cmesh Computational mesh
    void DomainDecompositionMETIS(CompMesh* cmesh); 

    /// @brief Searchs the node correspondence in a geometric mesh
    /// @param x point
    /// @param cmesh the geometric mesh
    /// @param elCorr correpondence element
    /// @param xsiCorr correspondence Xsi
    /// @param elemsearch the first guess element
    bool searchNodeCorrespondence(VecDouble &x, CompMesh *cmesh, int &elCorr, VecDouble &xsiCorr, Element *elemsearch);
}


#endif