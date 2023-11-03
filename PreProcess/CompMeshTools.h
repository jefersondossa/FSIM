#ifndef COMPMESHTOOLS_H
#define COMPMESHTOOLS_H

#include "Element.h"
#include "CompMesh.h"
class CompMesh;


namespace CompMeshTools{
    void DomainDecompositionMETIS(CompMesh* cmesh); 
    void InitialSolution(CompMesh* cmesh);
}


#endif