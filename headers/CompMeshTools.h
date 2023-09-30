#ifndef COMPMESHTOOLS_H
#define COMPMESHTOOLS_H

#include "Geometry.h"
#include "Element.h"
#include "Boundary.h"
#include "CompMesh.h"
class CompMesh;


namespace CompMeshTools{
    void DomainDecompositionMETIS(CompMesh* cmesh); 
    

    
}


#endif