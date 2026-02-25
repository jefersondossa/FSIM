#ifndef MixedCompMesh_H
#define MixedCompMesh_H

#include "CompMesh.h"

class MixedCompMesh : public CompMesh
{
private:
    int fNSpaces = 1;

public:
    MixedCompMesh() = default;

    MixedCompMesh(ApproxType approxType, int nspaces) : CompMesh(approxType) {
        fNSpaces = nspaces;
    };

    ~MixedCompMesh() = default;

    


};

#endif