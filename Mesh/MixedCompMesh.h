#ifndef MixedCompMesh_H
#define MixedCompMesh_H

#include "CompMesh.h"

class MixedCompMesh : public CompMesh
{
private:
    int fNSpaces = 1;

    std::vector<CompMesh* > fMeshVector;

public:
    MixedCompMesh(std::vector<CompMesh* > meshvector);

    ~MixedCompMesh() = default;

    void AutoBuild() override;

    void CreateMixedConnects();


};

#endif