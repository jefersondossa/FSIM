#include "Element.h"
#include <cassert>

VecDouble Element::GetGeometricCenter(CompMesh& mesh)
{
    VecDouble res{3};
    res.setZero();

    const auto number_of_nodes = NElNodes();

    VecInt connec = getConnectivity();

    for (int i = 0; i < connec.size(); i++){
        const auto i_node = connec[i];
        assert(i_node < mesh.NodeVec().size());
        res += (mesh.NodeVec()[i_node]->getCoordinates()) * (1.0 / number_of_nodes);
    }

    return res;
}