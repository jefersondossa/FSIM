#include "GeoElement.h"


VecDouble GeoElement::GetGeometricCenter()
{
    VecDouble res{3};
    res.setZero();

    const auto number_of_nodes = fGeoNodes.size();

    for (int i = 0; i < fGeoNodes.size(); i++){
        const auto i_node = fGeoNodes[i];
        assert(i_node < fMesh->NodeVec().size());
        res += (fMesh->NodeVec()[i_node]->getCoordinates()) * (1.0 / number_of_nodes);
    }

    return res;
}

