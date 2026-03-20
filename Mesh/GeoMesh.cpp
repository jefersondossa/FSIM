#include "GeoMesh.h"

void GeoMesh::Print(std::string filename){
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    file << "Geometric Mesh Information\n";
    file << "===============================\n";
    file << "Number of Nodes: " << NNodes() << "\n";
    file << "Number of Elements: " << NElements() << "\n";
    file << "--------------------------------\n";

    file << "--------------------------------\n";
    file << "Nodes Information:\n";
    for (int64_t i = 0; i < NNodes(); i++)
    {
        Node *node = fNodeVector[i];
        file << "Node " << i << ": ";
        file << "Coordinates = [";
        VecDouble coord = node->getCoordinates();
        for (size_t j = 0; j < coord.size(); j++) {
            file << coord[j];
            if (j < coord.size() - 1) file << ", ";
        }
        file << "]\n";
    }


    file << "--------------------------------\n";
    file << "Elements Information:\n";
    for (int64_t i = 0; i < NElements(); i++)
    {
        GeoElement *el = fElementVector[i];
        file << "Element " << i << ": ";
        file << "Material ID = " << el->Material() << ", ";
        file << "Element Type = " << el->PrintType() << ", ";
        file << "Geometric Nodes = [";
        VecInt &geoNodes = el->getGeometricNodes();
        for (size_t j = 0; j < geoNodes.size(); j++) {
            file << geoNodes[j];
            if (j < geoNodes.size() - 1) file << ", ";
        }
        file << "]\n";
    }      

    file.close();
}