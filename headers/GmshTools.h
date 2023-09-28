#ifndef GMSHTOOLS_H
#define GMSHTOOLS_H

#include "Geometry.h"
#include "Element.h"
#include "Boundary.h"
#include "CompMesh.h"
class CompMesh;

namespace GmshTools{

    void MeshReading(Geometry* &geometry_, const std::string& mshfile, CompMesh* cmesh);
    std::vector<std::string> SplitLine(std::string str, std::string delim);
    void ReadNodes(std::ifstream &file, CompMesh * cmesh);
    void ReadElements(Geometry* &geometry_, std::ifstream &file, std::unordered_map<int, std::string> &physicalEntities, CompMesh * cmesh);
    void RenumberConnectivity(CompMesh * cmesh);
    void BoundaryConstrains(CompMesh * cmesh);
    void BoundarySides(CompMesh * cmesh);
};


#endif