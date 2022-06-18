#include "ElementMesh.h"

///----------------------------------------------------------------------------
///-------------------------------IMPLEMENTATION-------------------------------
///----------------------------------------------------------------------------

ElementMesh::ElementMesh(const int& index, const std::vector<NodeMesh*>& nodes, Material* material, const double& thickness)
{
	index_ = index;
	nodes_ = nodes;
	material_ = material;
	thickness_ = thickness;
}

ElementMesh::~ElementMesh() {}
