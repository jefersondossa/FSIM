#include "NodeMesh.h"


///----------------------------------------------------------------------------
///-------------------------------IMPLEMENTATION-------------------------------
///----------------------------------------------------------------------------

NodeMesh::NodeMesh(const int& index, double* coord)
{
	index_ = index;
	for(int i=0; i<3; i++) coord_[i] = coord[i];
}

NodeMesh::~NodeMesh() {}

int NodeMesh::getIndex()
{
	return index_;
}

double NodeMesh::getX()
{
	return coord_[0];
}

double NodeMesh::getY()
{
	return coord_[1];
}

double NodeMesh::getZ()
{
	return coord_[2];
}

double* NodeMesh::getCoordinates()
{
	return coord_;
}

void NodeMesh::setIndex(const int& index)
{
	index_ = index;
}

void NodeMesh::setX(const double& x)
{
	coord_[0] = x;
}

void NodeMesh::setY(const double& y)
{
	coord_[1] = y;
}

void NodeMesh::setZ(const double& z)
{
	coord_[2] = z;
}