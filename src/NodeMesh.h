#pragma once

class NodeMesh
{
public:
	NodeMesh(const int& index, double* coord);

	~NodeMesh();

	int getIndex();

	double* getCoordinates();

	double getX();

	double getY();

	void setIndex(const int& index);

	void setCoordinates(double* coord);

	void setX(const double& x);

	void setY(const double& y);

private:
	int index_;
	double coord_[2];
};

///----------------------------------------------------------------------------
///-------------------------------IMPLEMENTATION-------------------------------
///----------------------------------------------------------------------------

NodeMesh::NodeMesh(const int& index, double* coord)
{
	index_ = index;
	coord_[0] = coord[0];
	coord_[1] = coord[1];
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