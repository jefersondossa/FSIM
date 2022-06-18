#ifndef NODE_MESH_H
#define NODE_MESH_H

#include <vector>
#include <boost/numeric/ublas/vector.hpp>

using namespace boost::numeric::ublas;

class NodeMesh
{
public:
	NodeMesh(const int& index, double* coord);

	~NodeMesh();

	int getIndex();

	double* getCoordinates();

	double getX();

	double getY();

	double getZ();

	void setIndex(const int& index);

	void setCoordinates(double* coord);

	void setX(const double& x);

	void setY(const double& y);

	void setZ(const double& z);

private:
	int index_;
	double coord_[3];
};

#endif
