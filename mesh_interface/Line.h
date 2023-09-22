#ifndef LINE_H
#define LINE_H

#include "Point.h"
#include "ElementMesh.h"
#include "Material.h"

class Line
{
public:
	Line();

	Line(const int& index, const std::string& name, std::vector<Point*> points, const bool& discretization = true);

	~Line();

	Line* operator-();

	int getIndex();

	std::string getName();

	Point* getInitialPoint();

	Point* getEndPoint();

	bool getDiscretization();

	std::vector<NodeMesh*> NodeVec();

	// std::vector<ParametricLineElement*> getParametricElements();

	std::vector<ElementMesh*> getElements();

	int getNumberOfNodes() const;

	// int getNumberOfParametricElements() const;

	int getNumberOfElements() const;

	Material* getMaterial();

	std::string virtual getGmshCode();

	void setIndex(const int& index);

	void setName(const std::string& name);

	void setInitialPoint(Point& point);

	void setEndPoint(Point& point);

	void setDiscretization(const bool& discretization);

	void addNodes(const std::vector<NodeMesh*>& nodes);

	// void addParametricElement(ParametricLineElement* element);

	void addElement(ElementMesh* element);

	void setMaterial(Material* material);

protected:
	int index_;
	std::string name_;
	std::vector<Point*> points_;
	bool discretization_;
	std::vector<NodeMesh*> nodes_;
	// std::vector<ParametricLineElement*> parametricElements_;
	std::vector<ElementMesh*> elements_;
	Material* material_;
};

#endif


