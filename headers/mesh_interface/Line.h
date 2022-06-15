#pragma once

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

	std::vector<NodeMesh*> getNodes();

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







Line::Line() {}

Line::Line(const int& index, const std::string& name, std::vector<Point*> points, const bool& discretization)
{
	index_ = index;
	name_ = name;
	discretization_ = discretization;
	points_.reserve(2);
	for (Point* point : points)
		points_.push_back(point);
}

Line::~Line() {}

Line* Line::operator-()
{
	Line* copy = new Line(index_, name_, {points_});
	copy->setName("-" + name_);
	return copy;
}

int Line::getIndex()
{
	return index_;
}
std::string Line::getName()
{
	return name_;
}
Point* Line::getInitialPoint()
{
	return points_[0];
}
Point* Line::getEndPoint()
{
	int last = points_.size() -1;
	return points_[last];
}
bool Line::getDiscretization()
{
	return discretization_;
}
std::vector<NodeMesh*> Line::getNodes()
{
	return nodes_;
}

// std::vector<ParametricLineElement*> Line::getParametricElements()
// {
// 	return parametricElements_;
// }

std::vector<ElementMesh*> Line::getElements()
{
	return elements_;
}

int Line::getNumberOfNodes() const
{
	return nodes_.size();
}

// int Line::getNumberOfParametricElements() const
// {
// 	return parametricElements_.size();
// }

int Line::getNumberOfElements() const
{
	return elements_.size();
}

Material* Line::getMaterial()
{
	return material_;
}

std::string Line::getGmshCode()
{
	std::stringstream text;
	if (discretization_) {
		text << name_ << " = newl; Line(" << name_ << ") = {" << points_[0]->getName() ;
		for (int i = 1; i < points_.size(); ++i){
			text <<  ", " << points_[i]->getName();	
		}
		text << "}; Physical Line('" << name_ << "') = {" << name_ << "};\n//\n";
		return text.str();
	}
	else {
		text << name_ << " = newl; Line(" << name_ << ") = {" << points_[0]->getName() << ", " << points_[1]->getName()
			<< "};\n//\n";
		return text.str();
	}
}
void Line::setIndex(const int& index)
{
	index_ = index;
}
void Line::setName(const std::string& name)
{
	name_ = name;
}
void Line::setInitialPoint(Point& point)
{
	//points_[0] = point;
}
void Line::setEndPoint(Point& point)
{
	//points_[1] = point;
}
void Line::setDiscretization(const bool& discretization)
{
	discretization_ = discretization;
}
void Line::setMaterial(Material* const material)
{
	material_ = material;
}
void Line::addNodes(const std::vector<NodeMesh*>& nodes)
{
	for (NodeMesh* node1 : nodes)
	{
		bool notDuplicate = true;
		for (NodeMesh* node2 : nodes_)
		{
			if (node1->getIndex() == node2->getIndex())
			{
				notDuplicate = false;
				break;
			}
		}
		if (notDuplicate)
			nodes_.push_back(node1);
	}
}
// void Line::addParametricElement(ParametricLineElement* element)
// {
// 	parametricElements_.push_back(element);
// }
void Line::addElement(ElementMesh* element)
{
	elements_.push_back(element);
}














