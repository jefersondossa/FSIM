#pragma once

#include "Point.h"
#include "ElementMesh.h"

class ArcCircle
{
public:
	ArcCircle();

	ArcCircle(const int& index, const std::string& name, std::vector<Point*> points, const bool& discretization = true);

	~ArcCircle();

	ArcCircle* operator-();

	int getIndex();

	std::string getName();

	Point* getInitialPoint();

	Point* getEndPoint();

	bool getDiscretization();

	std::vector<NodeMesh*> getLineNodes();

	std::string getGmshCode();

	void setIndex(const int& index);

	void setName(const std::string& name);

	void setInitialPoint(Point& point);

	void setEndPoint(Point& point);

	void setDiscretization(const bool& discretization);

	void addNodesToLine(const std::vector<NodeMesh*>& nodes);

	void addElementsToLine(ElementMesh* element);

private:
	int index_;
	std::string name_;
	std::vector<Point*> points_;
	bool discretization_;
	std::vector<NodeMesh*> lineNodes_;
};

///----------------------------------------------------------------------------
///-------------------------------IMPLEMENTATION-------------------------------
///----------------------------------------------------------------------------

ArcCircle::ArcCircle() {}

ArcCircle::ArcCircle(const int& index, const std::string& name, std::vector<Point*> points, const bool& discretization)
{
	index_ = index;
	name_ = name;
	discretization_ = discretization;
	points_.reserve(2);
		for (Point* point : points)
			points_.push_back(point);
}

ArcCircle::~ArcCircle() {}

ArcCircle* ArcCircle::operator-()
{
//    if(name_[0] == '-') {
//	Line* copy = this;
//	copy->setName(name_.erase(0));
//	return *copy;
//    } else {
//	Line* copy = this;
//	copy->setName(name_.insert(0, "-"));
//	return *copy;
//    }
	ArcCircle* copy = new ArcCircle(index_, name_, {points_});
	copy->setName("-" + name_);
	return copy;
}

int ArcCircle::getIndex()
{
	return index_;
}
std::string ArcCircle::getName()
{
	return name_;
}
Point* ArcCircle::getInitialPoint()
{
	return points_[0];
}
Point* ArcCircle::getEndPoint()
{
	return points_[1];
}
bool ArcCircle::getDiscretization()
{
	return discretization_;
}
std::vector<NodeMesh*> ArcCircle::getLineNodes()
{
	return lineNodes_;
}
std::string ArcCircle::getGmshCode()
{
	std::stringstream text;
	if (discretization_) {
		text << name_ << " = newl; Circle(" << name_ << ") = {" << points_[0]->getName() << ", " << points_[1]->getName() << ", " << points_[2]->getName()
			<< "}; Physical Line('" << name_ << "') = {" << name_ << "};\n//\n";
		return text.str();
	}
	else {
		text << name_ << " = newl; Circle(" << name_ << ") = {" << points_[0]->getName() << ", " << points_[1]->getName() << ", " << points_[2]->getName()
			<< "};\n//\n";
		return text.str();
	}
}
void ArcCircle::setIndex(const int& index)
{
	index_ = index;
}
void ArcCircle::setName(const std::string& name)
{
	name_ = name;
}
void ArcCircle::setInitialPoint(Point& point)
{
	//points_[0] = point;
}
void ArcCircle::setEndPoint(Point& point)
{
	//points_[1] = point;
}
void ArcCircle::setDiscretization(const bool& discretization)
{
	discretization_ = discretization;
}
void ArcCircle::addNodesToLine(const std::vector<NodeMesh*>& nodes)
{
	for (NodeMesh* node : nodes)
	lineNodes_.push_back(node);
}
void ArcCircle::addElementsToLine(ElementMesh* element)
{

}