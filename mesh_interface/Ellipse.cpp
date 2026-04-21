#include "Ellipse.h"

Ellipse::Ellipse(const int& index, const std::string& name, std::vector<Point*> points, const bool& discretization)
{
    index_ = index;
	name_ = name;
	discretization_ = discretization;
	points_.reserve(3);
	for (Point* point : points)
		points_.push_back(point);
}

Ellipse::~Ellipse() {}

Ellipse* Ellipse::operator-()
{
	Ellipse* copy = new Ellipse(index_, name_, {points_}, discretization_);
	copy->setName("-" + name_);
	return copy;
}

std::string Ellipse::getGmshCode()
{
    std::stringstream text;
	if (discretization_) {
		text << name_ << " = newl; Ellipse(" << name_ << ") = {" << 
		points_[0]->getName() << ", " << points_[1]->getName() << ", " << 
		points_[2]->getName() <<  ", " << points_[3]->getName() << 
		"}; Physical Line('" << name_ << "') = {" << name_ << "};\n//\n";
		return text.str();
	}
	else {
		text << name_ << " = newl; Ellipse(" << name_ << ") = {" << 
		points_[0]->getName() << ", " << points_[1]->getName() << ", " << 
		points_[2]->getName() <<  ", " << points_[3]->getName() << 
		"}; Physical Line('" << name_ << "') = {" << name_ << "};\n//\n";
		return text.str();
	}
}

Point* Ellipse::getMiddlePoint()
{
    return points_[1];
}
