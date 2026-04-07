#include "Ellipse.h"

Ellipse::Ellipse(const int& index, const std::string& name, std::vector<Point*> points, double rx, double ry, double angle, const bool& discretization)
{
    index_ = index;
	name_ = name;
	discretization_ = discretization;
	points_.reserve(3);
	for (Point* point : points)
		points_.push_back(point);
	rx_ = rx;
	ry_ = ry;
	angle_ = angle;
}

Ellipse::~Ellipse() {}

Ellipse* Ellipse::operator-()
{
	Ellipse* copy = new Ellipse(index_, name_, {points_}, rx_, ry_, angle_, discretization_);
	copy->setName("-" + name_);
	return copy;
}

std::string Ellipse::getGmshCode()
{
    std::stringstream text;
	if (discretization_) {
		text << name_ << " = newl; Ellipse(" << name_ << ") = {" << points_[0]->getX() << ", " << points_[0]->getY() << ", " << points_[0]->getZ()
				<<  ", " << rx_ << ", " << ry_ << ", 0, 2*Pi}; Physical Line('" << name_ << "') = {" << name_ << "};\n//\n";
		text << "Rotate {" << points_[0]->getX() << ", " << points_[0]->getY() << ", " << points_[0]->getZ() << "}, {0, 0, 1}, " << angle_ << "*Pi/180} {" << name_ << "};\n//\n";
		return text.str();
	}
	else {
		text << name_ << " = newl; Ellipse(" << name_ << ") = {" << points_[0]->getX() << ", " << points_[0]->getY() << ", " << points_[0]->getZ()
				<<  ", " << rx_ << ", " << ry_ << ", 0, 2*Pi}; Physical Line('" << name_ << "') = {" << name_ << "};\n//\n";
		text << "Rotate {" << points_[0]->getX() << ", " << points_[0]->getY() << ", " << points_[0]->getZ() << "}, {0, 0, 1}, " << angle_ << "*Pi/180} {" << name_ << "};\n//\n";
		return text.str();
	}
}

Point* Ellipse::getMiddlePoint()
{
    return points_[1];
}
