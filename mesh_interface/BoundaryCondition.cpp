#include "BoundaryCondition.h"



///----------------------------------------------------------------------------
///-------------------------------IMPLEMENTATION-------------------------------
///----------------------------------------------------------------------------

BoundaryCondition::BoundaryCondition(const int& index,
	const std::string& object,
	const std::vector<double>& componentX,
	const std::vector<double>& componentY,
	const std::vector<double>& componentZ,
	const std::string& referenceSystem,
	const std::string& method,
	const double& penaltyParameter)
{
	index_ = index;
	(object[0] == 'p') ? point_ = object : line_ = object;
	componentX_ = componentX;
	componentY_ = componentY;
	componentZ_ = componentZ;
	referenceSystem_ = referenceSystem;
	method_ = method;
	penaltyParameter_ = penaltyParameter;
}

BoundaryCondition::~BoundaryCondition() {}

int BoundaryCondition::getIndex()
{
	return index_;
}

std::string BoundaryCondition::getPointName()
{
	return point_;
}

std::string BoundaryCondition::getLineName()
{
	return line_;
}

std::string BoundaryCondition::getReferenceSystem()
{
	return referenceSystem_;
}

std::vector<double> BoundaryCondition::getComponentX()
{
	return componentX_;
}

std::vector<double> BoundaryCondition::getComponentY()
{
	return componentY_;
}

std::vector<double> BoundaryCondition::getComponentZ()
{
	return componentZ_;
}

std::string BoundaryCondition::getMethod()
{
	return method_;
}

double BoundaryCondition::getPenaltyParameter()
{
	return penaltyParameter_;
}

void BoundaryCondition::setIndex(const int& index)
{
	index_ = index;
}

void BoundaryCondition::setPointName(const std::string& name)
{
	point_ = name;
}

void BoundaryCondition::setLineName(const std::string& name)
{
	line_ = name;
}

void BoundaryCondition::setReferenceSystem(const std::string& referenceSystem)
{
	referenceSystem_ = referenceSystem;
}

void BoundaryCondition::setComponentX(std::vector<double> componentX)
{
	componentX_ = componentX;
}

void BoundaryCondition::setComponentY(std::vector<double> componentY)
{
	componentY_ = componentY;
}

void BoundaryCondition::setComponentZ(std::vector<double> componentZ)
{
	componentZ_ = componentZ;
}

void BoundaryCondition::setMethod(const std::string& method)
{
	method_ = method;
}

void BoundaryCondition::setPenaltyParameter(const double& penaltyParameter)
{
	penaltyParameter_ = penaltyParameter;
}