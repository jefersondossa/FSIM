#include "PlaneSurface.h"




PlaneSurface::PlaneSurface() {}

PlaneSurface::PlaneSurface(const int& index, const std::string& name, LineLoop* lineLoop)
	: Surface(index, name, lineLoop) {}

PlaneSurface::PlaneSurface(const int& index, const std::string& name, std::vector<LineLoop*> lineLoops)
	: Surface(index, name, lineLoops) {}

PlaneSurface::~PlaneSurface() {}

std::string PlaneSurface::getGmshCode()
{
	std::stringstream text;
	if (lineLoop_ != nullptr) {
		text << name_ << " = news; Plane Surface(" << name_ << ") = {" << lineLoop_->getName() << "}; Physical Surface('" << name_ << "') = {" << name_ << "};\n//\n";
	} else if (lineLoops_.size() != 0){
		text << name_ << " = news; Plane Surface(" << name_ << ") = {";
		for (size_t i = 0; i < lineLoops_.size(); i++)
		{
			text << lineLoops_[i]->getName();
			if (i != (lineLoops_.size() - 1))
				text << ", ";
		}
		text << "}; Physical Surface('" << name_ << "') = {" << name_ << "};\n//\n";
	}
	// text << name_ << " = news; Plane Surface(" << name_ << ") = {" << lineLoop_->getName() << "}; Physical Surface('" << name_ << "') = {" << name_ << "};\n//\n";
	return text.str();
}

