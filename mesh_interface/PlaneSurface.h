#ifndef PLANE_SURFACE_H
#define PLANE_SURFACE_H

#include "Surface.h"

class PlaneSurface : public Surface
{
public:
	PlaneSurface();

	PlaneSurface(const int& index, const std::string& name, LineLoop* lineLoop);

	~PlaneSurface();

	std::string getGmshCode() override;
};

#endif

