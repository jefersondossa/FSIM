#include "solidDomain.h"

int main()
{
	Geometry* solid = new Geometry(0);

	Point* p0 = solid->addPoint({ 0.0, 0.0 });
	Point* p1 = solid->addPoint({ 0.5, 0.0 });
	Point* p2 = solid->addPoint({ 0.5, 0.15 });
	Point* p3 = solid->addPoint({ 0.0, 0.35 });

	Line* l0 = solid->addLine({ p0, p1 });
	Line* l1 = solid->addLine({ p1, p2 });
	Line* l2 = solid->addLine({ p2, p3 });
	Line* l3 = solid->addLine({ p3, p0 });

	//LineLoop* ll0 = solid->addLineLoop({ l0, l1, l2, l3 });

	PlaneSurface* s0 = solid->addPlaneSurface({ l0, l1, l2, l3 });
	solid->transfiniteLine({ l0 }, 51);
	solid->transfiniteLine({ l1 }, 16);
	solid->transfiniteLine({ l2 }, 55);
	solid->transfiniteLine({ l3 }, 36);
	//solid->transfiniteSurface({ s0 }, "Left", {p0,p1,p2,p3});

	solid->addBoundaryCondition("DIRICHLET", l0, {}, { 0.0 }, "GLOBAL");
	solid->addBoundaryCondition("DIRICHLET", l1, {0.0}, {}, "GLOBAL");
	solid->addBoundaryCondition("DIRICHLET", l3, {0.0}, {}, "GLOBAL");
	
	SolidDomain* problem = new SolidDomain(solid);
	problem->addSurfaceMaterial({ s0 }, 1.0, 1.0, 1.0, "PLANE_STRESS");
	problem->generateMesh("T3", "DELAUNAY", "exemplo", "", true, true);


	return 0;

}