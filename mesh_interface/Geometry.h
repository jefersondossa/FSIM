#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "Volume.h"
#include "GeometricDirichlet.h"
#include "GeometricNeumann.h"
#include "BoundaryCondition.h"
#include <unordered_map>
#include "Ellipse.h"

class Geometry
{
public:
	Geometry();

	Geometry(const int& index);

	~Geometry();

	int getIndex();

	int getNumberOfPoints();

	int getNumberOfLines();

	int getNumberOfLineLoops();

	int getNumberOfSurfaces();

	int getNumberOfSurfaceLoops();

	int getNumberOfVolumes();

	int getNumberOfBoundaryConditions(const std::string& type);

	Point* getPoint(const std::string& name);

	Line* getLine(const std::string& name);

	LineLoop* getLineLoop(const std::string& name);

	Surface* getSurface(const std::string& name);

	SurfaceLoop* getSurfaceLoop(const std::string& name);

	Volume* getVolume(const std::string& name);

	std::unordered_map<std::string, Line*> getLines();

	std::unordered_map<std::string, Surface*> getSurfaces();

	std::unordered_map<std::string, Volume*> getVolumes();

	std::vector<GeometricDirichlet*> getDirichletBoundaryConditions();

	std::vector<GeometricNeumann*> getNeumannBoundaryConditions();

    std::vector<BoundaryCondition*> getBoundaryCondition(const std::string& type);

	std::string getGmshCode();

	Point* addPoint(std::vector<double> coordinates, const double& lcar = 1.0, const bool& discretization = true);

	Line* addLine(std::vector<Point*> points, const bool& discretization = true);

	Circle* addCircle(std::vector<Point*> points, const bool& discretization = true);

	Ellipse* addEllipse(std::vector<Point*> points, double rx, double ry, double angle = 0, const bool& discretization = true);

	Spline* addSpline(std::vector<Point*> points, double(*function)(double), const int& ndiv, const bool& discretization = true);

	LineLoop* addLineLoop(std::vector<Line*> lines);

	Surface* addSurface(LineLoop* lineLoop);

	Surface* addSurface(std::vector<Line*> lines);

	PlaneSurface* addPlaneSurface(LineLoop* lineLoop);

	PlaneSurface* addPlaneSurface(std::vector<Line*> lines);

	SurfaceLoop* addSurfaceLoop(std::vector<Surface*> surfaces);

	Volume* addVolume(SurfaceLoop* surfaceLoop);

	Volume* addVolume(std::vector<Surface*> surfaces);

	void appendGmshCode(std::string text);

	void transfiniteLine(std::vector<Line*> lines, const int& divisions, const double& progression = 1);

	void transfiniteSurface(std::vector<Surface*> surfaces, std::string oientation = "Left", std::vector<Point*> points = std::vector<Point*>());

	void transfiniteVolume(std::vector<Volume*> volumes);

	void addDirichletBoundaryCondition(const std::vector<Point*>& points, const ConstrainedDOF dof, const double value);
	
	void addDirichletBoundaryCondition(const std::vector<Line*>& lines, const ConstrainedDOF dof, const double value);

	void addDirichletBoundaryCondition(const std::vector<Surface*>& surfaces, const ConstrainedDOF dof, const double value);

	void addDirichletBoundaryCondition(const std::vector<Volume*>& volumes, const ConstrainedDOF dof, const double value);

	void addNeumannBoundaryCondition(const std::vector<Point*>& points, const double valueX, const double valueY, const double valueZ);
	
	void addNeumannBoundaryCondition(const std::vector<Line*>& lines, const double valueX, const double valueY, const double valueZ);

	void addNeumannBoundaryCondition(const std::vector<Surface*>& surfaces, const double valueX, const double valueY, const double valueZ);


	void addBoundaryCondition(const std::string& type, Point* point, const std::vector<double>& componentX = std::vector<double>(), const std::vector<double>& componentY = std::vector<double>(), const std::vector<double>& componentZ = std::vector<double>(), const std::string& referenceSystem = "GLOBAL", const std::string& method = "STRONG", const double& penaltyParameter = 1.0e6);
	
	void addBoundaryCondition(const std::string& type, Line* line, const std::vector<double>& componentX = std::vector<double>(), const std::vector<double>& componentY = std::vector<double>(), const std::vector<double>& componentZ = std::vector<double>(), const std::string& referenceSystem = "GLOBAL", const std::string& method = "STRONG", const double& penaltyParameter = 1.0e6);

	void addBoundaryCondition(const std::string& type, Surface* surface, const std::vector<double>& componentX = std::vector<double>(), const std::vector<double>& componentY = std::vector<double>(), const std::vector<double>& componentZ = std::vector<double>(), const std::string& referenceSystem = "GLOBAL", const std::string& method = "STRONG", const double& penaltyParameter = 1.0e6);


private:
	int index_;
	std::unordered_map<std::string, Point*> points_;
	std::unordered_map<std::string, Line*> lines_;
	std::unordered_map<std::string, LineLoop*> lineLoops_;
	std::unordered_map<std::string, Surface*> surfaces_;
	std::unordered_map<std::string, SurfaceLoop*> surfaceLoops_;
	std::unordered_map<std::string, Volume*> volumes_;
	std::vector<GeometricDirichlet*> dirichlet_;
	std::vector<GeometricNeumann*> neumann_;
	std::unordered_map<std::string, std::vector<BoundaryCondition*>> boundaryConditions_;
	std::string gmshCode_;
};


#endif