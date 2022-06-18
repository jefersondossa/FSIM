#ifndef FLUID_DOMAIN_H
#define FLUID_DOMAIN_H

#include "Geometry.h"
#include "Mesh.h"
#include "DirichletBoundaryCondition.h"
#include "NeumannBoundaryCondition.h"
#include "BoundaryCondition.h"
#include <unordered_map>
#include <petscksp.h>
#include <metis.h>
#include <algorithm>

class FluidDomain
{
public:
	FluidDomain(Geometry* geometry, const int& index = 0);

	~FluidDomain();

	void setNumberOfSteps(const int numberOfSteps);

	void setMaxNonlinearIterations(const int maxNonlinearIterations);

	void setNonlinearTolerance(const double nonlinearTolerance);

	void setDeltat(const double deltat);

	void setGravity(const double x, const double y, const double z);

	void setSpectralRadius(const double rhoInf);

	void setAlpha(const double& alpha);

	void setMeshLength(const double& h);

	void setInitialAcceleration(const bool initialAccel);

	// void setReferenceConfiguration(const ReferenceConfiguration reference);

	void applyMaterial(const std::vector<Line*> lines, Material*& material);

	void applyMaterial(const std::vector<Surface*> surfaces, Material*& material);

	void applyMaterial(const std::vector<Volume*> volumes, Material*& material);

	void generateMesh(const PartitionOfUnity& elementType, const MeshAlgorithm& algorithm = AUTO, std::string geofile = std::string(), 
					  const std::string& gmshPath = std::string(), const bool& plotMesh = true, const bool& showInfo = false);

	void solveStaticProblem();

	void solveTransientProblem();

	void solvePFEMProblem();

private:
	NodeMesh* getNode(const int& index);

	ElementMesh* getElement(const int& index);

	Material* getMaterial(const int& index);

	std::vector<DirichletBoundaryCondition*> getDirichletBoundaryConditions();

	std::vector<NeumannBoundaryCondition*> getNeumannBoundaryConditions();

	double getInitialPositionNorm() const;

	void setPastVariables();

	void computeCurrentVariables();

	void computeIntermediateVariables();

	void getConstrainedDOFs(int& ndofs, int*& constrainedDOFs);

	void getExternalForces(int& ndofs, int*& indexes, double*& externalForces);

	void assembleStaticLinearSystem(Mat& mat, Vec& vec);

	void assembleTransientLinearSystem(Mat& mat, Vec& vec);

	void assemblePFEMLinearSystem(Mat& mat, Vec& vec);

	void applyNeummanConditions(Vec& vec, int& ndofs, int*& indexes, double*& externalForces, const double& loadFactor);

	void solveLinearSystem(KSP& ksp, Mat& mat, Vec& rhs, Vec& solution);

	void updateVariables(Vec& solution, double& positionNorm, double& pressureNorm);

	void computeCauchyStress();

	void computeInitialAccel();

	void exportToParaview(const int& step);

	void readInput(const std::string& inputFile, const bool& deleteFiles, const PartitionOfUnity elementType);

	void transferGeometricBoundaryConditions();

	void domainDecompositionMETIS(const PartitionOfUnity& elementType);

	void remesh();

	void alphaShape();

	void identifyFreeSurfaces();

	void computeNodesCloudArea();

private:
	int numberOfSteps_;
	int maxNonlinearIterations_;
	double nonlinearTolerance_;
	double deltat_;
	double gravity_[3];
	double rhoInf_;
	double alpha_;
	double meshLength_;
	bool initialAccel_;

	int index_;
	int dimension_;
	int numberOfDOFs_;
	Geometry * geometry_;
	std::vector<NodeMesh*> nodes_;
	std::vector<NodeMesh*> isolatedParticles_;
	std::vector<ElementMesh*> elements_;
	std::vector<Material*> materials_;
	std::vector<DirichletBoundaryCondition*> dirichletBoundaryConditions_;
	std::vector<NeumannBoundaryCondition*> neumannBoundaryConditions_;
	idx_t* elementPartition_;
	idx_t* nodePartition_;
};

#endif



