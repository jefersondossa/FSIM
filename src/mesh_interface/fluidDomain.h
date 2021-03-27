#pragma once

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








std::vector<std::string> split2(std::string str, std::string delim)
{
	std::istringstream is(str);
	std::vector<std::string> values;
	std::string token;
	while (getline(is, token, ' '))
		values.push_back(token);
	return values;
}

// Public methods
FluidDomain::FluidDomain(Geometry* geometry, const int& index)
	: geometry_(geometry), index_(index), numberOfSteps_(1), maxNonlinearIterations_(3), nonlinearTolerance_(1.0e-6), deltat_(1.0), gravity_{0.0, -9.81, 0.0}, 
	  rhoInf_(1.0), alpha_(1.5), meshLength_(1.0), initialAccel_(false)
{
	system("mkdir -p ./results");
	system("rm ./results/*.vtu 2> /dev/null");
}

FluidDomain::~FluidDomain() {}

void FluidDomain::setNumberOfSteps(const int numberOfSteps)
{
	numberOfSteps_ = numberOfSteps;
}

void FluidDomain::setMaxNonlinearIterations(const int maxNonlinearIterations)
{
	maxNonlinearIterations_ = maxNonlinearIterations;
}

void FluidDomain::setNonlinearTolerance(const double nonlinearTolerance)
{
	nonlinearTolerance_ = nonlinearTolerance;
}

void FluidDomain::setDeltat(const double deltat)
{
	deltat_ = deltat;
}

void FluidDomain::setGravity(const double x, const double y, const double z)
{
	gravity_[0] = x; gravity_[1] = y; gravity_[2] = z;
}

void FluidDomain::setSpectralRadius(const double rhoInf)
{
	rhoInf_ = rhoInf;
}

void FluidDomain::setAlpha(const double& alpha)
{
	alpha_ = alpha;
}

void FluidDomain::setMeshLength(const double& h)
{
	meshLength_ = h;
}

void FluidDomain::setInitialAcceleration(const bool initialAccel)
{
	initialAccel_ = initialAccel;
}

// void FluidDomain::setReferenceConfiguration(const ReferenceConfiguration reference)
// {
// 	elements_[0]->setReferenceConfiguration(reference);
// }

void FluidDomain::applyMaterial(const std::vector<Line*> lines, Material*& material)
{
	materials_.push_back(material);
	for (Line* line : lines)
		line->setMaterial(material);
}

void FluidDomain::applyMaterial(const std::vector<Surface*> surfaces, Material*& material)
{
	materials_.push_back(material);
	for (Surface* surface : surfaces)
		surface->setMaterial(material);
}

void FluidDomain::applyMaterial(const std::vector<Volume*> volumes, Material*& material)
{
	materials_.push_back(material);
	for (Volume* volume : volumes)
		volume->setMaterial(material);
}

void FluidDomain::generateMesh(const PartitionOfUnity& elementType, const MeshAlgorithm& algorithm, std::string geofile, const std::string& gmshPath,
							   const bool& plotMesh, const bool& showInfo)
{	
	int rank, size;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	MPI_Comm_size(PETSC_COMM_WORLD, &size);

	std::pair<std::string, bool> pair; pair.second = false;
	if (rank == 0)
	{
		pair = createMesh(geometry_, elementType, algorithm, geofile, gmshPath, plotMesh, showInfo);

		for (int i = 1; i < size; i++)
		{
			MPI_Send(pair.first.c_str(), pair.first.length()+1, MPI_CHAR, i, 0, PETSC_COMM_WORLD);
			if (i == size-1)
			{
				MPI_Send(&pair.second, 1, MPI_C_BOOL, i, 0, PETSC_COMM_WORLD);
				pair.second = false;
			}
		}
	}
	else
	{
		MPI_Status status;
		MPI_Probe(0, 0, PETSC_COMM_WORLD, &status);
		int count;
		MPI_Get_count(&status, MPI_CHAR, &count);
		char buf[count+1];
		MPI_Recv(&buf, count+1, MPI_CHAR, 0, 0, PETSC_COMM_WORLD, &status);
		pair.first = buf;
		if (rank == size-1)
			MPI_Recv(&pair.second, 1, MPI_C_BOOL, 0, 0, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	// readInput(pair.first, pair.second, elementType);
	// transferGeometricBoundaryConditions();
	// domainDecompositionMETIS(elementType);
}

// Private methods
NodeMesh* FluidDomain::getNode(const int& index)
{
	return nodes_[index];
}

ElementMesh* FluidDomain::getElement(const int& index)
{
	return elements_[index];
}

Material* FluidDomain::getMaterial(const int& index)
{
	return materials_[index];
}

std::vector<DirichletBoundaryCondition*> FluidDomain::getDirichletBoundaryConditions()
{
	return dirichletBoundaryConditions_;
}

std::vector<NeumannBoundaryCondition*> FluidDomain::getNeumannBoundaryConditions()
{
	return neumannBoundaryConditions_;
}

void FluidDomain::readInput(const std::string& inputFile, const bool& deleteFiles, const PartitionOfUnity elementType)
{
	int rank;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

	//defyning the maps that are used to store the elements information
	const std::unordered_map<int, PartitionOfUnity> gmshElement = 
	{
		//vertex
		{15, VERTEX},
		//lines
		{1, L2}, {8, L3}, {26, L4},
		//triangles
		{2, T3}, {9, T6}, {21, T10},
		//quadrilaterals
		{3, Q4}, {10, Q9}, {36, Q16},
		//tetrahedrons
		{4, TET4}, {11, TET10}, {29, TET20}
	};
	const std::unordered_map<PartitionOfUnity, int> numNodes = 
	{
		//vertex
		{VERTEX, 1},
		//lines
		{L2, 2}, {L3, 3}, {L4, 4},
		//triangles
		{T3, 3}, {T6, 6}, {T10, 10},
		//quadrilaterals
		{Q4, 4}, {Q9, 9}, {Q16, 16},
		//tetrahedrons
		{TET4, 4}, {TET10, 10}, {TET20, 20}
	};
	const std::unordered_map<PartitionOfUnity, int> dimension =
	{
		//lines
		{L2, 3}, {L3, 3}, {L4, 3},
		//triangles
		{T3, 2}, {T6, 2}, {T10, 2},
		//quadrilaterals
		{Q4, 2}, {Q9, 2}, {Q16, 2},
		//tetrahedrons
		{TET4, 3}, {TET10, 3}, {TET20, 3}
	};

	dimension_ = dimension.at(elementType);

	bool mixed = false;
	if (materials_[0]->getType() == MaterialType::ELASTIC_INCOMPRESSIBLE_SOLID ||
		materials_[0]->getType() == MaterialType::NEWTONIAN_INCOMPRESSIBLE_FLUID ) mixed = true;

	//opening the .msh file
	std::ifstream file(inputFile);
	std::string line;
	std::getline(file, line); std::getline(file, line); std::getline(file, line); std::getline(file, line);

	//reading physical entities
	int nEntities;
	file >> nEntities;
	std::getline(file, line);
	std::unordered_map<int, std::string> physicalEntities;
	physicalEntities.reserve(nEntities);
	for (int i = 0; i < nEntities; i++)
	{
		std::getline(file, line);
		std::vector<std::string> tokens = split2(line, " ");
		int index;
		std::istringstream(tokens[1]) >> index;
		physicalEntities[index] = tokens[2].substr(1, tokens[2].size() - 2);
	}
	std::getline(file, line); std::getline(file, line);

	//reading nodes
	int nNodes;
	file >> nNodes;
	nodes_.reserve(nNodes);
	std::getline(file, line);
	for (int i = 0; i < nNodes; i++)
	{
		std::getline(file, line);
		std::vector<std::string> tokens = split2(line, " ");
		std::vector<DegreeOfFreedom*> degreesOfFreedom;
        mixed? degreesOfFreedom.reserve(dimension_+1) : degreesOfFreedom.reserve(dimension_);
        double coord_[3] = {};
		for (int j = 0; j < dimension_; j++)
		{
			double coord;
			std::istringstream(tokens[j+1]) >> coord;
			coord_[j] = coord;
            degreesOfFreedom.emplace_back(new DegreeOfFreedom(DOFType::POSITION, coord));
		}
		nodes_.push_back(new NodeMesh(i, coord_));
	}
	std::getline(file, line); std::getline(file, line);

	//reading elements
	int nElements;
	file >> nElements;
	std::getline(file, line);
	int contLine = 0;
	int contSurface = 0;
	int contVolume = 0;
	elements_.reserve(nElements);
	for (int i = 0; i < nElements; i++)
	{
		std::getline(file, line);
		std::vector<std::string> tokens = split2(line, " ");
		std::vector<int> values(tokens.size(), 0);
		for (size_t j = 0; j < tokens.size(); j++)
			std::istringstream(tokens[j]) >> values[j];
		PartitionOfUnity type = gmshElement.at(values[1]);
		int numberOfNodes = numNodes.at(type);
		std::vector<NodeMesh*> elementNodes;
		elementNodes.reserve(numberOfNodes);
		for (size_t j = 5 ; j < values.size(); j++)
			elementNodes.push_back(nodes_[values[j]-1]);
		std::string name = physicalEntities[values[3]];
		
		//Adding 3D elements to volumes
		if (name[0] == 'v')
		{
			// Volume* object = geometry_->getVolume(name);
			// ParametricVolumeElement* elem = new ParametricVolumeElement(contVolume, elementNodes, type);
			// object->addParametricElement(elem);
			// object->addNodes(elementNodes);
			// switch (elementType)
			// {
			// 	case TET4:
			// 	case TET10:
			// 	case TET20:
			// 		elem->setPlot(true);
			// 		Material* mat = object->getMaterial();
					// ElementMesh* element = new VolumeElement(contVolume, elem, mat);
					// elements_.push_back(element);
					// object->addElement(element);
					// break;
			// }
			// contVolume++;
		}
		//Adding 2D elements to surfaces
		else if (name[0] == 's')
		{	
			// Surface* object = geometry_->getSurface(name);
			// ParametricSurfaceElement* elem = new ParametricSurfaceElement(contSurface, elementNodes, type);
			// object->addParametricElement(elem);
			// object->addNodes(elementNodes);
			// switch (elementType)
			// {
			// 	case T3:
			// 	case T6:
			// 	case T10:
			// 	case Q4:
			// 	case Q9:
			// 	case Q16:
			// 		elem->setPlot(true);
			// 		Material* mat = object->getMaterial();
			// 		// ElementMesh* element = new PlaneElement(contSurface, elem, mat);
			// 		// elements_.push_back(element);
			// 		// object->addElement(element);
			// 		break;
			// }
			// contSurface++;
		}
		//Adding 1D elements to lines
		else if (name[0] == 'l')
		{
			// Line* object = geometry_->getLine(name);
			// ParametricLineElement* elem = new ParametricLineElement(contLine, elementNodes, type);
			// object->addParametricElement(elem);
			// object->addNodes(elementNodes);
			// switch (elementType)
			// {
			// 	case L2:
			// 	case L3:
			// 	case L4:
			// 		elem->setPlot(true);
			// 		Material* mat = object->getMaterial();
			// 		// ElementMesh* element = new LineElement(contLine, elem, mat);
			// 		// elements_.push_back(element);
			// 		// object->addElement(element);
			// 		break;
			// }
			// contLine++;
		}
		//Adding a node to point
		else
		{
			Point* object = geometry_->getPoint(name);
			object->addNode(getNode(elementNodes[0]->getIndex()));
		}
	}

	//Closing the file
	file.close();
	if (deleteFiles)
		system((rm + inputFile).c_str());
}

