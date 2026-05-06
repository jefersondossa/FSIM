//------------------------------------------------------------------------------
// 
//                   Jeferson W D Fernandes and Rodolfo A K Sanches
//                             University of Sao Paulo
//                           (C) 2017 All Rights Reserved
//
// <LicenseText>
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//------------------------------------FLUID-------------------------------------
//------------------------------------------------------------------------------

#ifndef FLUID_H
#define FLUID_H

#include "Element.h"
#include "Boundary.h"
#include "fluidDomain.h"

#include<cstdlib>
#include<fstream>
#include<iostream>
#include <petscviewerhdf5.h>

// PETSc libraries
#include <metis.h>
#include <petscksp.h> 




/// Mounts the incompressible flow problem

template<int DIM, int DEG>
class Fluid{
public:
    /// Defines the class Element locally
    typedef Element<DIM,DEG> Elements;

    /// Defines the class Node locally
    typedef typename Elements::Nodes  Node;

    /// Defines the class Boundary locally
    typedef Boundary<DIM,DEG> Boundaries;

    /// Defines the class Fluid Parameters locally
    typedef FluidParameters<DIM,DEG> Parameters;

    /// Defines the class Numerical integration locally
    typedef DomainIntegration<DIM,DEG> DIntegration;

    /// Defines the vector of fluid nodes
    std::vector<Node *>       nodes_;

    /// Defines the vector of fluid elements
    std::vector<Elements *>   elements_;
 
    /// Defines the vector of fluid boundaries mesh nodes
    std::vector<Boundaries *> boundary_;

public:
    //FLUID VARIABLES
    std::string inputFile; //Fluid input file
    int numElem;           //Number of elements in fluid mesh 
    int numNodes;          //Number of nodes in velocity/quadratic mesh
    int numBoundaries;     //Number of fluid boundaries
    int numBoundElems;     //Number of elements in fluid boundaries
    int numDOF;
    double pressInf;       //Undisturbed pressure 
    double rhoInf;         //Density
    double tempInf;        //Temperature
    double viscInf;        //Viscosity
    double ktermInf;       //Thermal condutivity
    VecDouble velocityInf; //Undisturbed velocity
    VecDouble fieldForces; //Field forces (constant)
    PetscInt* part_elem;      //Fluid Domain Decomposition - Elements
    PetscInt* part_nodes;     //Fluid Domain Decomposition - Nodes
    int numTimeSteps;      //Number of Time Steps
    int printFreq;         //Printing frequence of output files
    double dTime;          //Time Step
    int rank, size;
    int numFSIInterfaces;
    int iAux;
    bool computeDragAndLift;
    int iTimeStep;
    Mat               A;
    int nElNodes = (3+(DIM-2)*DEG)*(2+3*DEG+DEG*DEG)/6;
    int nLocDOF = -8*DIM -21*DEG + 15*DIM*DEG + 16;
    int nBdNodes = 3*(1-DEG)+DIM*(2*DEG-1);

public:
    std::vector<int> dragAndLiftBoundary;
    int numberOfLines;

    double glueZoneThickness;
    double arlequinEpsilon;
    int weightFunctionBehavior;
    bool printVelocity;
    bool printAcceleration;
    bool printRealVelocity;
    bool printLagrangeMultipliers;
    bool printElementCorrespondence;
    bool printDistFunction;
    bool printEnergyWeightFunction;
    bool printPressure;
    bool printRealPressure;
    bool printVorticity;
    bool printInnerNormal;
    bool printMeshVelocity;
    bool printMeshDisplacement;
    bool printGlueZone;
    bool printJacobian;
    bool printProcess;
    bool printLines;
    double integScheme;    //Time Integration Scheme
    Parameters fluidParameters;
    DIntegration* numIntegration; //Numerical integration

    ProblemType fProbType = ProblemType::ENavierStokes;

public:
    void setProblemType(ProblemType ptype){fProbType = ptype;};
    ProblemType &getProblemType(){return fProbType;}

    void meshReading(Geometry* &geometry_, const std::string& inputFile, const std::string& inputMesh, const std::string& mirror, const bool& deleteFiles);

    void readInitialValues(const std::string& inputPrev, const std::string& inputCurr);
    void readInputFile(const std::string& inputFile, std::ofstream& mirrorData);
    
    /// Reads the mesh nodes from a .msh file
    /// @param std::string input .msh file @param std::string mirror file
    void readNodes(std::ifstream &file, std::ofstream& mirrorData);

    /// Reads the mesh elements from a .msh file
    /// @param Geometry* mesh geometry
    /// @param std::string input .msh file @param std::string mirror file
    /// @param std::vector<Elements*> auxiliary vector of Elements
    /// @param std::unordered_map<int, std::string> mesh physical entities
    void readElements(Geometry* &geometry_,std::ifstream &file, std::ofstream& mirrorData, std::vector<Elements*> &elementsAux_, std::unordered_map<int, std::string> &physicalEntities);
    void renumberConnectivity();
    void setBoundaryConstrains();
    void setBoundarySides();

    Parameters &getFluidParameters(){
        return fluidParameters;
    }

    /// Performs the domain decomposition for parallel processing
    void domainDecompositionMETIS(); 

    int getNumberOfElements(){return numElem;}
    int getNumberOfNodes(){return numNodes;}

    /// Export the domain decomposition 
    /// @return pair with the elements and nodes domain decompositions
    std::pair<idx_t*,idx_t*> getDomainDecomposition(){
        return std::make_pair(part_elem,part_nodes);};

    /// Gets the number of time steps
    /// @return number of time steps
    int &getNumberOfTimeSteps(){return numTimeSteps;};

    /// Gets the time step size
    /// @return time step size
    double getTimeStep(){return dTime;};

    /// Gets the number of fluid-structure interfaces
    /// @return number of fluid boundaries which composes the 
    /// fluid structure interface
    int getNumberofFSIInterfaces(){return numFSIInterfaces;};
    
    /// Gets the flag for computing the Drag and Lift coefficients in a 
    /// specific boundary
    /// @return bool flag for computing Drag and Lift coefficients
    bool getComputeDragAndLift() {return computeDragAndLift;};

    /// Gets the number of the boundary for computing the Drag and Lift 
    /// coefficients
    /// @return int boundary number
    int getDragAndLiftBoundary(int index) {return dragAndLiftBoundary[index];};
 
    /// Mounts and solve the transient incompressible flow problem    
    /// @param int maximum number of Newton-Raphson's iterations
    /// @param double tolerance of the Newton-Raphson's process
    /// @param int problem type: 1 - Stokes problem; 2 - Navier-Stokes problem.
    int solveTransientProblem(int iterNumber,double tolerance, bool accel);

    /// Mounts and solve the transient incompressible flow problem for moving
    /// domain problems
    /// @param int maximum number of Newton-Raphson's iterations
    /// @param double tolerance of the Newton-Raphson's process
    /// @param int problem type: 1 - Stokes problem; 2 - Navier-Stokes problem.
    int solveTransientProblemMoving(int iterNumber,double tolerance, bool accel);

    /// Mounts and solve the steady Laplace problem
    /// @param int maximum number of Newton-Raphson's iterations
    /// @param double tolerance of the Newton-Raphson's process
    int solveSteadyLaplaceProblem(int iterNumber, double tolerance);

    /// Mounts and solve the transient incompressible flow problem for 
    /// for fluid structure interaction problems
    /// @param int maximum number of Newton-Raphson's iterations
    /// @param double tolerance of the Newton-Raphson's process
    /// @param int problem type: 1 - Stokes problem; 2 - Navier-Stokes problem.
    int solveFSIFluid(int iterNumber,
                      double tolerance,
                      int problem_type);
    int solvePoisson();

    /// Print the results for Paraview post-processing
    /// @param int time step
    // void printResults(int step);
    void printResultsPoisson();

    /// Compute and print drag and lift coefficients
    void dragAndLiftCoefficients(std::ofstream& dragLift);

    /// Gets the fluid model nodes and export for solving the overlapping
    /// mesh problem with the Arlequin method
    /// @return fluid model nodes information
    std::vector<Node *> getNodesVelocity(){return nodes_;}

    /// Gets the fluid model elements and export for solving the overlapping
    /// mesh problem with the Arlequin method
    /// @return fluid model elements information
    std::vector<Elements *> getElements(){return elements_;}

    std::vector<std::string> split2(std::string str, std::string delim)
    {
        std::istringstream is(str);
        std::vector<std::string> values;
        std::string token;
        while (getline(is, token, ' '))
            values.push_back(token);
        return values;
    }

    void computeError(VecDouble &errorsTotal);
};


#endif
