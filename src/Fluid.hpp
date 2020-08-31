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

#include "Element.hpp"
#include "Boundary.hpp"
#include "fluidDomain.h"

#include<cstdlib>
#include<fstream>
#include<iostream>


// PETSc libraries
#include <metis.h>
#include <petscksp.h> 

/// Mounts the incompressible flow problem

template<int DIM>
class Fluid{
public:
    /// Defines the class Element locally
    typedef Element<DIM> Elements;

    /// Defines the class Node locally
    typedef typename Elements::Nodes  Node;

    /// Defines the class Boundary locally
    typedef Boundary<DIM> Boundaries;

    typedef FluidParameters<DIM> Parameters;

    /// Defines the vector of fluid nodes
    std::vector<Node *>       nodes_;

    /// Defines the vector of fluid elements
    std::vector<Elements *>   elements_;
 
    /// Defines the vector of fluid boundaries mesh nodes
    std::vector<Boundaries *> boundary_;

private:
    //FLUID VARIABLES
    std::string inputFile; //Fluid input file
    int numElem;           //Number of elements in fluid mesh 
    int numNodes;          //Number of nodes in velocity/quadratic mesh
    int numBoundaries;     //Number of fluid boundaries
    int numBoundElems;     //Number of elements in fluid boundaries
    double pressInf;       //Undisturbed pressure 
    double rhoInf;         //Density
    double tempInf;        //Temperature
    double viscInf;        //Viscosity
    double ktermInf;       //Thermal condutivity
    double velocityInf[3]; //Undisturbed velocity
    double fieldForces[3]; //Field forces (constant)
    idx_t* part_elem;      //Fluid Domain Decomposition - Elements
    idx_t* part_nodes;     //Fluid Domain Decomposition - Nodes
    int numTimeSteps;      //Number of Time Steps
    int printFreq;         //Printing frequence of output files
    double dTime;          //Time Step
    int rank;
    int numFSIInterfaces;
    int iAux;
    bool computeDragAndLift;
    int iTimeStep;
    Mat               A;
    
    


    
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

public:


    void meshReading(Geometry* &geometry_, const std::string& inputFile, const std::string& inputMesh, const std::string& mirror, const bool& deleteFiles);

    void readInitialValues(const std::string& inputPrev, const std::string& inputCurr);


    /// Performs the domain decomposition for parallel processing
    void domainDecompositionMETIS(); 

    /// Export the domain decomposition 
    /// @return pair with the elements and nodes domain decompositions
    std::pair<idx_t*,idx_t*> getDomainDecomposition(){
        return std::make_pair(part_elem,part_nodes);};

    /// Gets the number of time steps
    /// @return number of time steps
    int getNumberOfTimeSteps(){return numTimeSteps;};

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

    int getInitialAcceleration();

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

    /// Print the results for Paraview post-processing
    /// @param int time step
    void printResults(int step);

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
};

//------------------------------------------------------------------------------
//--------------------------------IMPLEMENTATION--------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//---------------------SUBDIVIDES THE FINITE ELEMENT DOMAIN---------------------
//------------------------------------------------------------------------------
template<>
void Fluid<2>::domainDecompositionMETIS() {
    
    std::string mirror2;
    mirror2 = "domain_decomposition.txt";
    std::ofstream mirrorData(mirror2.c_str());
    
    int size;

    MPI_Comm_size(PETSC_COMM_WORLD, &size);

    idx_t objval;
    idx_t numEl = numElem;
    idx_t numNd = numNodes;
    idx_t dd = 2;
    idx_t ssize = size;
    idx_t one = 1;
    idx_t elem_start[numEl+1], elem_connec[(4*dd-2)*numEl];
    part_elem = new idx_t[numEl];
    part_nodes = new idx_t[numNd];


    for (idx_t i = 0; i < numEl+1; i++){
        elem_start[i]=(4*dd-2)*i;
    };
    for (idx_t jel = 0; jel < numEl; jel++){
        int *connec=elements_[jel]->getConnectivity();        
        
        for (idx_t i=0; i<(4*dd-2); i++){
        elem_connec[(4*dd-2)*jel+i] = connec[i];
        };
    };

    //Performs the domain decomposition
    METIS_PartMeshDual(&numEl, &numNd, elem_start, elem_connec, \
                              NULL, NULL, &one, &ssize, NULL, NULL,    \
                              &objval, part_elem, part_nodes);

    mirrorData << std::endl \
               << "FLUID MESH DOMAIN DECOMPOSITION - ELEMENTS" << std::endl;
    for(int i = 0; i < numElem; i++){
        mirrorData << "process = " << part_elem[i] \
                   << ", element = " << i << std::endl;
    };

    mirrorData << std::endl \
               << "FLUID MESH DOMAIN DECOMPOSITION - NODES" << std::endl;
    for(int i = 0; i < numNodes; i++){
        mirrorData << "process = " << part_nodes[i] \
                   << ", node = " << i << std::endl;
    };
    
    return;

};

//------------------------------------------------------------------------------
//----------------------------PRINT VELOCITY RESULTS----------------------------
//------------------------------------------------------------------------------
template<>
void Fluid<2>::printResults(int step) {

    //    std::cout << "Printing Velocity Results" << std::endl;
    std::string result;
    std::ostringstream convert;

    convert << step+100000;
    result = convert.str();
    std::string s = "saidaVel"+result+".vtu";
    
    std::fstream output_v(s.c_str(), std::ios_base::out);

    output_v << "<?xml version=\"1.0\"?>" << std::endl
             << "<VTKFile type=\"UnstructuredGrid\">" << std::endl
             << "  <UnstructuredGrid>" << std::endl
             << "  <Piece NumberOfPoints=\"" << numNodes
             << "\"  NumberOfCells=\"" << numElem
             << "\">" << std::endl;

    //WRITE NODAL COORDINATES
    output_v << "    <Points>" << std::endl
             << "      <DataArray type=\"Float64\" "
             << "NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

    for (int i=0; i<numNodes; i++){
        double* x = nodes_[i] -> getCoordinates();
        output_v << x[0] << " " << x[1] << " " << 0.0 << std::endl;        
    };
    output_v << "      </DataArray>" << std::endl
             << "    </Points>" << std::endl;
    
    //WRITE ELEMENT CONNECTIVITY
    output_v << "    <Cells>" << std::endl
             << "      <DataArray type=\"Int32\" "
             << "Name=\"connectivity\" format=\"ascii\">" << std::endl;
    
    for (int i=0; i<numElem; i++){
        int *connec = elements_[i] -> getConnectivity();
        output_v << connec[0] << " " << connec[1] << " " << connec[2] << " " \
                 << connec[3] << " " << connec[4] << " " << connec[5] << \
            std::endl;
    };
    output_v << "      </DataArray>" << std::endl;
  
    //WRITE OFFSETS IN DATA ARRAY
    output_v << "      <DataArray type=\"Int32\""
             << " Name=\"offsets\" format=\"ascii\">" << std::endl;
    
    int aux = 0;
    for (int i=0; i<numElem; i++){
        output_v << aux + 6 << std::endl;
        aux += 6;
    };
    output_v << "      </DataArray>" << std::endl;
  
    //WRITE ELEMENT TYPES
    output_v << "      <DataArray type=\"UInt8\" Name=\"types\" "
             << "format=\"ascii\">" << std::endl;
    
    for (int i=0; i<numElem; i++){
        output_v << 22 << std::endl;
    };

    output_v << "      </DataArray>" << std::endl
             << "    </Cells>" << std::endl;

    //WRITE NODAL RESULTS
    output_v << "    <PointData>" << std::endl;

    if (printVelocity){
        output_v<< "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                << "Name=\"Velocity\" format=\"ascii\">" << std::endl;
        for (int i=0; i<numNodes; i++){
            output_v << nodes_[i] -> getVelocity(0) << " "             
                     << nodes_[i] -> getVelocity(1) << " " << 0. << std::endl;
        };
        output_v << "      </DataArray> " << std::endl;
    };

    if (printMeshVelocity){
        output_v<< "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                << "Name=\"Mesh Velocity\" format=\"ascii\">" << std::endl;
        
        for (int i=0; i<numNodes; i++){
            output_v << nodes_[i] -> getMeshVelocity(0) << " "    
                     << nodes_[i] -> getMeshVelocity(1) << " " 
                     << 0. << std::endl;
        };
        output_v << "      </DataArray> " << std::endl;
    };

    if (printVorticity){
        output_v <<"      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                 << "Name=\"Vorticity\" format=\"ascii\">" << std::endl;
        for (int i=0; i<numNodes; i++){
            output_v << nodes_[i] -> getVorticity() << std::endl;
        };
        output_v << "      </DataArray> " << std::endl;
    };

    if (printPressure){
        output_v <<"      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                 << "Name=\"Pressure\" format=\"ascii\">" << std::endl;
        for (int i=0; i<numNodes; i++){
            output_v << 0. << " " << 0. << " " 
                     << nodes_[i] -> getPressure() << std::endl;
        };
        output_v << "      </DataArray> " << std::endl;
    };


    if (printMeshDisplacement){
        output_v <<"      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                 << "Name=\"Mesh Displacement\" format=\"ascii\">" << std::endl;
        for (int i=0; i<numNodes; i++){       
            double* x = nodes_[i] -> getCoordinates();
            double* xp = nodes_[i] -> getInitialCoordinates();
            
            output_v << x[0]-xp[0] << " " << x[1]-xp[1] << " " 
                     << 0. << std::endl;
        };
        output_v << "      </DataArray> " << std::endl;
    };

    output_v << "    </PointData>" << std::endl; 

    //WRITE ELEMENT RESULTS
    output_v << "    <CellData>" << std::endl;
    
    if (printProcess){
        output_v <<"      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                 << "Name=\"Process\" format=\"ascii\">" << std::endl;
        for (int i=0; i<numElem; i++){
            output_v << part_elem[i] << std::endl;
        };
        output_v << "      </DataArray> " << std::endl;
    };

    output_v <<"      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
             << "Name=\"Lines\" format=\"ascii\">" << std::endl;
    for (int i=0; i<numElem; i++){
        int res = 0;
        for (int j=0; j<numBoundElems; j++){
           if (boundary_[j] -> getElement() == i) res = boundary_[j] -> getBoundaryGroup();
        }
        output_v << res << std::endl;
    };
    output_v << "      </DataArray> " << std::endl;

    if (printJacobian){
        output_v <<"      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                 << "Name=\"Jacobian\" format=\"ascii\">" << std::endl;
        for (int i=0; i<numElem; i++){
            output_v << elements_[i] -> getJacobian() << std::endl;
        };
        output_v << "      </DataArray> " << std::endl;
    };
    
    output_v << "    </CellData>" << std::endl; 

    //FINALIZE OUTPUT FILE
    output_v << "  </Piece>" << std::endl;
    
    // output_v << "  <FieldData>" << std::endl;

    // output_v << "      <DataArray type=\"Float64\" Name=\"Time\" NumberOfTuples=\"1\" "
    //          << " format=\"ascii\">" << std::endl;
    // output_v << step << std::endl;
    // output_v << "      </DataArray> " << std::endl;

    // output_v << "      <DataArray type=\"Float64\" Name=\"LiftCoefficient\" NumberOfTuples=\"1\" "
    //          << " format=\"ascii\">" << std::endl;
    // output_v << liftCoefficient << std::endl;
    // output_v << "      </DataArray> " << std::endl;

    // // output_v << "      <DataSet type=\"Float64\" Name=\"Drag Coefficient\" NumberOfTuples=\"1\" "
    // //          << " format=\"ascii\">" << std::endl;
    // // output_v << dragCoefficient << std::endl;
    // // output_v << "      </DataSet> " << std::endl;

    // output_v << "  </FieldData>" << std::endl
    output_v << "  </UnstructuredGrid>" << std::endl
             << "</VTKFile>" << std::endl;

};


//------------------------------------------------------------------------------
//----------------------COMPUTES DRAG AND LIFT COEFFICIENTS---------------------
//------------------------------------------------------------------------------
template<>
void Fluid<2>::dragAndLiftCoefficients(std::ofstream& dragLift){

    double dragCoefficient = 0.;
    double liftCoefficient = 0.;
    double pressureDragCoefficient = 0.;
    double pressureLiftCoefficient = 0.;
    double frictionDragCoefficient = 0.;
    double frictionLiftCoefficient = 0.;
    
    for (int jel = 0; jel < numBoundElems; jel++){   
        
        double rhoInf = 1.0;
        double velocityInf[2];
        velocityInf[0] = -1.;
        velocityInf[1] = 0.;
        
        double dForce = 0.;
        double lForce = 0.;
        double pDForce = 0.;
        double pLForce = 0.;
        double fDForce = 0.;
        double fLForce = 0.;
        double aux_Mom = 0.;
        double aux_Per = 0.;
        
       for (int i=0; i<numberOfLines; i++){
            //std::cout << "Bound group " << boundary_[jel] -> getBoundaryGroup() << std::endl;
            if (boundary_[jel] -> getBoundaryGroup() == dragAndLiftBoundary[i]){
                //std::cout << "AQUI " << numberOfLines<< " " << i << " " << dragAndLiftBoundary[i] << std::endl;
                int iel = boundary_[jel] -> getElement();
                elements_[iel] -> computeDragAndLiftForces(pDForce, pLForce, fDForce, fLForce, dForce, lForce, aux_Mom, aux_Per);
                elements_[iel] -> computeSeparationAngle();
            };
        };
        
        pressureDragCoefficient += pDForce ;/// 
            //(0.5 * rhoInf * velocityInf[0] * velocityInf[0]);
        pressureLiftCoefficient += pLForce;// / 
            //(0.5 * rhoInf * velocityInf[0] * velocityInf[0]);
        
        frictionDragCoefficient += fDForce / 
            (0.5 * rhoInf * velocityInf[0] * velocityInf[0]);
        frictionLiftCoefficient += fLForce / 
            (0.5 * rhoInf * velocityInf[0] * velocityInf[0]);
        
        dragCoefficient += dForce / 
            (0.5 * rhoInf * velocityInf[0] * velocityInf[0]);
        liftCoefficient += lForce / 
            (0.5 * rhoInf * velocityInf[0] * velocityInf[0]);
        
    };
    // std::cout << "vazao " << pressureDragCoefficient << " " << pressureLiftCoefficient << std::endl;
    if (rank == 0) {
        const int timeWidth = 11;
        const int numWidth = 11;
        dragLift << std::setprecision(3) << std::scientific;
        dragLift << std::left << std::setw(timeWidth) << iTimeStep * dTime;
        dragLift << std::setw(numWidth) << pressureDragCoefficient;
        dragLift << std::setw(numWidth) << pressureLiftCoefficient;
        dragLift << std::setw(numWidth) << frictionDragCoefficient;
        dragLift << std::setw(numWidth) << frictionLiftCoefficient;
        dragLift << std::setw(numWidth) << dragCoefficient;
        dragLift << std::setw(numWidth) << liftCoefficient;
        dragLift << std::endl;
    }
}


template<>
void Fluid<2>::meshReading(Geometry* &geometry_, const std::string& inputFile, const std::string& inputMesh, const std::string& mirror, const bool& deleteFiles) {

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //+++++++++++++++++++++++++++++OPPENING FILES+++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //defyning the maps that are used to store the elements information
    std::unordered_map<int, std::string> gmshElement = { {1, "line"}, {2, "triangle"}, {3, "quadrilateral"}, {8, "line3"}, {9, "triangle6"}, {10, "quadrilateral9"}, {15, "vertex"}, {16, "quadrilateral8"}, {20, "triangle9"}, {21, "triangle10"}, {26, "line4"}, {36, "quadrilateral16"}, {39, "quadrilateral12"} };
    std::unordered_map<std::string, int> numNodes2 = { {"vertex", 1}, {"line", 2}, {"triangle", 3}, {"quadrilateral", 4}, {"line3", 3}, {"triangle6", 6}, {"quadrilateral8", 8}, {"quadrilateral9", 9}, {"line4", 4}, {"triangle", 9}, {"triangle10", 10}, {"quadrilateral12", 12}, {"quadrilateral16", 16}};
    std::unordered_map<std::string, std::string> supportedElements = { {"triangle", "T3"}, {"triangle6", "T6"}, {"triangle10", "T10"}, {"quadrilateral", "Q4"}, {"quadrilateral8", "Q8"}, {"quadrilateral9", "Q9"}, {"quadrilateral12", "Q12"}, {"quadrilateral16", "Q16"} };
    std::unordered_map<Line*, std::vector< std::vector<int> >> lineElements;

    //opening the .msh file
    std::ifstream inputData(inputFile.c_str());
    std::ofstream mirrorData(mirror.c_str());
    std::ifstream file(inputMesh);
    std::string line;
    std::getline(file, line); std::getline(file, line); std::getline(file, line); std::getline(file, line);
  

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++READING PROBLEM VARIABLES+++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    
    //Read number of nodes, elements, time steps and printing frequence
    inputData >> numTimeSteps >> printFreq;
    mirrorData << "Number of Time Steps   = " << numTimeSteps << std::endl;
    mirrorData << "Printing Frequence     = " << printFreq << std::endl;
    
    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    
    //Read undisturbed velocity and pressure components
    inputData >> velocityInf[0] >> velocityInf[1] >> velocityInf[2] >> pressInf;

    mirrorData << "Undisturbed Velocity x = " << velocityInf[0] << std::endl;
    mirrorData << "Undisturbed Velocity y = " << velocityInf[1] << std::endl;
    mirrorData << "Undisturbed Velocity z = " << velocityInf[2] << std::endl;
    mirrorData << "Undisturbed Pressure   = " << pressInf << std::endl;

    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);getline(inputData,line);

    //Read undisturbed density, temperature, viscosity and thermal condutivity
    inputData >> rhoInf >> viscInf;

    mirrorData << "Undisturbed Density    = " << rhoInf << std::endl;
    mirrorData << "Undisturbed Viscosity  = " << viscInf << std::endl;

    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);getline(inputData,line);

    //Read time step lenght
    inputData >> dTime >> integScheme;

    mirrorData << "Time Step              = " << dTime << std::endl;
    mirrorData << "Time Integration Scheme= " << integScheme << std::endl;

    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);getline(inputData,line);getline(inputData,line);

    //Read field forces
    inputData >> fieldForces[0] >> fieldForces[1] >> fieldForces[2];

    mirrorData << "Field Forces x         = " << fieldForces[0] << std::endl;
    mirrorData << "Field Forces y         = " << fieldForces[1] << std::endl;
    mirrorData << "Field Forces z         = " << fieldForces[2] << std::endl \
               << std::endl;

    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);getline(inputData,line);getline(inputData,line);




    //Read Arlequin variables
    double k1,k2;
    inputData >> glueZoneThickness >> arlequinEpsilon >> weightFunctionBehavior
              >> k1 >> k2;

    mirrorData << "Glue Zone Thickness    = " << glueZoneThickness << std::endl;
    mirrorData << "Epsilon                = " << arlequinEpsilon << std::endl;
    mirrorData << "Energy Weight Function = " << weightFunctionBehavior
               << std::endl;
    mirrorData << "K1                     = " << k1 << std::endl;
    mirrorData << "K2                     = " << k2 << std::endl << std::endl;

    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);getline(inputData,line);

    fluidParameters.setViscosity(viscInf);
    fluidParameters.setDensity(rhoInf);
    fluidParameters.setTimeStep(dTime);
    fluidParameters.setSpectralRadius(integScheme);
    fluidParameters.setFieldForce(fieldForces);
    fluidParameters.setArlequinOperatorConstants(k1,k2);

    //Drag and lift
    inputData >> computeDragAndLift >> numberOfLines; 
    dragAndLiftBoundary.reserve(numberOfLines);
    for (int i = 0; i < numberOfLines; ++i)
    {
        int aux;
        inputData >> aux; 
        dragAndLiftBoundary.push_back(aux);
    }
    

    mirrorData << "Compute Drag and Lift  = " << computeDragAndLift<< std::endl;
    mirrorData << "Number of Lines  = " << numberOfLines << std::endl;
    for (int i = 0; i < numberOfLines; ++i)
    {
        mirrorData << "Lines  = " << dragAndLiftBoundary[i] << std::endl;
    }

    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);getline(inputData,line);

    //Printing results
    inputData >> printVelocity;              getline(inputData,line);
    inputData >> printAcceleration;              getline(inputData,line);
    inputData >> printRealVelocity;          getline(inputData,line);
    inputData >> printLagrangeMultipliers;   getline(inputData,line);
    inputData >> printElementCorrespondence; getline(inputData,line);
    inputData >> printDistFunction;          getline(inputData,line);
    inputData >> printEnergyWeightFunction;  getline(inputData,line);
    inputData >> printPressure;              getline(inputData,line);
    inputData >> printRealPressure;          getline(inputData,line);
    inputData >> printVorticity;             getline(inputData,line);
    inputData >> printInnerNormal;           getline(inputData,line);
    inputData >> printMeshVelocity;          getline(inputData,line);
    inputData >> printMeshDisplacement;      getline(inputData,line);
    inputData >> printGlueZone;              getline(inputData,line);
    inputData >> printJacobian;              getline(inputData,line);
    inputData >> printProcess;               getline(inputData,line);
    inputData >> printLines;              

    mirrorData << "PrintVelocity              = " << printVelocity << std::endl;
    mirrorData << "printAcceleration          = " << printAcceleration << std::endl;
    mirrorData << "PrintRealVelocity          = " << printRealVelocity << std::endl;
    mirrorData << "PrintLagrangeMultipliers   = " << printLagrangeMultipliers << std::endl;
    mirrorData << "PrintElementCorrespondence = " << printElementCorrespondence << std::endl;
    mirrorData << "PrintDistFunction          = " << printDistFunction << std::endl;
    mirrorData << "PrintEnergyWeightFunction  = " << printEnergyWeightFunction << std::endl;
    mirrorData << "PrintPressure              = " << printPressure << std::endl;
    mirrorData << "PrintRealPressure          = " << printRealPressure << std::endl;
    mirrorData << "PrintVorticity             = " << printVorticity << std::endl;
    mirrorData << "PrintInnerNormal           = " << printInnerNormal << std::endl;
    mirrorData << "PrintMeshVelocity          = " << printMeshVelocity << std::endl;
    mirrorData << "PrintMeshDisplacement      = " << printMeshDisplacement << std::endl;
    mirrorData << "PrintGlueZone              = " << printGlueZone << std::endl;
    mirrorData << "PrintJacobian              = " << printJacobian << std::endl;
    mirrorData << "PrintProcess               = " << printProcess << std::endl;
    mirrorData << "PrintLines                 = " << printLines << std::endl << std::endl;

    int dimension = 2;

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++READIN MESH+++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //+++++++++++++++++++++++++++PHYSICAL ENTITIES++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    int number_physical_entities;
    file >> number_physical_entities;
    std::getline(file, line);
    std::unordered_map<int, std::string> physicalEntities;
    physicalEntities.reserve(number_physical_entities);

    for (int i = 0; i < number_physical_entities; i++)
    {
        std::getline(file, line);
        std::vector<std::string> tokens = split(line, " ");
        int index;
        std::istringstream(tokens[1]) >> index;
        physicalEntities[index] = tokens[2].substr(1, tokens[2].size() - 2);
    }
    std::getline(file, line); std::getline(file, line);

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //+++++++++++++++++++++++++++++++++NODES++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    file >> numNodes;
    nodes_.reserve(numNodes);
    std::getline(file, line);
    int index = 0;
    if (rank == 0) std::cout << "Number of Nodes " << " " << numNodes << std::endl;
    for (int i = 0; i < numNodes; i++)
    {
        double x[2];
        std::getline(file, line);
        std::vector<std::string> tokens = split(line, " ");
        std::istringstream(tokens[1]) >> x[0];
        std::istringstream(tokens[2]) >> x[1];
        //addNode(i, coord);
         Node *node = new Node(x, index++);
         nodes_.push_back(node);
    }
    std::getline(file, line); std::getline(file, line);

    mirrorData << "Nodal Coordinates " << numNodes << std::endl;
    for (int i = 0 ; i<numNodes; i++){
        
        double* x = nodes_[i]->getCoordinates();       
        for (int j=0; j<2; j++){
            mirrorData << x[j] << " ";
        };
        mirrorData << std::endl;
        nodes_[i] -> setVelocity(velocityInf);
        nodes_[i] -> setPreviousVelocity(velocityInf);
        double u[2];
        u[0] = 0.; u[1] = 0.;
        nodes_[i] -> setMeshVelocity(u);

        nodes_[i] -> setPreviousMeshVelocityComponent(0,0.);
        nodes_[i] -> setPreviousMeshVelocityComponent(1,0.);

        nodes_[i] -> setPreviousCoordinates(0,x[0]);
        nodes_[i] -> setPreviousCoordinates(1,x[1]);

    };


    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++ELEMENTS++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    int number_elements;
    file >> number_elements;
    elements_.reserve(number_elements);
    boundary_.reserve(number_elements/10);
    index = 0;
    std::getline(file, line);
    int cont = 0;

    numBoundElems = 0;
    numElem = 0;
    numFSIInterfaces = 0;//É melhor setar esse valor no main

    std::vector<BoundaryCondition*> dirichlet, neumann, glue, FSinterface;
    dirichlet = geometry_->getBoundaryCondition("DIRICHLET"); 
    neumann = geometry_->getBoundaryCondition("NEUMANN"); 
    glue = geometry_->getBoundaryCondition("GLUE"); 
    FSinterface = geometry_->getBoundaryCondition("FSINTERFACE"); 

    numFSIInterfaces = FSinterface.size();

    for (int i = 0; i < number_elements; i++)
    {
        std::getline(file, line);
        std::vector<std::string> tokens = split(line, " ");
        std::vector<int> values(tokens.size(), 0);
        for (size_t j = 0; j < tokens.size(); j++)
            std::istringstream(tokens[j]) >> values[j];
        std::string elementType = gmshElement[values[1]];
        int number_nodes_per_element = numNodes2[elementType];
        std::vector<int> elementNodes;
        elementNodes.reserve(number_nodes_per_element);

        for (size_t j = 5 ; j < values.size(); j++)
            elementNodes.push_back(values[j]-1);
 
        std::string name = physicalEntities[values[3]];

        //Adding 2D elements to surfaces
        if (name[0] == 's'){
            if (supportedElements.find(elementType) == supportedElements.end()){
                std::cout << elementType << " is not supported.\n";
                exit(EXIT_FAILURE);
            }

            PlaneSurface* object = geometry_ -> getPlaneSurface(name);
            int materialIndex = object -> getMaterial() -> getIndex();
            double thickness = object -> getThickness();
            numElem++;

            int *connect;
            connect = new int[6];
            for (int j = 0 ; j < 6; j++) connect[j] = elementNodes[j];
           


            Elements *el = new Elements(index++,connect,nodes_,fluidParameters);
            elements_.push_back(el);
            for (int k = 0; k<6; k++){
                nodes_[connect[k]] -> pushInverseIncidence(index);
            };
        }
        else if (name[0] == 'l')
        {
            int connectB[3];

            connectB[0] = elementNodes[0];
            connectB[1] = elementNodes[1];
            connectB[2] = elementNodes[2];

            int ibound;
            std::string::size_type sz;   // alias of size_t
            ibound = std::stoi (&name[1],nullptr,10);
            
            int constrain[3];
            double value[3];


            for (int i = 0; i < dirichlet.size(); i++){
                if (name == dirichlet[i] -> getLineName()){
                    if ((dirichlet[i] -> getComponentX()).size() == 0){
                        constrain[0] = 0; value[0] = 0;
                    }else{
                        std::vector<double> c = dirichlet[i] -> getComponentX();
                        constrain[0] = 1;
                        value[0] = c[0];
                    }
                    if ((dirichlet[i] -> getComponentY()).size() == 0){
                        constrain[1] = 0; value[1] = 0;
                    }else{
                        std::vector<double> c = dirichlet[i] -> getComponentY();
                        constrain[1] = 1;
                        value[1] = c[0];
                    }
                }
            }

            for (int i = 0; i < neumann.size(); i++){
                if (name == neumann[i] -> getLineName()){
                    if ((neumann[i] -> getComponentX()).size() == 0){
                        constrain[0] = 0; value[0] = 0;
                    }else{
                        std::vector<double> c = neumann[i] -> getComponentX();
                        constrain[0] = 0;
                        value[0] = c[0];
                    }
                    if ((neumann[i] -> getComponentY()).size() == 0){
                        constrain[1] = 0; value[1] = 0;
                    }else{
                        std::vector<double> c = neumann[i] -> getComponentX();
                        constrain[1] = 0;
                        value[1] = c[0];
                    }
                }
            }  
            
            for (int i = 0; i < glue.size(); i++){
                if (name == glue[i] -> getLineName()){
                    if ((glue[i] -> getComponentX()).size() == 0){
                        constrain[0] = 2; value[0] = 0;
                    }else{
                        std::vector<double> c = glue[i] -> getComponentX();
                        constrain[0] = 2;
                        value[0] = c[0];
                    }
                    if ((glue[i] -> getComponentY()).size() == 0){
                        constrain[1] = 2; value[1] = 0;
                    }else{
                        std::vector<double> c = glue[i] -> getComponentY();
                        constrain[1] = 2;
                        value[1] = c[0];
                    }
                }
            }              

            for (int i = 0; i < FSinterface.size(); i++){
                if (name == FSinterface[i] -> getLineName()){
                    if ((FSinterface[i] -> getComponentX()).size() == 0){
                        constrain[0] = 3; value[0] = 0;
                    }else{
                        std::vector<double> c = FSinterface[i] -> getComponentX();
                        constrain[0] = 3;
                        value[0] = c[0];
                    }
                    if ((FSinterface[i] -> getComponentY()).size() == 0){
                        constrain[1] = 3; value[1] = 0;
                    }else{
                        std::vector<double> c = FSinterface[i] -> getComponentY();
                        constrain[1] = 3;
                        value[1] = c[0];
                    }
                }
            }        
            Boundaries * bound = new Boundaries(connectB, numBoundElems++, constrain, value, ibound);
            boundary_.push_back(bound); 
        }   
    }

    if (rank == 0) std::cout << "Number of elements " << number_elements << " " 
                             << numElem << " " << numBoundElems << std::endl;
    mirrorData << std::endl << "Element Connectivity" << std::endl;        
    
    for (int jel=0; jel<numElem; jel++){
        int *connec = elements_[jel] -> getConnectivity();       
        for (int i=0; i<4*dimension-2; i++){
            mirrorData << connec[i] << " ";
        };
        mirrorData << std::endl;
    };


    // for (int i = 0; i < numNodes; ++i){
    //     typename Node::VecLocD x;
    //     x = nodes_[i]->getCoordinates();    
    //     if ((x(0) == 0) && (x(1) == 0)) {
    //         nodes_[i] -> setConstrainsLaplace(0,1,0);
    //         nodes_[i] -> setConstrainsLaplace(0,1,0);
    //     }
    // }
    




    //Sets boundary constrains
    for (int ibound = 0; ibound < numBoundElems; ibound++){
        
        int* connectB = boundary_[ibound] -> getBoundaryConnectivity();
        int no1 = connectB[0];
        int no2 = connectB[1];
        int no3 = connectB[2];
        //if ((boundary_[ibound] -> getConstrain(0) != 2) && (boundary_[ibound] -> getConstrain(0) != 0)){
        if (boundary_[ibound] -> getConstrain(0) == 3){
            nodes_[no1] -> setConstrainsLaplace(0,1,0);
            nodes_[no2] -> setConstrainsLaplace(0,1,0);
            nodes_[no3] -> setConstrainsLaplace(0,1,0);
        };
        // if ((boundary_[ibound] -> getConstrain(1) != 2) && (boundary_[ibound] -> getConstrain(1) != 0)){
        if (boundary_[ibound] -> getConstrain(1) == 3){
            nodes_[no1] -> setConstrainsLaplace(1,1,0);
            nodes_[no2] -> setConstrainsLaplace(1,1,0);
            nodes_[no3] -> setConstrainsLaplace(1,1,0);
        };


        
        if ((boundary_[ibound] -> getConstrain(0) == 1) || (boundary_[ibound] -> getConstrain(0) == 3)){

            //Desfazer primeira parte do if para voltar a cond. cont. constante
            // if (boundary_[ibound] -> getConstrainValue(0) <= 1.){
                
            //     typename Node::VecLocD x;
            //     x = nodes_[no1]->getCoordinates();                 
                
            //     nodes_[no1] -> setConstrains(0,boundary_[ibound] -> 
            //                                  getConstrain(0),
            //                                  x(1) * boundary_[ibound] -> 
            //                                  getConstrainValue(0));
                
            //     x = nodes_[no2]->getCoordinates();                 
                
            //     nodes_[no2] -> setConstrains(0,boundary_[ibound] -> 
            //                                  getConstrain(0),
            //                                  x(1) * boundary_[ibound] -> 
            //                                  getConstrainValue(0));
                
            //     x = nodes_[no3]->getCoordinates();                 
                
            //     nodes_[no3] -> setConstrains(0,boundary_[ibound] -> 
            //                                  getConstrain(0),
            //                                  x(1) * boundary_[ibound] -> 
            //                                  getConstrainValue(0));
            //     //ate aqui
            // } else {
            nodes_[no1] -> setConstrains(0,boundary_[ibound] -> getConstrain(0),
                                     boundary_[ibound] -> getConstrainValue(0));
            nodes_[no2] -> setConstrains(0,boundary_[ibound] -> getConstrain(0),
                                     boundary_[ibound] -> getConstrainValue(0));
            nodes_[no3] -> setConstrains(0,boundary_[ibound] -> getConstrain(0),
                                     boundary_[ibound] -> getConstrainValue(0));
             // };
        };

        if((boundary_[ibound] -> getConstrain(1) == 1) || (boundary_[ibound] -> getConstrain(1) == 3)){
            nodes_[no1] -> setConstrains(1,boundary_[ibound] -> getConstrain(1),
                                     boundary_[ibound] -> getConstrainValue(1));
            nodes_[no2] -> setConstrains(1,boundary_[ibound] -> getConstrain(1),
                                     boundary_[ibound] -> getConstrainValue(1));
            nodes_[no3] -> setConstrains(1,boundary_[ibound] -> getConstrain(1),
                                     boundary_[ibound] -> getConstrainValue(1));
        };     
    };

























        //Print nodal constrains
    for (int i=0; i<numNodes; i++){

        mirrorData<< "Constrains " << i
                  << " " << nodes_[i] -> getConstrains(0)
                  << " " << nodes_[i] -> getConstrainValue(0)
                  << " " << nodes_[i] -> getConstrains(1)
                  << " " << nodes_[i] -> getConstrainValue(1) << std::endl;
    }; 

    for (int i=0; i<numBoundElems; i++){

        mirrorData<< "Bound Elements " << i
                  << " " << boundary_[i] -> getBoundaryGroup() << std::endl;
    }; 


    //Sets fluid elements and sides on interface boundaries
    for (int i=0; i<numBoundElems; i++){
        if ((boundary_[i] -> getConstrain(0) > 0) ||
            (boundary_[i] -> getConstrain(1) > 0)) {

            int* connectB = boundary_[i] -> getBoundaryConnectivity();

            for (int j=0; j<numElem; j++){
                int *connect = elements_[j] -> getConnectivity();
                
                int flag = 0;
                int side[3];
                for (int k=0; k<6; k++){
                    if ((connectB[0] == connect[k]) || 
                        (connectB[1] == connect[k]) ||
                        (connectB[2] == connect[k])){
                        side[flag] = k;
                        flag++;
                    };
                };
                
                if (flag == 3){
                    boundary_[i] -> setElement(j);
                    //Sets element index and side
                    if ((side[0]==4) || (side[1]==4) || (side[2]==4)){
                        boundary_[i] -> setElementSide(0);
                        elements_[boundary_[i]->getElement()] -> setElemSideInBoundary(0);
                    };
                    if ((side[0]==5) || (side[1]==5) || (side[2]==5)){
                        boundary_[i] -> setElementSide(1);
                        elements_[boundary_[i]->getElement()] -> setElemSideInBoundary(1);
                    };
                    if ((side[0]==3) || (side[1]==3) || (side[2]==3)){
                        boundary_[i] -> setElementSide(2);
                        elements_[boundary_[i]->getElement()] -> setElemSideInBoundary(2);
                    };
                };
            };
        };
    };

    domainDecompositionMETIS();

    iAux = 0;


    //Closing the file
    file.close();
    if (deleteFiles)
        system((remove + inputFile).c_str());


return;
};

//------------------------------------------------------------------------------
//-------------------------SOLVE TRANSIENT FLUID PROBLEM------------------------
//------------------------------------------------------------------------------
template<>
int Fluid<2>::getInitialAcceleration() {

    Vec               b, u, All;
    PetscErrorCode    ierr;
    PetscInt          Ii, Ione, iterations;
    KSP               ksp;
    PC                pc;
    VecScatter        ctx;
    PetscScalar       val;
    //IS             rowperm       = NULL,colperm = NULL;
    //    MatNullSpace      nullsp;
   
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    double &alpha_f = fluidParameters.getAlphaF();
    double &alpha_m = fluidParameters.getAlphaM();
    double &gamma = fluidParameters.getGamma();

    if (rank == 0) {std::cout << "---------------- COMPUTING INITIAL ACCELERATION "
                              << "----------------"
                              << std::endl;}
        
    ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                        2*numNodes+numNodes, 2*numNodes+numNodes,
                        60,NULL,100,NULL,&A); CHKERRQ(ierr);
    
    //Create PETSc vectors
    ierr = VecCreate(PETSC_COMM_WORLD,&b); CHKERRQ(ierr);
    ierr = VecSetSizes(b,PETSC_DECIDE,2*numNodes+numNodes); CHKERRQ(ierr);
    ierr = VecSetFromOptions(b); CHKERRQ(ierr);
    ierr = VecDuplicate(b,&u); CHKERRQ(ierr);
    ierr = VecDuplicate(b,&All); CHKERRQ(ierr);
    
    for (int jel = 0; jel < numElem; jel++){   
        
        //if (part_elem[jel] == rank) {
            //Compute Element matrix
            int *connec = elements_[jel] -> getConnectivity();

            double **matrix;
            double rhs[18] = {};
            matrix = new double*[18]();
            for (int i = 0; i < 18; i++) matrix[i] = new double[18]();

            elements_[jel] -> getTransientNavierStokesInitial(rhs,matrix);
            //Disperse local contributions into the global matrix
            //Matrix K and C
            for (int i=0; i<6; i++){
                for (int j=0; j<6; j++){
                    int dof_i = 2 * connec[i];
                    int dof_j = 2 * connec[j];
                    MatSetValues(A,1,&dof_i,1, &dof_j,&matrix[2*i  ][2*j  ],ADD_VALUES);
                    
                    dof_i = 2 * connec[i] + 1;
                    dof_j = 2 * connec[j];
                    MatSetValues(A, 1, &dof_i, 1, &dof_j,&matrix[2*i+1][2*j  ],ADD_VALUES);
                    
                    dof_i = 2 * connec[i];
                    dof_j = 2 * connec[j] + 1;
                    MatSetValues(A, 1, &dof_i, 1, &dof_j,&matrix[2*i  ][2*j+1],ADD_VALUES);
                    
                    dof_i = 2 * connec[i] + 1;
                    dof_j = 2 * connec[j] + 1;
                    MatSetValues(A, 1, &dof_i, 1, &dof_j,&matrix[2*i+1][2*j+1],ADD_VALUES);
                
                    //Matrix Q and Qt
                    dof_i = 2 * connec[i];
                    dof_j = 2 * numNodes + connec[j];
                    MatSetValues(A, 1, &dof_i, 1, &dof_j,&matrix[2*i  ][12+j],ADD_VALUES);
                    
                    dof_i = 2 * connec[i];
                    dof_j = 2 * numNodes + connec[j];
                    MatSetValues(A, 1, &dof_j, 1, &dof_i,&matrix[12+j][2*i  ],ADD_VALUES);
                    
                    dof_i = 2 * connec[i] + 1;
                    dof_j = 2 * numNodes + connec[j];
                    MatSetValues(A, 1, &dof_i, 1, &dof_j,&matrix[2*i+1][12+j],ADD_VALUES);
                    
                    dof_i = 2 * connec[i] + 1;
                    dof_j = 2 * numNodes + connec[j];
                    MatSetValues(A, 1, &dof_j, 1, &dof_i,&matrix[12+j][2*i+1],ADD_VALUES);
                    
                    dof_i = 2 * numNodes + connec[i];
                    dof_j = 2 * numNodes + connec[j];
                    MatSetValues(A, 1, &dof_i, 1, &dof_j,&matrix[12+i][12+j],ADD_VALUES);
                };
                
                //Rhs vector
                int dof_i = 2 * connec[i];
                VecSetValues(b, 1, &dof_i, &rhs[2*i  ],ADD_VALUES);

                dof_i = 2 * connec[i]+1;
                VecSetValues(b, 1, &dof_i, &rhs[2*i+1],ADD_VALUES);
                
                dof_i = 2 * numNodes + connec[i];
                VecSetValues(b, 1, &dof_i, &rhs[12+i],ADD_VALUES);
            };
            for (int i = 0; i < 18; ++i) delete [] matrix[i];
            delete [] matrix;
        //};
    }; //Elements
    
    //Assemble matrices and vectors
    ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    
    ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(b);CHKERRQ(ierr);
    

    // Mat Aperm;
    // MatGetOrdering(A,MATORDERINGRCM,&rowperm,&colperm);
    // MatPermute(A,rowperm,colperm,&Aperm);
    // VecPermute(b,colperm,PETSC_FALSE);
    // MatDestroy(&A);
    // A    = Aperm;    

    // MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    //MatView(A,PETSC_VIEWER_DRAW_WORLD);CHKERRQ(ierr);
    //ierr = VecView(b,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);


    // PetscViewer viewer;
    // PetscViewerASCIIOpen(PETSC_COMM_WORLD, "Amat.m", &viewer);
    // PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
    // MatView(A,viewer);
    // PetscViewerPopFormat(viewer);
    // PetscViewerDestroy(&viewer);
    
    //Create KSP context to solve the linear system
    ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
    
    ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
    
    // ierr = KSPSetTolerances(ksp,1.e-10,PETSC_DEFAULT,PETSC_DEFAULT,
    //                         500);CHKERRQ(ierr);
    
    // ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
    
    // // ierr = KSPGetPC(ksp,&pc);
    
    // // ierr = PCSetType(pc,PCNONE);
    
    // // ierr = KSPSetType(ksp,KSPDGMRES); CHKERRQ(ierr);

    // ierr = KSPGMRESSetRestart(ksp, 500); CHKERRQ(ierr);
    
    // //    ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
    

// //   //   ierr = MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,NULL,&nullsp);
// // // ierr = MatSetNullSpace(A, nullsp);
// // // ierr = MatNullSpaceDestroy(&nullsp);



#if defined(PETSC_HAVE_MUMPS)
    ierr = KSPSetType(ksp,KSPPREONLY);
    ierr = KSPGetPC(ksp,&pc);
    ierr = PCSetType(pc, PCLU);
#endif          
    ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
    ierr = KSPSetUp(ksp);



    ierr = KSPSolve(ksp,b,u);CHKERRQ(ierr);

    ierr = KSPGetTotalIterations(ksp, &iterations);            

    // VecPermute(u,rowperm,PETSC_TRUE);

    // ierr = VecView(u,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKERRQ(ierr);
    
    //Gathers the solution vector to the master process
    ierr = VecScatterCreateToAll(u, &ctx, &All);CHKERRQ(ierr);
    
    ierr = VecScatterBegin(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);
    CHKERRQ(ierr);
    
    ierr = VecScatterEnd(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);
    CHKERRQ(ierr);
    
    ierr = VecScatterDestroy(&ctx);CHKERRQ(ierr);
    
    //Updates nodal values
    double p_;
    double dpNorm = 0.;
    Ione = 1;
    double u_ = 0.;
    double uPrev = 0.;
    double pPrev = 0.;
    
    for (int i = 0; i < numNodes; ++i){
        
        Ii = 2*i;
        ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
        // duNorm += val*val;
        u_ = nodes_[i] -> getVelocity(0);
        uPrev = u_ - val * dTime * gamma;
        nodes_[i] -> setPreviousVelocityComponent(0,uPrev);
        nodes_[i] -> setAccelerationComponent(0,val);
        
        // std::cout << "Velocity " << i << " " << uPrev << " " << u_ << " " << val << std::endl;


        Ii = 2*i+1;
        ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
        // duNorm += val*val;
        u_ = nodes_[i] -> getVelocity(1);
        uPrev = u_ - val * dTime * gamma;
        nodes_[i] -> setPreviousVelocityComponent(1,uPrev);
        nodes_[i] -> setAccelerationComponent(1,val);
    };
    
    for (int i = 0; i<numNodes; i++){
        Ii = 2*numNodes+i;
        ierr = VecGetValues(All,Ione,&Ii,&val);CHKERRQ(ierr);
        dpNorm += val*val;
        // std::cout << "Pressure " << i << " " << val << std::endl;
        nodes_[i] -> setPressure(val);
    };
    
    //Computes the solution vector norm
    //ierr = VecNorm(u,NORM_2,&val);CHKERRQ(ierr);
 
              
    ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
    ierr = VecDestroy(&b); CHKERRQ(ierr);
    ierr = VecDestroy(&u); CHKERRQ(ierr);
    ierr = VecDestroy(&All); CHKERRQ(ierr);
    ierr = MatDestroy(&A); CHKERRQ(ierr);
        
    printResults(1000);

    
    return 0;
};




//------------------------------------------------------------------------------
//-------------------------SOLVE STEADY LAPLACE PROBLEM-------------------------
//------------------------------------------------------------------------------
template<>
int Fluid<2>::solveSteadyLaplaceProblem(int iterNumber, double tolerance) {

    Vec               b, u, All;
    PetscErrorCode    ierr;
    PetscInt          Istart, Iend, Ii, Ione, iterations;
    KSP               ksp;
    PC                pc;
    VecScatter        ctx;
    PetscScalar       val;
   
    int rank;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);        

    for (int inewton = 0; inewton < iterNumber; inewton++){

        ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                            2*numNodes, 2*numNodes,
                            100,NULL,100,NULL,&A); CHKERRQ(ierr);
        
        ierr = MatGetOwnershipRange(A, &Istart, &Iend);CHKERRQ(ierr);
        
        //Create PETSc vectors
        ierr = VecCreate(PETSC_COMM_WORLD,&b);CHKERRQ(ierr);
        ierr = VecSetSizes(b,PETSC_DECIDE,2*numNodes);CHKERRQ(ierr);
        ierr = VecSetFromOptions(b);CHKERRQ(ierr);
        ierr = VecDuplicate(b,&u);CHKERRQ(ierr);
        ierr = VecDuplicate(b,&All);CHKERRQ(ierr);
        
        //std::cout << "Istart = " << Istart << " Iend = " << Iend << std::endl;
        
        for (int jel = 0; jel < numElem; jel++){               
            //Compute Element matrix
            int *connec = elements_[jel] -> getConnectivity();

            double **matrix;
            double rhs[18] = {};
            matrix = new double*[18]();
            for (int i = 0; i < 18; i++) matrix[i] = new double[18]();

            elements_[jel] -> getSteadyLaplace2(matrix,rhs);
            
            //Disperse local contributions into the global matrix
            //Matrix K and C
            for (int i=0; i<6; i++){
                for (int j=0; j<6; j++){
                    if (fabs(matrix[2*i  ][2*j  ]) >= 1.e-8){
                        int dof_i = 2*connec[i];
                        int dof_j = 2*connec[j];
                        ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&matrix[2*i  ][2*j  ],ADD_VALUES);
                    };
                    if (fabs(matrix[2*i+1][2*j  ]) >= 1.e-8){
                        int dof_i = 2*connec[i]+1;
                        int dof_j = 2*connec[j];
                        ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&matrix[2*i+1][2*j  ],ADD_VALUES);
                    };
                    if (fabs(matrix[2*i  ][2*j+1]) >= 1.e-8){
                        int dof_i = 2*connec[i];
                        int dof_j = 2*connec[j]+1;
                        ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&matrix[2*i  ][2*j+1],ADD_VALUES);
                    };
                    if (fabs(matrix[2*i+1][2*j+1]) >= 1.e-8){
                        int dof_i = 2*connec[i]+1;
                        int dof_j = 2*connec[j]+1;
                        ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&matrix[2*i+1][2*j+1],ADD_VALUES);
                    };
                };
                                    
                //Rhs vector
                if (fabs(rhs[2*i  ]) >= 1.e-8){
                    int dof_i = 2*connec[i];
                    ierr = VecSetValues(b,1,&dof_i,&rhs[2*i  ],ADD_VALUES);
                };
                
                if (fabs(rhs[2*i+1]) >= 1.e-8){
                    int dof_i = 2*connec[i]+1;
                    ierr = VecSetValues(b,1,&dof_i,&rhs[2*i+1],ADD_VALUES);
                };
            };
            for (int i = 0; i < 18; ++i) delete [] matrix[i];
            delete [] matrix;
        };
        
        //Assemble matrices and vectors
        ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        
        ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
        ierr = VecAssemblyEnd(b);CHKERRQ(ierr);
        
        //MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
        //ierr = VecView(b,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
        
        //Create KSP context to solve the linear system
        ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
        
        ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
        
#if defined(PETSC_HAVE_MUMPS)
        ierr = KSPSetType(ksp,KSPPREONLY);
        ierr = KSPGetPC(ksp,&pc);
        ierr = PCSetType(pc, PCLU);
#endif
        
        ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
        ierr = KSPSetUp(ksp);
        
        
        
        ierr = KSPSolve(ksp,b,u);CHKERRQ(ierr);
        
        ierr = KSPGetTotalIterations(ksp, &iterations);

        //std::cout << "GMRES Iterations = " << iterations << std::endl;
        
        //ierr = VecView(u,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKERRQ(ierr);
        
        //Gathers the solution vector to the master process
        ierr = VecScatterCreateToAll(u, &ctx, &All);CHKERRQ(ierr);
        
        ierr = VecScatterBegin(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);CHKERRQ(ierr);
        
        ierr = VecScatterEnd(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);CHKERRQ(ierr);
        
        ierr = VecScatterDestroy(&ctx);CHKERRQ(ierr);
                
        //Updates nodal values
        double u_ [2];
        Ione = 1;

        for (int i = 0; i < numNodes; ++i){
            Ii = 2*i;
            ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
            u_[0] = val;
            Ii = 2*i+1;
            ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
            u_[1] = val;
            if (nodes_[i] -> getConstrainsLaplace(0) != 1) nodes_[i] -> incrementCoordinate(0,u_[0]);
            if (nodes_[i] -> getConstrainsLaplace(1) != 1) nodes_[i] -> incrementCoordinate(1,u_[1]);




            // double x = nodes_[i] -> getCoordinateValue(0);
            // double y = nodes_[i] -> getCoordinateValue(1);
            // double xp = nodes_[i] -> getPreviousCoordinateValue(0);
            // double yp = nodes_[i] -> getPreviousCoordinateValue(1);
            // double vx = nodes_[i] -> getMeshVelocity(0);
            // double vy = nodes_[i] -> getMeshVelocity(1);
            // double ax = nodes_[i] -> getMeshAcceleration(0);
            // double ay = nodes_[i] -> getMeshAcceleration(1);

            // double accelx = (x - xp) / (0.25 * dTime * dTime) - vx / (0.25 * dTime) - ax * (0.5/0.25 - 1.0);
            // double accely = (y - yp) / (0.25 * dTime * dTime) - vy / (0.25 * dTime) - ay * (0.5/0.25 - 1.0);

            // nodes_[i] -> setMeshAccelerationComponent(0,accelx);
            // nodes_[i] -> setMeshAccelerationComponent(1,accely);

            // double velx = 0.5 * dTime * accelx + vx + dTime * (1.0 - 0.5) * ax;
            // double vely = 0.5 * dTime * accely + vy + dTime * (1.0 - 0.5) * ay;
            // // std::cout << "asd as " << vx << " " << velx << " " << accelx << " " << ax << std::endl;
            // nodes_[i] -> setMeshVelocityComponent(0,velx);
            // nodes_[i] -> setMeshVelocityComponent(1,vely);

        };
        
        //Computes the solution vector norm
        ierr = VecNorm(u,NORM_2,&val);CHKERRQ(ierr);
        
        if(rank == 0){
            std::cout << "MESH MOVING - ERROR = " << val 
                      << std::scientific <<  std::endl;
        };

        ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
        ierr = VecDestroy(&b); CHKERRQ(ierr);
        ierr = VecDestroy(&u); CHKERRQ(ierr);
        ierr = VecDestroy(&All); CHKERRQ(ierr);
        ierr = MatDestroy(&A); CHKERRQ(ierr);

        if(val <= tolerance){
            break;
        };
    };
    
    // for (int i=0; i<numElem; i++){
    //     elements_[i] -> computeNodalGradient();            
    // };

    if (rank == 0) {
        //Computing velocity divergent
        //      printResults(1);
    };

    return 0;
};

//------------------------------------------------------------------------------
//-------------------------SOLVE TRANSIENT FLUID PROBLEM------------------------
//------------------------------------------------------------------------------
template<>
void Fluid<2>::readInitialValues(const std::string& inputPrev, const std::string& inputCurr) {

    int rank;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    hid_t filePrevious;
    hid_t fileCurrent;
    hid_t dataset;
    herr_t status;

    double *vecValues;
    double *scaValues;
    vecValues = new double[3*numNodes];
    scaValues = new double[numNodes];

    filePrevious = H5Fopen(inputPrev.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
    fileCurrent = H5Fopen(inputCurr.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);

    char datasetName[] = "/velocity";
    dataset = H5Dopen( filePrevious, datasetName, H5P_DEFAULT );
    status = H5Dread( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vecValues[0]);
    status = H5Dclose(dataset);
    for (int i = 0; i < numNodes; ++i){
        nodes_[i] -> setPreviousVelocityComponent(0,vecValues[3*i  ]);
        nodes_[i] -> setPreviousVelocityComponent(1,vecValues[3*i+1]);
    }
    dataset = H5Dopen( fileCurrent, datasetName, H5P_DEFAULT );
    status = H5Dread( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vecValues[0]);
    status = H5Dclose(dataset);
    for (int i = 0; i < numNodes; ++i){
        nodes_[i] -> setVelocityComponent(0,vecValues[3*i  ]);
        nodes_[i] -> setVelocityComponent(1,vecValues[3*i+1]);
    }

    char datasetName2[] = "/acceleration";
    dataset = H5Dopen( filePrevious, datasetName2, H5P_DEFAULT );
    status = H5Dread( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vecValues[0]);
    status = H5Dclose(dataset);
    for (int i = 0; i < numNodes; ++i){
        nodes_[i] -> setPreviousAccelerationComponent(0,vecValues[3*i  ]);
        nodes_[i] -> setPreviousAccelerationComponent(1,vecValues[3*i+1]);
    }
    dataset = H5Dopen( fileCurrent, datasetName2, H5P_DEFAULT );
    status = H5Dread( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vecValues[0]);
    status = H5Dclose(dataset);
    for (int i = 0; i < numNodes; ++i){
        nodes_[i] -> setAccelerationComponent(0,vecValues[3*i  ]);
        nodes_[i] -> setAccelerationComponent(1,vecValues[3*i+1]);
    }

    char datasetName3[] = "/lagrangeMultiplers";
    dataset = H5Dopen( fileCurrent, datasetName3, H5P_DEFAULT );
    status = H5Dread( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vecValues[0]);
    status = H5Dclose(dataset);
    for (int i = 0; i < numNodes; ++i){
        nodes_[i] -> setLagrangeMultiplier(0,vecValues[3*i  ]);
        nodes_[i] -> setLagrangeMultiplier(1,vecValues[3*i+1]);
    }

    char datasetName4[] = "/pressure";
    dataset = H5Dopen( fileCurrent, datasetName4, H5P_DEFAULT );
    status = H5Dread( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &scaValues[0]);
    status = H5Dclose(dataset);
    for (int i = 0; i < numNodes; ++i){
        nodes_[i] -> setPressure(vecValues[i]);
    }

    delete [] vecValues;
    delete [] scaValues;

    //End HDF5 file
    status = H5Fclose(filePrevious);

    return;
}


//------------------------------------------------------------------------------
//-------------------------SOLVE TRANSIENT FLUID PROBLEM------------------------
//------------------------------------------------------------------------------
template<>
int Fluid<2>::solveFSIFluid(int iterNumber, double tolerance, int problem_type){

    Vec               b, u, All, Allu;
    PetscErrorCode    ierr;
    PetscInt          Ii, Ione, iterations;
    KSP               ksp;
    PC                pc;
    VecScatter        ctx;
    PetscScalar       val;
    //    MatNullSpace      nullsp;

    int rank;

    iAux++;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    //Check if the problem type can be computed
    if ((problem_type > 2) || (problem_type < 1)){
        std::cout << "WRONG PROBLEM TYPE." << std::endl;
        return 0;
    };
            
    double duNorm=100.;

    double &alpha_f = fluidParameters.getAlphaF();
    double &alpha_m = fluidParameters.getAlphaM();
    double &gamma = fluidParameters.getGamma();
         
    for (int inewton = 0; inewton < iterNumber; inewton++){
        
        std::clock_t t1 = std::clock();
        
        ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                            2*numNodes+numNodes, 2*numNodes+numNodes,
                            100,NULL,300,NULL,&A); 
        CHKERRQ(ierr);
                
        //Create PETSc vectors
        ierr = VecCreate(PETSC_COMM_WORLD,&b);CHKERRQ(ierr);
        ierr = VecSetSizes(b,PETSC_DECIDE,2*numNodes+numNodes);
        CHKERRQ(ierr);
        ierr = VecSetFromOptions(b);CHKERRQ(ierr);
        ierr = VecDuplicate(b,&u);CHKERRQ(ierr);
        ierr = VecDuplicate(b,&All);CHKERRQ(ierr);
        
        //std::cout << "Istart = " << Istart << " Iend = " << Iend << std::endl;

        for (int jel = 0; jel < numElem; jel++){   
            
            if (part_elem[jel] == rank) {
                //Compute Element matrix
                int *connec = elements_[jel] -> getConnectivity();

                double **matrix;
                double rhs[18] = {};
                matrix = new double*[18]();
                for (int i = 0; i < 18; i++) matrix[i] = new double[18]();

                elements_[jel] -> getTransientNavierStokes(matrix,rhs);
                
                //Disperse local contributions into the global matrix
                //Matrix K and C
                for (int i=0; i<6; i++){
                    for (int j=0; j<6; j++){
                        int dof_i = 2 * connec[i];
                        int dof_j = 2 * connec[j];
                        MatSetValues(A, 1, &dof_i,1, &dof_j, &matrix[2*i  ][2*j  ], ADD_VALUES);
                        
                        dof_i = 2 * connec[i] + 1;
                        dof_j = 2 * connec[j];
                        MatSetValues(A, 1, &dof_i, 1, &dof_j, &matrix[2*i+1][2*j  ], ADD_VALUES);
                        
                        dof_i = 2 * connec[i];
                        dof_j = 2 * connec[j] + 1;
                        MatSetValues(A, 1, &dof_i, 1, &dof_j, &matrix[2*i  ][2*j+1], ADD_VALUES);

                        dof_i = 2 * connec[i] + 1;
                        dof_j = 2 * connec[j] + 1;
                        MatSetValues(A, 1, &dof_i, 1, &dof_j, &matrix[2*i+1][2*j+1], ADD_VALUES);
                        
                        //Matrix Q and Qt
                        dof_i = 2 * connec[i];
                        dof_j = 2 * numNodes + connec[j];
                        MatSetValues(A, 1, &dof_i, 1, &dof_j, &matrix[2*i  ][12+j], ADD_VALUES);
                    
                        dof_i = 2 * connec[i];
                        dof_j = 2 * numNodes + connec[j];
                        MatSetValues(A, 1, &dof_j, 1, &dof_i, &matrix[12+j][2*i  ], ADD_VALUES);
                        
                        dof_i = 2 * connec[i] + 1;
                        dof_j = 2 * numNodes + connec[j];
                        MatSetValues(A, 1, &dof_i, 1, &dof_j, &matrix[2*i+1][12+j], ADD_VALUES);
                        
                        dof_i = 2 * connec[i] + 1;
                        dof_j = 2 * numNodes + connec[j];
                        MatSetValues(A, 1, &dof_j, 1, &dof_i, &matrix[12+j][2*i+1], ADD_VALUES);
                        
                        dof_i = 2 * numNodes + connec[i];
                        dof_j = 2 * numNodes + connec[j];
                        MatSetValues(A, 1, &dof_i, 1, &dof_j, &matrix[12+i][12+j], ADD_VALUES);
                    };
                    
                    //Rhs vector
                    int dof_i = 2 * connec[i];
                    VecSetValues(b, 1, &dof_i, &rhs[2*i  ], ADD_VALUES);
                    
                    dof_i = 2 * connec[i]+1;
                    VecSetValues(b, 1, &dof_i, &rhs[2*i+1], ADD_VALUES);
                    
                    dof_i = 2 * numNodes + connec[i];
                    VecSetValues(b, 1, &dof_i, &rhs[12+i], ADD_VALUES);
                };
                for (int i = 0; i < 18; ++i) delete [] matrix[i];
                delete [] matrix;
            };
        }; //Elements
        
        //Assemble matrices and vectors
        ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        
        ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
        ierr = VecAssemblyEnd(b);CHKERRQ(ierr);
        
        // MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
        // ierr = VecView(b,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
        
        //Create KSP context to solve the linear system
        ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
        
        ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
        



    //     ierr = KSPSetTolerances(ksp,1.e-10,PETSC_DEFAULT,PETSC_DEFAULT,
    //                             500);CHKERRQ(ierr);
        
    //     ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
        
    //     ierr = KSPGetPC(ksp,&pc);
        
    //     ierr = PCSetType(pc,PCJACOBI);
        
    //     //ierr = KSPSetType(ksp,KSPBCGS); CHKERRQ(ierr);

    //     // ierr = KSPGMRESSetRestart(ksp, 10); CHKERRQ(ierr);
        
    //        //ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
        

    //   //   ierr = MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,NULL,&nullsp);
    // // ierr = MatSetNullSpace(A, nullsp);
    // // ierr = MatNullSpaceDestroy(&nullsp);



#if defined(PETSC_HAVE_MUMPS)
        ierr = KSPSetType(ksp,KSPPREONLY);
        ierr = KSPGetPC(ksp,&pc);
        ierr = PCSetType(pc, PCLU);
#endif          
        ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
        ierr = KSPSetUp(ksp);



        ierr = KSPSolve(ksp,b,u);CHKERRQ(ierr);

        ierr = KSPGetTotalIterations(ksp, &iterations);            

        //ierr = VecView(u,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKERRQ(ierr);
        
        //Gathers the solution vector to the master process
        ierr = VecScatterCreateToAll(u, &ctx, &All);CHKERRQ(ierr);
        ierr = VecScatterBegin(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);CHKERRQ(ierr);
        ierr = VecScatterEnd(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);CHKERRQ(ierr);
        ierr = VecScatterDestroy(&ctx);CHKERRQ(ierr);

        ierr = VecScatterCreateToAll(b, &ctx, &Allu);CHKERRQ(ierr);
        ierr = VecScatterBegin(ctx, b, Allu, INSERT_VALUES, SCATTER_FORWARD);CHKERRQ(ierr);
        ierr = VecScatterEnd(ctx, b, Allu, INSERT_VALUES, SCATTER_FORWARD);CHKERRQ(ierr);
        ierr = VecScatterDestroy(&ctx);CHKERRQ(ierr);
        
        //Updates nodal values
        double p_;
        duNorm = 0.;
        double dpNorm = 0.;
        Ione = 1;
        
        for (int i = 0; i < numNodes; ++i){
            Ii = 2*i;
            ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
            nodes_[i] -> incrementAcceleration(0,val);
            nodes_[i] -> incrementVelocity(0,val*gamma*dTime);

            ierr = VecGetValues(Allu, Ione, &Ii, &val);CHKERRQ(ierr);
            duNorm += val*val;
        
            Ii = 2*i+1;
            ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
            nodes_[i] -> incrementAcceleration(1,val);
            nodes_[i] -> incrementVelocity(1,val*gamma*dTime);

            ierr = VecGetValues(Allu, Ione, &Ii, &val);CHKERRQ(ierr);
            duNorm += val*val;
        };
        
        for (int i = 0; i<numNodes; i++){
            Ii = 2*numNodes+i;
            ierr = VecGetValues(All,Ione,&Ii,&val);CHKERRQ(ierr);
            p_ = val;
            nodes_[i] -> incrementPressure(p_);

            ierr = VecGetValues(Allu,Ione,&Ii,&val);CHKERRQ(ierr);
            dpNorm += val*val;
        };
        
        //Computes the solution vector norm
        //ierr = VecNorm(u,NORM_2,&val);CHKERRQ(ierr);

        std::clock_t t2 = std::clock();
     
        if(rank == 0){

            std::cout << "Iteration = " << inewton 
                      << " (" << iterations << ")"  
                      << "   Du Norm = " << std::scientific << sqrt(duNorm) 
                      << " " << sqrt(dpNorm)
                      << "  Time (s) = " << std::fixed
                      << 1000.*(t2-t1)/CLOCKS_PER_SEC/1000. << std::endl;
        };
                  
        ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
        ierr = VecDestroy(&b); CHKERRQ(ierr);
        ierr = VecDestroy(&u); CHKERRQ(ierr);
        ierr = VecDestroy(&All); CHKERRQ(ierr);
        ierr = VecDestroy(&Allu); CHKERRQ(ierr);
        ierr = MatDestroy(&A); CHKERRQ(ierr);

        if (sqrt(duNorm) <= tolerance) {
            break;
        };



    };//Newton-Raphson
    
    return 0;
};

#endif
