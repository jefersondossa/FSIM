#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include "OpenFOAMWriter.h"
#include "DataTypes.h"

int main() {
   
    OpenFOAMWriter openFOAMWriter("../ed1.msh");
    
    

    VecDouble internalField(3);
    internalField.setZero();
    internalField[0] = 1;
    double puniform = 0.0;
    double nutuniform = 0.00001;
    double kuniform = 0.00107334375;
    double omegauniform = 100.0;
    
    if (!openFOAMWriter.WriteInitialConditions(internalField,puniform,nutuniform,kuniform,omegauniform)) {
        std::cerr << "Failed to write initial conditions." << std::endl;
        return 1;
    };

    double visc = 0.000001;
    double dens = 1.;
    if (!openFOAMWriter.WriteConstant(dens,visc)) {
        std::cerr << "Failed to write the constant folder." << std::endl;
        return 1;
    };

    double dInlet = 60.;
    double dOutlet = 120.;
    double cellSize = 5.;
    double dt = 0.01;
    double endTime = 3.;
    double writeInterval = 0.1; 
    double refProportion = 0.2; // proportion of inlet and outlet distances to be refined from the building
    if (!openFOAMWriter.WriteSystem(dInlet, dOutlet, cellSize, refProportion, dt, endTime, writeInterval)) {
        std::cerr << "Failed to write the system folder." << std::endl;
        return 1;
    };

    //Load OpenFoam RunFunctions
    system("cd OpenFOAMRun && . $WM_PROJECT_DIR/bin/tools/RunFunctions");
    //Create building surface mesh
    system("cd OpenFOAMRun && surfaceFeatures");
    //Create volume mesh
    system("cd OpenFOAMRun && blockMesh");
    //Create mesh with building refinement
    system("cd OpenFOAMRun && snappyHexMesh -overwrite");
    //Decompose mesh into processors
    system("cd OpenFOAMRun && decomposePar");
    //Renumber mesh to reduce bandwidth - improves performance
    std::string bandwidth = "cd OpenFOAMRun && mpirun -np " + std::to_string(openFOAMWriter.GetNSubdomains())+ " renumberMesh -overwrite -parallel";
    char const *bandwidthchar = bandwidth.c_str();
    system(bandwidthchar);
    //Run the simulation
    std::string run = "cd OpenFOAMRun && mpirun -np " + std::to_string(openFOAMWriter.GetNSubdomains())+ " foamRun -parallel";
    char const *runchar = run.c_str();
    system(runchar);
    //Reconstruct the parallel mesh and results to view in ParaView
    system("cd OpenFOAMRun && reconstructPar");



    std::cout << "Successfully run simmulation with OpenFOAM." << std::endl;
    return 0;
}
