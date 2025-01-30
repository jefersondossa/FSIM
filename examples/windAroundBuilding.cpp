#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include "OpenFOAMWriter.h"
#include "DataTypes.h"

bool RunChamine();
bool RunBuildingRafaelTCC();
bool RunBuildingMestrado();

int main() {
   
    //Create the mesh using the plane xy as the ground

    // RunBuildingRafaelTCC();
    RunBuildingMestrado();
    // RunChamine();

    std::cout << "Successfully run simmulation with OpenFOAM." << std::endl;
    return 0;
}


bool RunChamine(){

    OpenFOAMWriter openFOAMWriter("../chamine.msh");

    //If it is a new simulation, clears all previous simulation files
    openFOAMWriter.ClearAllFiles();
    //You can also start the simulation from the latest time of the previous simulation
    // openFOAMWriter.StartFromPreviousResults(dt,endtime,writeinterval);

    // Input velocity
    VecDouble internalField(3);
    internalField.setZero();
    internalField[0] = 42.5;//m/s
    
    //"standard" turbulence model parameters
    double puniform = 0.0;
    double nutuniform = 0.00001;
    double kuniform = 0.00107334375;
    double omegauniform = 100.0;
    
    if (!openFOAMWriter.WriteInitialConditions(internalField,puniform,nutuniform,kuniform,omegauniform)) {
        std::cerr << "Failed to write initial conditions." << std::endl;
        return 1;
    };

    //Values for 300 Kelvin air
    double visc = 1.85e-5;//kg/m-s
    double dens = 1.1774;//kg/m^3
    if (!openFOAMWriter.WriteConstant(dens,visc)) {
        std::cerr << "Failed to write the constant folder." << std::endl;
        return 1;
    };

    //Distance from inlet and outlet boundaries to the building - must be given in meters
    double dInlet = 3.;
    double dOutlet = 30.;
    //Cell size must be given in meters - this represents the max size of elements.
    double cellSizex = 0.5;
    double cellSizey = 0.5;
    double cellSizez = 4.;
    //Refinement levels from boundary to the building - default is 2 levels
    int nRefinements = 4;
    //Proportion of inlet and outlet distances to be refined from the building [0,1]
    double refProportion = 0.2; 
    //Time variables  - must be given in seconds
    double dt = 0.01;
    double endTime = 60.;
    double writeInterval = 0.1; 
    if (!openFOAMWriter.WriteSystem(dInlet, dOutlet, cellSizex, cellSizey, cellSizez, refProportion, nRefinements, dt, endTime, writeInterval)) {
        std::cerr << "Failed to write the system folder." << std::endl;
        return 1;
    };

    //Load OpenFoam RunFunctions
    system("cd OpenFOAMRun && . $WM_PROJECT_DIR/bin/tools/RunFunctions");
    //Create building surface mesh
    system("cd OpenFOAMRun && surfaceFeatures");
    //Create volume mesh
    system("cd OpenFOAMRun && blockMesh");
    //  system("cd OpenFOAMRun && decomposePar");
    //Create mesh with building refinement
    // std::string snappy = "cd OpenFOAMRun && mpirun -np " + std::to_string(openFOAMWriter.GetNSubdomains())+ " snappyHexMesh -overwrite -parallel";
    std::string snappy = "cd OpenFOAMRun && snappyHexMesh -overwrite ";
    system(snappy.c_str());
    // system("cd OpenFOAMRun && reconstructParMesh -constant");
    system("cd OpenFOAMRun && paraFoam");

    //Decompose mesh into processors
    system("cd OpenFOAMRun && decomposePar");
    //Renumber mesh to reduce bandwidth - improves performance
    std::string bandwidth = "cd OpenFOAMRun && mpirun -np " + std::to_string(openFOAMWriter.GetNSubdomains())+ " renumberMesh -overwrite -parallel";
    system(bandwidth.c_str());
    //Run the simulation
    std::string run = "cd OpenFOAMRun && mpirun -np " + std::to_string(openFOAMWriter.GetNSubdomains())+ " foamRun -parallel";
    system(run.c_str());
    //Reconstruct the parallel mesh and results to view in ParaView
    system("cd OpenFOAMRun && reconstructPar");

}

bool RunBuildingRafaelTCC(){

    OpenFOAMWriter openFOAMWriter("../Projeto_c.msh");
    // OpenFOAMWriter openFOAMWriter("../ed1.msh");
    // OpenFOAMWriter openFOAMWriter("../ed2.msh");

    //If it is a new simulation, clears all previous simulation files
    openFOAMWriter.ClearAllFiles();
    //You can also start the simulation from the latest time of the previous simulation
    // openFOAMWriter.StartFromPreviousResults(dt,endtime,writeinterval);

    // Input velocity
    VecDouble internalField(3);
    internalField.setZero();
    internalField[0] = 30.;//m/s
    
    //"standard" turbulence model parameters
    double puniform = 0.0;
    double nutuniform = 0.00001;
    double kuniform = 0.00107334375;
    double omegauniform = 100.0;
    
    if (!openFOAMWriter.WriteInitialConditions(internalField,puniform,nutuniform,kuniform,omegauniform)) {
        std::cerr << "Failed to write initial conditions." << std::endl;
        return 1;
    };

    //Values for 300 Kelvin air
    double visc = 1.85e-5;//kg/m-s
    double dens = 1.1774;//kg/m^3
    if (!openFOAMWriter.WriteConstant(dens,visc)) {
        std::cerr << "Failed to write the constant folder." << std::endl;
        return 1;
    };

    //Distance from inlet and outlet boundaries to the building - must be given in meters
    double dInlet = 40.;
    double dOutlet = 200.;
    //Cell size must be given in meters - this represents the max size of elements.
    double cellSizex = 5.;
    double cellSizey = 5.;
    double cellSizez = 5.;
    //Refinement levels from boundary to the building - default is 2 levels
    int nRefinements = 2;
    //Proportion of inlet and outlet distances to be refined from the building [0,1]
    double refProportion = 0.2; 
    //Time variables  - must be given in seconds
    double dt = 0.01;
    double endTime = 1.;
    double writeInterval = 0.1; 
    if (!openFOAMWriter.WriteSystem(dInlet, dOutlet, cellSizex, cellSizey, cellSizez, refProportion, nRefinements, dt, endTime, writeInterval)) {
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
    std::string snappy = "cd OpenFOAMRun && snappyHexMesh -overwrite ";
    system(snappy.c_str());

    //Write initial velocity field
    system("cd OpenFOAMRun && funkySetFields -time 0");

    //Decompose mesh into processors
    system("cd OpenFOAMRun && decomposePar");

    //Renumber mesh to reduce bandwidth - improves performance
    std::string bandwidth = "cd OpenFOAMRun && mpirun -np " + std::to_string(openFOAMWriter.GetNSubdomains())+ " renumberMesh -overwrite -parallel";
    system(bandwidth.c_str());
    //Run the simulation
    std::string run = "cd OpenFOAMRun && mpirun -np " + std::to_string(openFOAMWriter.GetNSubdomains())+ " foamRun -parallel";
    system(run.c_str());
    //Reconstruct the parallel mesh and results to view in ParaView
    system("cd OpenFOAMRun && reconstructPar");
    system("cd OpenFOAMRun && paraFoam");
}