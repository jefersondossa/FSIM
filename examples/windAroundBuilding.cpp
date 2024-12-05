#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include "OpenFOAMWriter.h"
#include "DataTypes.h"

int main() {
   
    OpenFOAMWriter openFOAMWriter("../ed1.msh");
    
    double dInlet = 10.;
    double dOutlet = 30.;
    double cellSize = 1.;

    VecDouble internalField(3);
    internalField.setZero();
    internalField[0] = 0.5;
    double puniform = 0.0;
    double nutuniform = 0.00001;
    double kuniform = 0.00107334375;
    double omegauniform = 100.0;
    
    if (!openFOAMWriter.WriteInitialConditions(internalField,puniform,nutuniform,kuniform,omegauniform)) {
        std::cerr << "Failed to write initial conditions." << std::endl;
        return 1;
    };

    openFOAMWriter.WriteConstant();

    if (!openFOAMWriter.WriteBlockMeshDict(dInlet, dOutlet, cellSize)) {
        std::cerr << "Failed to write the blockMeshDict file." << std::endl;
        return 1;
    };



    std::cout << "Successfully create OPENFoamFiles." << std::endl;
    return 0;
}
