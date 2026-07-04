#include "XDMFGenerator.h"

void XDMFGenerator::ReadInitialSolution(const std::string& inputPrev, const std::string& inputCurr){

    int rank;

    // MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    // hid_t filePrevious;
    // hid_t fileCurrent;
    // hid_t dataset;
    // herr_t status;

    // REAL *vecValues;
    // REAL *scaValues;
    // vecValues = new double[3*numNodes];
    // scaValues = new double[numNodes];

    // filePrevious = H5Fopen(inputPrev.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
    // fileCurrent = H5Fopen(inputCurr.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);

    // char datasetName[] = "/velocity";
    // dataset = H5Dopen( filePrevious, datasetName, H5P_DEFAULT );
    // status = H5Dread( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vecValues[0]);
    // status = H5Dclose(dataset);
    // for (int i = 0; i < numNodes; ++i){
    //     NodeVec()[i] -> setPreviousVelocityComponent(0,vecValues[3*i  ]);
    //     NodeVec()[i] -> setPreviousVelocityComponent(1,vecValues[3*i+1]);
    // }
    // dataset = H5Dopen( fileCurrent, datasetName, H5P_DEFAULT );
    // status = H5Dread( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vecValues[0]);
    // status = H5Dclose(dataset);
    // for (int i = 0; i < numNodes; ++i){
    //     NodeVec()[i] -> setVelocityComponent(0,vecValues[3*i  ]);
    //     NodeVec()[i] -> setVelocityComponent(1,vecValues[3*i+1]);
    // }

    // char datasetName2[] = "/acceleration";
    // dataset = H5Dopen( filePrevious, datasetName2, H5P_DEFAULT );
    // status = H5Dread( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vecValues[0]);
    // status = H5Dclose(dataset);
    // for (int i = 0; i < numNodes; ++i){
    //     NodeVec()[i] -> setPreviousAccelerationComponent(0,vecValues[3*i  ]);
    //     NodeVec()[i] -> setPreviousAccelerationComponent(1,vecValues[3*i+1]);
    // }
    // dataset = H5Dopen( fileCurrent, datasetName2, H5P_DEFAULT );
    // status = H5Dread( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vecValues[0]);
    // status = H5Dclose(dataset);
    // for (int i = 0; i < numNodes; ++i){
    //     NodeVec()[i] -> setAccelerationComponent(0,vecValues[3*i  ]);
    //     NodeVec()[i] -> setAccelerationComponent(1,vecValues[3*i+1]);
    // }

    // char datasetName3[] = "/lagrangeMultiplers";
    // dataset = H5Dopen( fileCurrent, datasetName3, H5P_DEFAULT );
    // status = H5Dread( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vecValues[0]);
    // status = H5Dclose(dataset);
    // for (int i = 0; i < numNodes; ++i){
    //     NodeVec()[i] -> setLagrangeMultiplier(0,vecValues[3*i  ]);
    //     NodeVec()[i] -> setLagrangeMultiplier(1,vecValues[3*i+1]);
    // }

    // char datasetName4[] = "/pressure";
    // dataset = H5Dopen( fileCurrent, datasetName4, H5P_DEFAULT );
    // status = H5Dread( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &scaValues[0]);
    // status = H5Dclose(dataset);
    // for (int i = 0; i < numNodes; ++i){
    //     NodeVec()[i] -> setPressure(vecValues[i]);
    // }

    // delete [] vecValues;
    // delete [] scaValues;

    // //End HDF5 file
    // status = H5Fclose(filePrevious);

}