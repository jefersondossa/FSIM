#include "CompMeshTools.h"
#include "fstream"

#ifdef HAS_METIS
#include <metis.h>
#endif

void CompMeshTools::DomainDecompositionMETIS(CompMesh *cmesh){

    PanicButton();//Need refactor
    // std::string mirror2;
    // mirror2 = "domain_decomposition.txt";
    // std::ofstream mirrorData(mirror2.c_str());
    
    // int size;

    // MPI_Comm_size(PETSC_COMM_WORLD, &size);

    // idx_t objval;
    // idx_t numEl = cmesh->NElements();
    // idx_t numNd = cmesh->NNodes();
    // idx_t ssize = size;
    // idx_t one = 1;
    // int nElNodes = cmesh->NElNodes();
    // idx_t elem_start[numEl+1], elem_connec[nElNodes*numEl];
    // cmesh->part_elem = new int[numEl];
    // cmesh->part_nodes = new int[numNd];


    // for (int i = 0; i < numEl+1; i++){
    //     elem_start[i]=nElNodes*i;
    // };
    // for (int jel = 0; jel < numEl; jel++){
    //     auto connec=cmesh->ElementVec()[jel]->getConnectivity();        
        
    //     for (int i=0; i<nElNodes; i++){
    //     elem_connec[nElNodes*jel+i] = connec[i];
    //     };
    // };

    // //Performs the domain decomposition
    // if (size == 1){
    //     for (int i = 0; i < numNd; i++) cmesh->part_nodes[i] = 0;
    //     for (int i = 0; i < numEl; i++) cmesh->part_elem[i] = 0;
    // } else {
    //     METIS_PartMeshDual(&numEl, &numNd, elem_start, elem_connec, \
    //                             NULL, NULL, &one, &ssize, NULL, NULL,    \
    //                             &objval, cmesh->part_elem, cmesh->part_nodes);
    // }
    
    // mirrorData << std::endl \
    //            << "FLUID MESH DOMAIN DECOMPOSITION - ELEMENTS" << std::endl;
    // for(int i = 0; i < numEl; i++){
    //     mirrorData << "process = " << cmesh->part_elem[i] \
    //                << ", element = " << i << std::endl;
    // };

    // mirrorData << std::endl \
    //            << "FLUID MESH DOMAIN DECOMPOSITION - NODES" << std::endl;
    // for(int i = 0; i < numNd; i++){
    //     mirrorData << "process = " << cmesh->part_nodes[i] \
    //                << ", node = " << i << std::endl;
    // };


    
}


bool CompMeshTools::searchNodeCorrespondence(VecDouble &x, CompMesh *cmesh, int64_t &elCorr, VecDouble &xsiCorr, Element *elemsearch){
    
    int DIM = cmesh->Dimension();
    int DEG = cmesh->GetDefaultOrder();
    VecDouble xsiCC(DIM);
    std::pair<VecDouble,VecDouble> XK;

    elCorr = 150000;
    VecDouble xsi(DIM);
    VecDouble x_(DIM);
    VecDouble deltaX(DIM);
    VecDouble deltaXsi(DIM);
    xsi.setZero(); x_.setZero(); deltaX.setZero(); deltaXsi.setZero();
    
    xsiCC.fill(1.e10);
    xsiCorr.fill(1.e50);
    xsi.fill(0.);
    x_.fill(0.);
    
    // Element *elemsearch = nullptr;
    // if (cmesh->ElementVec()[elSearch]->Dimension() != cmesh->Dimension()){
    //     for (int i=0; i<cmesh->NElements(); i++){
    //        if(cmesh->ElementVec()[i]->Dimension() != cmesh->Dimension()) continue;
    //        elemsearch = cmesh->ElementVec()[i];
    //        break;
    //    }
    // } else {
    //     elemsearch = cmesh->ElementVec()[elSearch];
    // }
    VecInt connec = elemsearch -> Reference() -> getGeometricNodes();
    
    auto &integdata = elemsearch->IntegrationData();
    integdata.fAdimCoord = xsi;
    elemsearch-> Reference() -> ComputeJacobianSearch(integdata);
    int nElNodes = integdata.fPhi.size();

    for (int i = 0; i < nElNodes; i++){
        VecDouble xint = cmesh->Reference()->NodeVec()[connec[i]] -> getCoordinates();
        for (int k = 0; k < DIM; k++){
            x_[k] += xint[k] * integdata.fPhi[i];
        }        
    };

    double error = 1.e6;
    int iterations = 0;

    while ((error > 1.e-8) && (iterations < 4)) {
        
        iterations++;
        
        for (int k = 0; k < DIM; k++) deltaX[k] = x[k] - x_[k];
        deltaXsi.setZero();
        
        elemsearch -> Reference() -> ComputeJacobianSearch(elemsearch->IntegrationData());

        deltaXsi = elemsearch->IntegrationData().fA0Inv*deltaX;    

        xsi += deltaXsi;
        x_.setZero();
        
        integdata.fAdimCoord = xsi;
        
        for (int i=0; i<nElNodes; i++){
            elemsearch -> Reference() -> ComputeJacobianSearch(elemsearch->IntegrationData());
            VecDouble xint = cmesh->Reference()->NodeVec()[connec[i]] -> getCoordinates();
            for (int k = 0; k < DIM; k++)x_[k] += xint[k] * integdata.fPhi[i];
        };

        error = std::sqrt(deltaXsi[0]*deltaXsi[0] + deltaXsi[1]*deltaXsi[1]);
    };
    
    double t1 = -1.e-2;
    double t2 =  1. - t1;
    
    xsiCC[0] = xsi[0];
    xsiCC[1] = xsi[1];


    switch (elemsearch->Type())
    {
    case ETriangle:
        if ((xsiCC[0] >= t1) && (xsiCC[1] >= t1) && ((1. - xsiCC[0] - xsiCC[1]) >= t1) &&
            (xsiCC[0] <= t2) && (xsiCC[1] <= t2) && ((1. - xsiCC[0] - xsiCC[1]) <= t2)){

            xsiCorr[0] = xsi[0]; xsiCorr[1] = xsi[1];
            elCorr = elemsearch->Index();
            return true;
        }
        break;
    case EQuadrilateral:
        if ((xsiCC[0] >= t1-1.) && (xsiCC[1] >= t1-1.) &&
            (xsiCC[0] <= t2) && (xsiCC[1] <= t2)){

            xsiCorr[0] = xsi[0]; xsiCorr[1] = xsi[1];
            elCorr = elemsearch->Index();
            return true;
        }
        break;
    case EOneD:
        if ((xsiCC[0] >= t1-1.) && (xsiCC[0] <= 1-t1)){
            if (fabs(xsiCC[1]) >= fabs(t1)) return false;
            xsiCorr[0] = xsi[0];
            elCorr = elemsearch->Index();
            return true;
        }
        break;

    default:
        PanicButton();
        return false;
        break;
    }

    if (fabs(xsi[0]) > 2.) {
        //std::cout << "PROBLEM SEARCHING NODE CORRESPONDENCE " << std::endl;
        //PanicButton();
        return false;
    }
    if (elCorr == 150000){
        //PanicButton();
        return false;
    } 
};