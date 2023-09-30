#include "FSInteraction.h"

//------------------------------------------------------------------------------
//--------------------------------IMPLEMENTATION--------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//-------------------COMPUTE NODAL CORRESPONDECE WITH ELEMENTS------------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void FSInteraction<DIM,DEG>::searchSolidNodeCorrespondence(int interface, int iSol){
    
    for (int isolid = 0; isolid < numNodesSolid; isolid++){

        VecInt connec;
        ShapeFunction shapeQuad(DIM,DEG);
        int nElNodes = fluidModel->NElNodes();
        VecDouble phi_(nElNodes);
        
        MatrixDouble ainv(DIM,DIM);

        double xsiCC[3];
        // std::pair<double*,double*> XK;
        int elemC;

        VecDouble x = nodesSolid_[iSol][isolid] -> getCoordinates();
        
        elemC = 150000;
        VecDouble xsiC(DIM);
        for (int k = 0; k<DIM; k++) xsiC[k] = 1.e50;
        VecDouble xsi(DIM);
        VecDouble x_(DIM);
        VecDouble deltaX(DIM);
        VecDouble deltaXsi(DIM);

        nodesSolid_[iSol][isolid] -> setNodalCorrespondence(elemC,xsiC);
        
        for (int ibound = 0; ibound < numElemFluidBoundary; ibound++){
            
            if (boundaryFluid_[ibound] -> getBoundaryGroup() == interface){
                
                int jel = boundaryFluid_[ibound] -> getElement();
                
                connec = elementsFluid_[jel] -> getConnectivity();
                
                //get boxes information        
                // XK = elementsFluid_[jel] -> getXIntersectionParameter();
                
                //Chech if the node is inside the element box
                // if ((x(0) < XK.first(0)) || (x(0) > XK.second(0)) ||
                //     (x(1) < XK.first(1)) || (x(1) > XK.second(1))) continue;
                
                //Compute nodal correspondence
                xsiCC[0] = 1.e10;
                xsiCC[1] = 1.e10;
                xsiCC[2] = 1.e10;
                
                x_.setZero();
                xsi.fill(1./3.);
                
                shapeQuad.evaluate(xsi,phi_);
                
                for (int i = 0; i < nElNodes; i++){
                    VecDouble xint = nodesFluid_[connec[i]] -> getCoordinates();
                    x_ += xint * phi_[i];
                };
                
                double error = 1.e6;                
                int iterations = 0;
                
                while ((error > 1.e-8) && (iterations < 4)) {
                    
                    iterations++;
                
                    deltaX = x - x_;                    
                    deltaXsi.setZero();
                    
                    double djac_ = 0.;
                    elementsFluid_[jel] -> getJacobianMatrix(xsi,ainv,djac_,0);
                    
                    // noalias(deltaXsi) = prod(trans(ainv),deltaX);

                    for (int i = 0; i < DIM; i++)
                        for (int j = 0; j < DIM; j++)
                            deltaXsi[i] += ainv(j,i) * deltaX[j];

                    xsi += deltaXsi;
                    x_.setZero();
                    
                    shapeQuad.evaluate(xsi,phi_);
                    
                    for (int i=0; i<nElNodes; i++){
                        VecDouble xint = nodesFluid_[connec[i]] -> getCoordinates();
                        x_ += xint * phi_[i];  
                    }; 

                    error = std::sqrt(deltaXsi[0]*deltaXsi[0] + deltaXsi[1]*deltaXsi[1]);
                };
                
                double t1 = -1.e-1;
                double t2 =  1. - t1;
                
                xsiCC[0] = xsi[0];
                xsiCC[1] = xsi[1];       
                xsiCC[2] = 1. - xsiCC[0] - xsiCC[1];
                
                if ((xsiCC[0] >= t1) && (xsiCC[1] >= t1) && (xsiCC[2] >= t1) &&
                    (xsiCC[0] <= t2) && (xsiCC[1] <= t2) && (xsiCC[2] <= t2)){
                    
                    xsiC[0] = xsi[0]; xsiC[1] = xsi[1];
                    elemC = jel;
                    nodesSolid_[iSol][isolid] -> setNodalCorrespondence(elemC,xsiC);
                    break;
                };           
            };
        };

        // nodesSolid_[iSol][isolid] -> setNodalCorrespondence(elemC,xsiC);
        
        // std::cout << "isolid " << isolid << " " << interface << " " << elemC << " " << x_[0] << " " << x_[1] << " " << xsiC[0] << " " << xsiC[1] << std::endl;


    };

    return;

};

//------------------------------------------------------------------------------
//-------------------COMPUTE NODAL CORRESPONDECE WITH ELEMENTS------------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void FSInteraction<DIM,DEG>::searchSolidNodeCorrespondenceArlequin(int interface, int iSol){
    
    for (int isolid = 0; isolid < numNodesSolid; isolid++){
        int nElNodes = arlequinModel->MeshVec()[1]->NElNodes();
        VecInt connec;
        ShapeFunction shapeQuad(DIM,DEG);
        VecDouble phi_(nElNodes);
        
        MatrixDouble ainv(DIM,DIM);

        double xsiCC[3];
        // std::pair<typename Elements::DimVector,typename Elements::DimVector> XK;
        int elemC;

        VecDouble x = nodesSolid_[iSol][isolid] -> getCoordinates();
        
        elemC = 150000;
        VecDouble xsiC(DIM);
        for (int k = 0; k<DIM; k++) xsiC[k] = 1.e50;
        VecDouble xsi(DIM);
        double x_[DIM];
        double deltaX[DIM];
        double deltaXsi[DIM];
        nodesSolid_[iSol][isolid] -> setNodalCorrespondence(elemC,xsiC);
        
        for (int ibound = 0; ibound < numElemArlequinBoundaryFine; ibound++){
            
            if (boundaryArlequinFine_[ibound] -> getBoundaryGroup() == interface){
                
                int jel = boundaryArlequinFine_[ibound] -> getElement();
                
                connec = elementsArlequinFine_[jel] -> getConnectivity();
                
                //get boxes information        
                // XK = elementsArlequinFine_[jel] -> getXIntersectionParameter();
                
                //Chech if the node is inside the element box
                // if ((x(0) < XK.first(0)) || (x(0) > XK.second(0)) ||
                //     (x(1) < XK.first(1)) || (x(1) > XK.second(1))) continue;
                
                //Compute nodal correspondence
                xsiCC[0] = 1.e10;
                xsiCC[1] = 1.e10;
                xsiCC[2] = 1.e10;
                
                for (int k = 0; k<DIM; k++){
                    xsi[k] = 1. / 3.;
                    x_[k] = 0.;
                }
                
                shapeQuad.evaluate(xsi,phi_);
                                
                for (int i = 0; i < nElNodes; i++){
                    VecDouble xint = nodesArlequinFine_[connec[i]] -> getCoordinates();
                    for (int k = 0; k<DIM; k++) x_[k] += xint[k] * phi_[i];                
                };
                
                double error = 1.e6;
                
                int iterations = 0;
                
                while ((error > 1.e-8) && (iterations < 4)) {
                    
                    iterations++;
                    
                    for (int k = 0; k<DIM; k++){
                        deltaX[k] = x[k] - x_[k];
                        deltaXsi[k] = 0.;
                    }
                    
                    double djac_ = 0.;
                    elementsArlequinFine_[jel] -> getJacobianMatrix(xsi,ainv,djac_,0);
                    
                    // noalias(deltaXsi) = prod(trans(ainv),deltaX);

                    for (int i = 0; i < DIM; i++)
                        for (int j = 0; j < DIM; j++)
                            deltaXsi[i] += ainv(j,i) * deltaX[j];
                    
                    for (int k = 0; k<DIM; k++){
                        xsi[k] += deltaXsi[k];
                        x_[k] = 0.;
                    }
                    
                    shapeQuad.evaluate(xsi,phi_);
                    
                    for (int i=0; i<nElNodes; i++){
                        VecDouble xint = nodesArlequinFine_[connec[i]] -> getCoordinates();
                        for (int k = 0; k<DIM; k++) x_[k] += xint[k] * phi_[i];
                    };                   

                    error = std::sqrt(deltaXsi[0]*deltaXsi[0] + deltaXsi[1]*deltaXsi[1]);
                };
                
                double t1 = -1.e-1;
                double t2 =  1. - t1;
                
                xsiCC[0] = xsi[0];
                xsiCC[1] = xsi[1];       
                xsiCC[2] = 1. - xsiCC[0] - xsiCC[1];
                
                if ((xsiCC[0] >= t1) && (xsiCC[1] >= t1) && (xsiCC[2] >= t1) &&
                    (xsiCC[0] <= t2) && (xsiCC[1] <= t2) && (xsiCC[2] <= t2)){
                    
                    xsiC[0] = xsi[0]; xsiC[1] = xsi[1];
                    elemC = jel;
                    nodesSolid_[iSol][isolid] -> setNodalCorrespondence(elemC,xsiC);
                    break;
                };           
            };
        };

        // if (rank == 0) std::cout << "isolid " << isolid << " " << interface << " " << elemC << " " << x_[0] << " " << x_[1] << " " << xsiC[0] << " " << xsiC[1] << std::endl;
    };
};

//------------------------------------------------------------------------------
//-------------------COMPUTE NODAL CORRESPONDECE WITH ELEMENTS------------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void FSInteraction<DIM,DEG>::searchFluidNodeCorrespondence(int interface){
    
    for (int ibound = 0; ibound < numElemFluidBoundary; ibound++){

        if (boundaryFluid_[ibound] -> getBoundaryGroup() == interface){    
            int elemC;
            double xsiC;

            VecInt connec = boundaryFluid_[ibound] -> getBoundaryConnectivity();
            int nBdNodes = fluidModel->NBdNodes();
            for (int inode = 0; inode < nBdNodes; inode++){
                VecDouble x = nodesFluid_[connec[inode]] -> getCoordinates();
                
                searchcorrespondencefluid_(&x[0], &x[1], &elemC, &xsiC);
                VecDouble xsi(DIM);
                xsi[0] = xsiC;

                // std::cout << "asdasd " << elemC << " " << xsiC << std::endl; 

                nodesFluid_[connec[inode]] -> setNodalCorrespondence(elemC,xsi);
                
                // std::cout << "isolid " << connec(inode) << " " << elemC << " " << xsi(0) << std::endl;
            };
        };        
    };
};

//------------------------------------------------------------------------------
//-------------------COMPUTE NODAL CORRESPONDECE WITH ELEMENTS------------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void FSInteraction<DIM,DEG>::searchArlequinNodeCorrespondence(int interface){
    
    int flag = 0;

    for (int ibound = 0; ibound < numElemArlequinBoundaryFine; ibound++){

        if (boundaryArlequinFine_[ibound] -> getBoundaryGroup() == interface){
            flag++;
            int elemC;
            double xsiC;

            VecInt connec = boundaryArlequinFine_[ibound] -> getBoundaryConnectivity();
            int nBdNodes = arlequinModel->MeshVec()[1]->NBdNodes();
            for (int inode = 0; inode < nBdNodes; inode++){
                VecDouble x = nodesArlequinFine_[connec[inode]] -> getCoordinates();
                
                searchcorrespondencefluid_(&x[0], &x[1], &elemC, &xsiC);
                VecDouble xsi(DIM);
                xsi[0] = xsiC;

                // if (rank == 0) std::cout << "asdasd " << elemC << " " << xsiC << " " << " " << connec[inode] << " " << interface << " " << numElemArlequinFine << " " << flag << std::endl; 

                nodesArlequinFine_[connec[inode]] -> setNodalCorrespondence(elemC,xsi);
                
                // std::cout << "isolid " << connec(inode) << " " << elemC << " " << xsi(0) << std::endl;
            };
        };        
    };
};

//------------------------------------------------------------------------------
//-------------------------COMPUTE ELEMENT REGIONS/BOXES------------------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void FSInteraction<DIM,DEG>::setElementBoxes() {
    
    VecInt connec;
    VecDouble xk(DIM), Xk(DIM);
    double dCk[3], dck[3];

    //Compute element boxes for coarse model
    //Only function for straight elements
    for (int jel = 0; jel < numElemFluid; jel++){
        connec = elementsFluid_[jel] -> getConnectivity();
        VecDouble x1 = nodesFluid_[connec[0]] -> getCoordinates();
        VecDouble x2 = nodesFluid_[connec[1]] -> getCoordinates();
        VecDouble x3 = nodesFluid_[connec[2]] -> getCoordinates();      

        xk[0] = std::min(x1[0],std::min(x2[0], x3[0]));
        xk[1] = std::min(x1[1],std::min(x2[1], x3[1]));

        Xk[0] = std::max(x1[0],std::max(x2[0], x3[0]));
        Xk[1] = std::max(x1[1],std::max(x2[1], x3[1]));        
        
        elementsFluid_[jel] -> setIntersectionParameters(xk, Xk);
    };

};

//------------------------------------------------------------------------------
//---------SETS FLUID AND SOLID MODELS AND GETS ITS BASIC INFORMATIONS----------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void FSInteraction<DIM,DEG>::preProcessFluid(){
    
    
    numInterfaces = fluidModel->getNumberofFSIInterfaces();
    
    if(rank==0) std::cout << "Number of interfaces " 
                          << numInterfaces << std::endl;

    groupInterfaces.reserve(numInterfaces);
    nodesSolid_.reserve(numInterfaces);
    //    int groupInterfaces2[numInterfaces];

    int fl=0;
    //Sets fluid elements and sides on interface boundaries
    for (int i=0; i<numElemFluidBoundary; i++){
        if ((boundaryFluid_[i] -> getConstrain(0) == 3) ||
            (boundaryFluid_[i] -> getConstrain(1) == 3)) {

            VecInt connectB = boundaryFluid_[i] -> getBoundaryConnectivity();

            for (int j=0; j<numElemFluid; j++){
                VecInt connect = elementsFluid_[j] -> getConnectivity();
                int nBdNodes = fluidModel->NBdNodes();
                int flag = 0;
                int side[nBdNodes];
                int nElNodes = fluidModel->NElNodes();
                for (int k=0; k<nElNodes; k++){
                    if ((connectB[0] == connect[k]) || 
                        (connectB[1] == connect[k]) ||
                        (connectB[2] == connect[k])){
                        side[flag] = k;
                        flag++;
                    };
                };
                
                if (flag == nBdNodes){
                    boundaryFluid_[i] -> setElement(j);
                    elementsFluid_[j] -> setFSIInterface();
                    // Counts number of interfaces
                    if (fl == 0){
                        int aux = boundaryFluid_[i] -> getBoundaryGroup();
                        groupInterfaces.push_back(aux);
                        fl++;
                    }else{
                        int fl2 = 0;
                        for (int m=0; m<fl; m++){
                            if (boundaryFluid_[i] -> getBoundaryGroup() == groupInterfaces[m]) fl2++;
                        };
                        if (fl2 == 0){
                            int aux = boundaryFluid_[i] -> getBoundaryGroup();
                            groupInterfaces.push_back(aux);
                            fl++;
                        };
                    };

                    //Sets element index and side
                    if ((side[0]==4) || (side[1]==4) || (side[2]==4)){
                        boundaryFluid_[i] -> setElementSide(0);
                        elementsFluid_[boundaryFluid_[i]->getElement()] -> 
                            setElemSideInBoundary(0);
                    };
                    if ((side[0]==5) || (side[1]==5) || (side[2]==5)){
                        boundaryFluid_[i] -> setElementSide(1);
                        elementsFluid_[boundaryFluid_[i]->getElement()] -> 
                            setElemSideInBoundary(1);
                    };
                    if ((side[0]==3) || (side[1]==3) || (side[2]==3)){
                        boundaryFluid_[i] -> setElementSide(2);
                        elementsFluid_[boundaryFluid_[i]->getElement()] -> 
                            setElemSideInBoundary(2);
                    };
                };
            };
        };
    };   
    
    //Get Solid nodal Positions
    // for (int k=0; k<numInterfaces; k++){

    // };

    for (int k=0; k<numInterfaces; k++){    
        std::vector<Node *> j;

        nodesSolid_.push_back(j);
    
        nodesSolid_[k].reserve(numNodesSolid);
        int index = 0;
        
        for (int i=0; i<numNodesSolid; i++){
            VecDouble x(DIM);
            int inode = i+1;
            getsolidposition_(&inode,&x[0],&x[1]);
            
            Node *node = new Node(x,index++);
            nodesSolid_[k].push_back(node);
        };
    };

    setElementBoxes();
    
    // Search solid node correspondence into fluid elements
    for (int i = 0; i < numInterfaces; i++){
        int interf = groupInterfaces[i];
        searchSolidNodeCorrespondence(interf,i);
    };

    // Search fluid node correspondences into solid elements
    for (int i = 0; i < numInterfaces; i++){
        int interf = groupInterfaces[i];       
        searchFluidNodeCorrespondence(interf);
    };
    
    // Set element mesh moving parameters
    double vMax = 0., vMin = 1.e10;
    for (int i = 0; i < numElemFluid; i++){
        double v = elementsFluid_[i] -> getJacobian();
        if (v > vMax) vMax = v;
        if (v < vMin) vMin = v;
    };
    for (int i = 0; i < numElemFluid; i++){
        double v = elementsFluid_[i] -> getJacobian();
        double eta = 1 + (1. - vMin / vMax) / (v / vMax);
        elementsFluid_[i] -> setMeshMovingParameter(eta);

        // std::cout << "MESH MOVING PARAMETER " << i << " " << eta << " " << vMin << " " << vMax << std::endl;
    };
       
    domDecompFluid = fluidModel->getDomainDecomposition();

};

//------------------------------------------------------------------------------
//---------SETS FLUID AND SOLID MODELS AND GETS ITS BASIC INFORMATIONS----------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void FSInteraction<DIM,DEG>::preProcessArlequin(){
        
    // numInterfaces = arlequinModel->MeshVec()[1]->getNumberofFSIInterfaces();
    PanicButton();
    std::cout << "Number of interfaces " << numInterfaces << std::endl;

    groupInterfaces.reserve(numInterfaces);
    nodesSolid_.reserve(numInterfaces);

    int fl=0;
    //Sets fluid elements and sides on interface boundaries
    for (int i=0; i<numElemArlequinBoundaryFine; i++){
        if ((boundaryArlequinFine_[i] -> getConstrain(0) == 3)) {
            // Counts number of interfaces
            if (fl == 0){
                int aux = boundaryArlequinFine_[i] ->getBoundaryGroup();
                groupInterfaces.push_back(aux);
                fl++;
            }else{
                int fl2 = 0;
                for (int m=0; m<fl; m++){
                    if (boundaryArlequinFine_[i] -> 
                        getBoundaryGroup() == groupInterfaces[m]) fl2++;
                };
                if (fl2 == 0){
                    int aux = boundaryArlequinFine_[i] -> 
                        getBoundaryGroup();
                    groupInterfaces.push_back(aux);
                    fl++;
                };
            };
        };
    };   
    
    //Get Solid nodal Positions
    for (int k=0; k<numInterfaces; k++){    
        std::vector<Node *> j;

        nodesSolid_.push_back(j);
    
        nodesSolid_[k].reserve(numNodesSolid);
        int index = 0;
        
        for (int i=0; i<numNodesSolid; i++){
            VecDouble x(DIM);
            int inode = i+1;
            getsolidposition_(&inode,&x[0],&x[1]);
            
            Node *node = new Node(x,index++);
            nodesSolid_[k].push_back(node);
        };
    };

// Search solid node correspondence into fluid elements
    for (int i = 0; i < numInterfaces; i++){
        int interf = groupInterfaces[i];
        searchSolidNodeCorrespondenceArlequin(interf,i);
    };

    // Search fluid node correspondences into solid elements
    for (int i = 0; i < numInterfaces; i++){
        int interf = groupInterfaces[i];       
        searchArlequinNodeCorrespondence(interf);
    };

    // Set element mesh moving parameters
    double vMax = 0., vMin = 1.e10;
    for (int i = 0; i < numElemArlequinFine; i++){
        double v = elementsArlequinFine_[i] -> getJacobian();
        if (v > vMax) vMax = v;
        if (v < vMin) vMin = v;
    };

    for (int i = 0; i < numElemArlequinFine; i++){
        double v = elementsArlequinFine_[i] -> getJacobian();
        double eta = 1 + (1. - vMin / vMax) / (v / vMax);
        elementsArlequinFine_[i] -> setMeshMovingParameter(eta);

        // std::cout << "MESH MOVING PARAMETER " << i << " " << eta << " " << vMin << " " << vMax << std::endl;
    };

    PanicButton();
    // domDecompArlequinCoarse =arlequinModel->MeshVec()[0]->getDomainDecomposition();
    // domDecompArlequinFine = arlequinModel->MeshVec()[1]->getDomainDecomposition();

};


//------------------------------------------------------------------------------
//---------SETS FLUID AND SOLID MODELS AND GETS ITS BASIC INFORMATIONS----------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void FSInteraction<DIM,DEG>::setFluidAndSolidModels(Fluid &fluid, char *in_solid){
    
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &size);

    fluidModel = &fluid;

    numElemFluid = fluidModel->NElements();
    numNodesFluid = fluidModel->NNodes();
    numElemFluidBoundary = fluidModel->NBoundElements();

    nodesFluid_  = fluidModel->NodeVec();
    elementsFluid_ = fluidModel->ElementVec();
    boundaryFluid_ = fluidModel->BoundaryVec();

    // Reads Solid input file
    preprocessing_(in_solid);
    
    getnumberofnodessolid_(&numNodesSolid);
    getnumberofelementssolid_(&numElemSolid);

    dTime = fluidModel->getProblemParameters().GetTimeStep();

    if(rank == 0) std::cout << "NumElemSolid " << numElemSolid << std::endl;

    //Pre Processing data
    preProcessFluid();

};

//------------------------------------------------------------------------------
//---------SETS FLUID AND SOLID MODELS AND GETS ITS BASIC INFORMATIONS----------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void FSInteraction<DIM,DEG>::setArlequinAndSolidModels(ArlequinModel& arlq, char *in_solid){
    
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &size);

    arlequinModel = &arlq;

    numElemArlequinCoarse = arlequinModel->MeshVec()[0]->NElements();
    numElemArlequinFine = arlequinModel->MeshVec()[1]->NElements();
    numNodesArlequinCoarse = arlequinModel->MeshVec()[0]->NNodes();
    numNodesArlequinFine = arlequinModel->MeshVec()[1]->NNodes();
    numElemArlequinBoundaryCoarse = arlequinModel->MeshVec()[0]->NBoundElements();
    numElemArlequinBoundaryFine = arlequinModel->MeshVec()[1]->NBoundElements();
    
    nodesArlequinCoarse_ = arlequinModel->MeshVec()[0]->NodeVec();
    nodesArlequinFine_ = arlequinModel->MeshVec()[1]->NodeVec();
    elementsArlequinCoarse_ = arlequinModel->MeshVec()[0]->ElementVec();
    elementsArlequinFine_ = arlequinModel->MeshVec()[1]->ElementVec();
    boundaryArlequinCoarse_ = arlequinModel->MeshVec()[0]->BoundaryVec();
    boundaryArlequinFine_ = arlequinModel->MeshVec()[1]->BoundaryVec();

    // Reads Solid input file
    preprocessing_(in_solid);

    getnumberofnodessolid_(&numNodesSolid);
    getnumberofelementssolid_(&numElemSolid);

    dTime = arlequinModel->MeshVec()[1]->getProblemParameters().GetTimeStep();

    if(rank == 0) std::cout << "NumElemSolid " << numElemSolid << std::endl;

    //Pre Processing data
    preProcessArlequin();

};

//------------------------------------------------------------------------------
//---------SETS FLUID AND SOLID MODELS AND GETS ITS BASIC INFORMATIONS----------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void FSInteraction<DIM,DEG>::updateFluidMesh(){

    for (int i = 0; i < numInterfaces; i++){

        int interf = groupInterfaces[i];       

        for (int ibound = 0; ibound < numElemFluidBoundary; ibound++){          
            if (boundaryFluid_[ibound] -> getBoundaryGroup() == interf){
                
                VecInt connec = boundaryFluid_[ibound] -> getBoundaryConnectivity();
                int nBdNodes = fluidModel->NBdNodes();

                for (int k = 0; k < nBdNodes; k++){

                    VecDouble x(DIM);
                    
                    int elem = nodesFluid_[connec[k]] -> getNodalElemCorrespondence();
                    VecDouble xsi = nodesFluid_[connec[k]] -> getNodalXsiCorrespondence();
                    
                    
                    if (rank == 0) getupdatedcoordinates_(&x[0],&x[1],&elem,&xsi[0]);
              
                    MPI_Bcast(&x[0],1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
                    MPI_Bcast(&x[1],1,MPI_DOUBLE,0,PETSC_COMM_WORLD);

                    nodesFluid_[connec[k]] -> setUpdatedCoordinates(x);
                    nodesFluid_[connec[k]] -> setCoordinates(x);

                    //if(connec(k) == 29)std::cout << "Updated Coord " << x(0) << " " << x(1) << std::endl;
                };
                
            };//if interface
        };//ibound
    };//i


    //Solves Laplace Smoothing mesh moving scheme
    fluidModel->solveSteadyLaplaceProblem(5,1.e-6);
   
    for (int i = 0; i < numNodesFluid; i++){
        VecDouble u(DIM), up(DIM);
            
        VecDouble x = nodesFluid_[i] -> getCoordinates();
        VecDouble xp = nodesFluid_[i] -> getPreviousCoordinates();
        up[0] = nodesFluid_[i] -> getPreviousMeshVelocity(0);
        up[1] = nodesFluid_[i] -> getPreviousMeshVelocity(1);
        
        u[0] = (x[0] - xp[0]) / dTime;//2. * (xp(0) - x(0)) / dTime - up(0);
        u[1] = (x[1] - xp[1]) / dTime;//2. * (xp(1) - x(1)) / dTime - up(1);
        
        nodesFluid_[i] -> setMeshVelocity(u);
    };


};

//------------------------------------------------------------------------------
//---------SETS FLUID AND SOLID MODELS AND GETS ITS BASIC INFORMATIONS----------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void FSInteraction<DIM,DEG>::updateArlequinMesh(){

    // std::cout << "AQUI6.1 " << rank << std::endl;
    // MPI_Barrier(PETSC_COMM_WORLD);

    for (int i = 0; i < numInterfaces; i++){

        int interf = groupInterfaces[i];       

        for (int ibound = 0; ibound < numElemArlequinBoundaryFine; ibound++){
            if (boundaryArlequinFine_[ibound] -> getBoundaryGroup() == interf){

                VecInt connec = boundaryArlequinFine_[ibound] -> getBoundaryConnectivity();
                int nBdNodes = arlequinModel->MeshVec()[1]->NBdNodes();
                for (int k = 0; k < nBdNodes; k++){

                    VecDouble x(2);
                    
                    int elem = nodesArlequinFine_[connec[k]] -> getNodalElemCorrespondence();
                    VecDouble xsi = nodesArlequinFine_[connec[k]] -> getNodalXsiCorrespondence();

                    // if (rank == 0) std::cout << "AQUI6.1.1 " << rank << " " << elem << " " << xsi << " " << x(0) << " " << x(1) << std::endl;
                    
                    if (rank == 0) getupdatedcoordinates_(&x[0],&x[1],&elem,&xsi[0]);

                    MPI_Bcast(&x[0],1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
                    MPI_Bcast(&x[1],1,MPI_DOUBLE,0,PETSC_COMM_WORLD);

                    // if (rank == 0) std::cout << "AQUI6.1.2 " << rank << " " << elem << " " << xsi << " " << x(0) << " " << x(1) << std::endl;
                    MPI_Barrier(PETSC_COMM_WORLD);

                    nodesArlequinFine_[connec[k]] -> setUpdatedCoordinates(x);
                    nodesArlequinFine_[connec[k]] -> setCoordinates(x);

                    // std::cout << "Updated Coord " << x(0) << " " << x(1) << std::endl;
                };
                
            };//if interface
        };//ibound
    };//i

    // std::cout << "AQUI6.2 " << rank << std::endl;
    // MPI_Barrier(PETSC_COMM_WORLD);


    //Solves Laplace Smoothing mesh moving scheme
    PanicButton();
    // arlequinModel->MeshVec()[1]->solveSteadyLaplaceProblem(1,1.e-4);
   
    // std::cout << "AQUI6.3 " << rank << std::endl;
    // MPI_Barrier(PETSC_COMM_WORLD);

    for (int i = 0; i < numNodesArlequinFine; i++){
        VecDouble up(2);
        VecDouble u(2);
            
        VecDouble x = nodesArlequinFine_[i] -> getCoordinates();
        VecDouble xp = nodesArlequinFine_[i] -> getPreviousCoordinates();
        up[0] = nodesArlequinFine_[i] -> getPreviousMeshVelocity(0);
        up[1] = nodesArlequinFine_[i] -> getPreviousMeshVelocity(1);
        
        u[0] = (x[0] - xp[0]) / dTime;//2. * (xp(0) - x(0)) / dTime - up(0);
        u[1] = (x[1] - xp[1]) / dTime;//2. * (xp(1) - x(1)) / dTime - up(1);
        
        nodesArlequinFine_[i] -> setMeshVelocity(u);
    };


};

//------------------------------------------------------------------------------
//---------SETS FLUID AND SOLID MODELS AND GETS ITS BASIC INFORMATIONS----------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void FSInteraction<DIM,DEG>::transferSolidVelocity(){
    

    for (int i = 0; i < numInterfaces; i++){

        int interf = groupInterfaces[i];       

        for (int ibound = 0; ibound < numElemFluidBoundary; ibound++){          
            if (boundaryFluid_[ibound] -> getBoundaryGroup() == interf){
                
                VecInt connec = boundaryFluid_[ibound] -> getBoundaryConnectivity();
              
                VecDouble u(DIM);
                int nBdNodes = fluidModel->NBdNodes();
                for (int k = 0; k < nBdNodes; k++){
                    
                    int elem = nodesFluid_[connec[k]] -> getNodalElemCorrespondence();
                    VecDouble xsi = nodesFluid_[connec[k]] -> getNodalXsiCorrespondence();
                    
                    
                    if (rank == 0) 
                        getinterpolatedvelocity_(&u[0],&u[1],&elem,&xsi[0]);
                        
                    MPI_Bcast(&u[0],1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
                    MPI_Bcast(&u[1],1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
               
                    nodesFluid_[connec[k]] -> setVelocity(u);
  
                    // x = nodesFluid_[connec(k)] -> getCoordinates();
                    // Acc(0) = nodesFluid_[connec(k)] -> getAcceleration(0);
                    // Acc(1) = nodesFluid_[connec(k)] -> getAcceleration(1);

                    //if(rank == 0) std::cout << "Accel " << u[0] << " " << u[1] << " " << std::endl; 
                };
            };
        };
    };

};

//------------------------------------------------------------------------------
//---------SETS FLUID AND SOLID MODELS AND GETS ITS BASIC INFORMATIONS----------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void FSInteraction<DIM,DEG>::transferSolidVelocityArlequin(){
    

    for (int i = 0; i < numInterfaces; i++){

        int interf = groupInterfaces[i];       

        for (int ibound = 0; ibound < numElemArlequinBoundaryFine; ibound++){
            if (boundaryArlequinFine_[ibound] -> getBoundaryGroup() == interf){
                
                VecInt connec = boundaryArlequinFine_[ibound] -> getBoundaryConnectivity();
              
                VecDouble u(DIM);
                int nBdNodes = arlequinModel->MeshVec()[1]->NBdNodes();
                for (int k = 0; k < nBdNodes; k++){
                    
                    int elem = nodesArlequinFine_[connec[k]] -> getNodalElemCorrespondence();
                    VecDouble xsi = nodesArlequinFine_[connec[k]] ->  getNodalXsiCorrespondence();
                    
                    if (rank == 0) 
                        getinterpolatedvelocity_(&u[0],&u[1],&elem,&xsi[0]);
                        
                    MPI_Bcast(&u[0],1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
                    MPI_Bcast(&u[1],1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
               
                    nodesArlequinFine_[connec[k]] -> setVelocity(u);
  
                    // x = nodesFluid_[connec(k)] -> getCoordinates();
                    // Acc(0) = nodesFluid_[connec(k)] -> getAcceleration(0);
                    // Acc(1) = nodesFluid_[connec(k)] -> getAcceleration(1);

                    //if(rank == 0) std::cout << "Accel " << u[0] << " " << u[1] << " " << std::endl; 
                };
            };
        };
    };
 

};

//------------------------------------------------------------------------------
//---------SETS FLUID AND SOLID MODELS AND GETS ITS BASIC INFORMATIONS----------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void FSInteraction<DIM,DEG>::transferFluidLoad(){
    
    clearcouplingloads_();
    
    for (int iInterf = 0; iInterf < numInterfaces; iInterf++){
        for (int isolid = 0; isolid < numNodesSolid; isolid++){
            
            int ielem = nodesSolid_[iInterf][isolid] -> getNodalElemCorrespondence();
            VecDouble xsi = nodesSolid_[iInterf][isolid] -> getNodalXsiCorrespondence();

            VecDouble load(DIM);
            elementsFluid_[ielem] -> getBoundaryLoad(xsi,load);

            int inode = isolid+1;
            setcouplingload_(&load[0],&load[1],&inode);
                       
        };//isolid
    };//iInterf
    
};

//------------------------------------------------------------------------------
//---------SETS FLUID AND SOLID MODELS AND GETS ITS BASIC INFORMATIONS----------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void FSInteraction<DIM,DEG>::transferArlequinLoad(){
    
    clearcouplingloads_();
    // std::cout << "AAAAA 1 "<< std::endl;
    for (int iInterf = 0; iInterf < numInterfaces; iInterf++){
        for (int isolid = 0; isolid < numNodesSolid; isolid++){
            // std::cout << "AAAAA 2 "<< std::endl;
            int ielem = nodesSolid_[iInterf][isolid] -> getNodalElemCorrespondence();
            // std::cout << "AAAAA 3 " << ielem << " " << isolid << " " << std::endl;
            VecDouble xsi = nodesSolid_[iInterf][isolid] -> getNodalXsiCorrespondence();
            // std::cout << "AAAAA 4 "<< ielem << " " << xsi[0] << " " << xsi[1] << " " << std::sqrt(std::inner_product(xsi,xsi,xsi,0.0L)) << std::endl;
            VecDouble load(DIM);
            if (xsi[0] < 3) elementsArlequinFine_[ielem] -> getBoundaryLoad(xsi,load);
            // std::cout << "AAAAA 5 "<< std::endl;
            int inode = isolid+1;
            setcouplingload_(&load[0],&load[1],&inode);
                       
        };//isolid
    };//iInterf
    
};


//------------------------------------------------------------------------------
//----------------SOLVES THE FLUID-STRUCTURE INTERACTION PROBLEM----------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void FSInteraction<DIM,DEG>::solveFSIProblem(int numTimeSteps){

    std::string om = "omega.txt";
    std::ofstream saidaOmega(om.c_str());


    double sizeSolid = 3 * numNodesSolid;

    double *X_k;
    double *Y_k;
    double *deltaXi;
    double *deltaXii;  

    X_k = new double[3* numNodesSolid]();
    Y_k = new double[3* numNodesSolid]();
    deltaXi = new double[3* numNodesSolid]();
    deltaXii = new double[3* numNodesSolid]();
    

    double omega = 1.;
    double mu = 0.;

    double &alpha_f = fluidModel->fProbParameters.getAlphaF();
    double &alpha_m = fluidModel->fProbParameters.getAlphaM();
    double &gamma = fluidModel->fProbParameters.getGamma();

    for (int iTimeStep = 0; iTimeStep < numTimeSteps; iTimeStep++){  

        std::clock_t t1 = std::clock();

        if (rank == 0) {
            std::cout << std::endl;
            std::cout << "****************************************"
                      << "****************************************"
                      << std::endl;
            std::cout << 
                "                               TIME STEP = "
                      << iTimeStep << std::endl;
        };

        // if (iTimeStep == 20){
        //     double integ = 1.0;
        //     fluidModel->ProblemParameters.setSpectralRadius(integ);  
        // } 
        
        // //SOMENTE PARA EXEMPLO DA CAVIDADE - INICIO
        // for (int ibound = 0; ibound < numElemFluidBoundary; ibound++){
            
        //     Boundary::BoundConnect connectB;
        //     connectB = boundaryFluid_[ibound] -> getBoundaryConnectivity();
        //     int no1 = connectB(0);
        //     int no2 = connectB(1);
        //     int no3 = connectB(2);
            
        //     if (boundaryFluid_[ibound] -> getConstrain(0) == 1){
                
        //         double value = boundaryFluid_[ibound] -> getConstrainValue(0) * 
        //             (1. - cos(0.4 * pi * dTime * iTimeStep));
        //         nodesFluid_[no1] -> setConstrains(0,boundaryFluid_[ibound] -> 
        //                                           getConstrain(0),value);
        //         nodesFluid_[no2] -> setConstrains(0,boundaryFluid_[ibound] -> 
        //                                           getConstrain(0),value);
        //         nodesFluid_[no3] -> setConstrains(0,boundaryFluid_[ibound] ->
        //                                           getConstrain(0),value);
        //     };
        // };
        // //SOMENTE PARA EXEMPLO DA CAVIDADE - FIM

        for (int i = 0; i < numNodesFluid; i++){
            VecDouble accel(DIM), u(DIM), uprev(DIM);
            //Compute acceleration
            u[0] = nodesFluid_[i] -> getVelocity(0);
            u[1] = nodesFluid_[i] -> getVelocity(1);
            
            nodesFluid_[i] -> setPreviousVelocity(u);

            accel[0] = nodesFluid_[i] -> getAcceleration(0);
            accel[1] = nodesFluid_[i] -> getAcceleration(1);
            
            nodesFluid_[i] -> setPreviousAcceleration(accel);
            accel[0] *= (gamma - 1.) / gamma;
            accel[1] *= (gamma - 1.) / gamma;

            nodesFluid_[i] -> setAcceleration(accel);
        };

        if (rank == 0) updateqsrs_();

        for (int i = 0; i < numNodesFluid; i++){
            VecDouble x = nodesFluid_[i] -> getCoordinates();
            nodesFluid_[i] -> setPreviousCoordinates(0,x[0]);
            nodesFluid_[i] -> setPreviousCoordinates(1,x[1]);
        };

        for (int i = 0; i < numNodesSolid; i++){
            int dof = 3*i+1;
            if (rank == 0) getposition_(&dof,&Y_k[3*i  ]);
            dof++;
            if (rank == 0) getposition_(&dof,&Y_k[3*i+1]);
            MPI_Bcast(&Y_k[3*i  ],1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
            MPI_Bcast(&Y_k[3*i+1],1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
        };
      
        saidaOmega << std::endl << "Passo de tempo " << iTimeStep << std::endl;

        //COMPUTING PREDICTOR
        for (int i = 0; i < numNodesSolid; i++){
            if (rank == 0) {
                double v_ = 0.;
                double v_prev = 0.;
                int dof = 3*i+1;
                getposition_(&dof,&Y_k[3*i  ]);
                getvelocity_(&dof,&v_);
                getpreviousvelocity_(&dof,&v_prev);

                Y_k[3*i  ] += dTime * (1.5 * v_ - 0.5 * v_prev);

                setposition_(&dof,&Y_k[3*i  ]);
                
                dof++;

                getposition_(&dof,&Y_k[3*i+1]);
                getvelocity_(&dof,&v_);
                getpreviousvelocity_(&dof,&v_prev);

                Y_k[3*i+1] += dTime * (1.5 * v_ - 0.5 * v_prev);

                setposition_(&dof,&Y_k[3*i+1]);
            };
            MPI_Bcast(&Y_k[3*i  ],1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
            MPI_Bcast(&Y_k[3*i+1],1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
        };
           
        //X_k = Y_k;

        double residual = 1.e10;
        int iterations = 0;


        //Gauss-Seidel iterations
        while ((residual > 1.e-4) && (iterations < 30)){
            if(rank == 0) {std::cout << 
                    "........................... GAUSS-SEIDEL ITERATION "
                                     << iterations + 1 << 
                    " ..........................." << std:: endl;};

            for (int i=0; i<sizeSolid; i++) X_k[i] = Y_k[i];
            
            updateFluidMesh();
             
            transferSolidVelocity();
            
            fluidModel->solveFSIFluid(3, 1.e-5, 2);
            
            if (rank == 0) transferFluidLoad();
            
            if (rank == 0) solveframestructure_(&iTimeStep);
            
            
            for (int i = 0; i < numNodesSolid; i++){
                int dof = 3*i+1;
                if (rank == 0) getposition_(&dof,&Y_k[3*i  ]);
                dof++;
                if (rank == 0) getposition_(&dof,&Y_k[3*i+1]);
                MPI_Bcast(&Y_k[3*i  ],1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
                MPI_Bcast(&Y_k[3*i+1],1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
            };
            
            for (int i=0; i<sizeSolid; i++) deltaXii[i] = X_k[i] - Y_k[i];
            
            //Aitken Relaxation
            if (iterations > 0){
                double aux1 = 0.;
                double aux2 = 0.;

                for (int i=0; i<sizeSolid; i++){
                    aux1 += (deltaXi[i]-deltaXii[i]) * deltaXii[i];
                    aux2 += (deltaXi[i]-deltaXii[i]) * (deltaXi[i]-deltaXii[i]);
                }

                mu = mu + (mu - 1.0) * aux1 / sqrt(aux2);
            };        
            
            residual = 0.;

            for (int i=0; i<sizeSolid; i++) {
                deltaXi[i] = deltaXii[i];
                residual += deltaXii[i] * deltaXii[i];
            }
            
            residual = sqrt(residual);

            omega = 1. - mu;

            for (int i=0; i<sizeSolid; i++) Y_k[i] = (1. - omega) * X_k[i] + omega * Y_k[i];
            

            for (int i = 0; i < numNodesSolid; i++){
                if (rank == 0) {
                    int dof = 3*i+1;
                    setposition_(&dof,&Y_k[3*i  ]);
                    dof++;
                    setposition_(&dof,&Y_k[3*i+1]);
                };
                MPI_Bcast(&Y_k[3*i  ],1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
                MPI_Bcast(&Y_k[3*i+1],1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
            };
       
            saidaOmega << std::scientific << omega << " " << std::scientific << residual << " " << iterations << std::endl;

            iterations++;

            if(rank == 0) std::cout << "GAUSS-SEIDEL OMEGA = " 
                                    << omega <<  std::endl;

            if(rank == 0) std::cout << "GAUSS-SEIDEL RESIDUAL = " 
                                    << std::scientific << residual << std::endl;

        };
        
        if (rank == 0) updatesolid_(&iTimeStep);

        //Updates SUPG Parameter
        // for (int i = 0; i < numElemFluid; i++){
        //     elementsFluid_[i] -> getParameterSUPG();
        // };
                
        std::clock_t t2 = std::clock();

        if (rank == 0) {
            std::cout << "****************************************"
                      << "****************************************"
                      << std::endl;
            std::cout << "********************** PROCESSING TIME = " << 
                std::fixed << 1000.*(t2-t1)/CLOCKS_PER_SEC/1000.
                      << " seconds **********************" 
                       << std::endl;
            std::cout << "****************************************"
                      << "****************************************"
                      << std::endl;
        };

        // Printing Results
        if (rank == 0) {
            if(iTimeStep % 1 == 0){
                // fluidModel->printResults(iTimeStep);
                printstructure_();
            };
        };


    };//Time Steps

    delete [] X_k;
    delete [] Y_k;
    delete [] deltaXi;
    delete [] deltaXii;

    return;
};



//------------------------------------------------------------------------------
//----------------SOLVES THE FLUID-STRUCTURE INTERACTION PROBLEM----------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void FSInteraction<DIM,DEG>::solveFSIProblemArlequin(int numTimeSteps){


    std::string om = "omega.txt";
    std::ofstream saidaOmega(om.c_str());


    double sizeSolid = 3 * numNodesSolid;

    double *X_k;
    double *Y_k;
    double *deltaXi;
    double *deltaXii;  

    X_k = new double[3* numNodesSolid]();
    Y_k = new double[3* numNodesSolid]();
    deltaXi = new double[3* numNodesSolid]();
    deltaXii = new double[3* numNodesSolid]();

    double omega = 1.;
    double mu = 0.;

    if (rank == 0) {
        arlequinModel->printResultsCoarse(0);
        arlequinModel->printResultsFine(0);
        printstructure_();
    };

    double &alpha_f = arlequinModel->MeshVec()[1]->fProbParameters.getAlphaF();
    double &alpha_m = arlequinModel->MeshVec()[1]->fProbParameters.getAlphaM();
    double &gamma = arlequinModel->MeshVec()[1]->fProbParameters.getGamma();

    for (int iTimeStep = 0; iTimeStep < numTimeSteps; iTimeStep++){  

        std::clock_t t1 = std::clock();

        if (rank == 0) {
            std::cout << std::endl;
            std::cout << "****************************************"
                      << "****************************************"
                      << std::endl;
            std::cout << 
                "                               TIME STEP = "
                      << iTimeStep << std::endl;
        };
        
        // //SOMENTE PARA EXEMPLO DA CAVIDADE - INICIO
        // for (int ibound = 0; ibound < numElemFluidBoundary; ibound++){
            
        //     Boundary::BoundConnect connectB;
        //     connectB = boundaryFluid_[ibound] -> getBoundaryConnectivity();
        //     int no1 = connectB(0);
        //     int no2 = connectB(1);
        //     int no3 = connectB(2);
            
        //     if (boundaryFluid_[ibound] -> getConstrain(0) == 1){
                
        //         double value = boundaryFluid_[ibound] -> getConstrainValue(0) * 
        //             (1. - cos(0.4 * pi * dTime * iTimeStep));
        //         nodesFluid_[no1] -> setConstrains(0,boundaryFluid_[ibound] -> 
        //                                           getConstrain(0),value);
        //         nodesFluid_[no2] -> setConstrains(0,boundaryFluid_[ibound] -> 
        //                                           getConstrain(0),value);
        //         nodesFluid_[no3] -> setConstrains(0,boundaryFluid_[ibound] ->
        //                                           getConstrain(0),value);
        //     };
        // };
        // //SOMENTE PARA EXEMPLO DA CAVIDADE - FIM

        // std::cout << "AQUI1 " << rank << " " << gamma << std::endl;
        // MPI_Barrier(PETSC_COMM_WORLD);

        // if (iTimeStep == 10){
        //     double spec = 0.0;
        //     arlequinModel->MeshVec()[1]->ProblemParameters.setSpectralRadius(spec);
        //     arlequinModel->MeshVec()[0]->ProblemParameters.setSpectralRadius(spec);
        //     std::cout << "AQUI " << rank << std::endl;
        // }

        for (int i = 0; i < numNodesArlequinCoarse; i++){
            VecDouble accel(DIM), u(DIM), uprev(DIM);
            
            //Compute acceleration
            u[0] = nodesArlequinCoarse_[i] -> getVelocity(0);
            u[1] = nodesArlequinCoarse_[i] -> getVelocity(1);

            nodesArlequinCoarse_[i] -> setPreviousVelocity(u);
            
            accel[0] = nodesArlequinCoarse_[i] -> getAcceleration(0);
            accel[1] = nodesArlequinCoarse_[i] -> getAcceleration(1);
            
            nodesArlequinCoarse_[i] -> setPreviousAcceleration(accel);

            accel[0] *= (gamma - 1.) / gamma;
            accel[1] *= (gamma - 1.) / gamma;
            
            nodesArlequinCoarse_[i] -> setAcceleration(accel);
        };

        for (int i = 0; i < numNodesArlequinFine; i++){
            VecDouble accel(DIM), u(DIM), uprev(DIM), lag(DIM);
            
            //Compute acceleration
            u[0] = nodesArlequinFine_[i] -> getVelocity(0);
            u[1] = nodesArlequinFine_[i] -> getVelocity(1);

            nodesArlequinFine_[i] -> setPreviousVelocity(u);
            
            accel[0] = nodesArlequinFine_[i] -> getAcceleration(0);
            accel[1] = nodesArlequinFine_[i] -> getAcceleration(1);
            
            nodesArlequinFine_[i] -> setPreviousAcceleration(accel);

            accel[0] *= (gamma - 1.) / gamma;
            accel[1] *= (gamma - 1.) / gamma;
            
            nodesArlequinFine_[i] -> setAcceleration(accel);
        };

        // std::cout << "AQUI2 " << rank << std::endl;
        // MPI_Barrier(PETSC_COMM_WORLD);

        if (rank == 0) updateqsrs_();

        // std::cout << "AQUI3 " << rank << std::endl;
        // MPI_Barrier(PETSC_COMM_WORLD);

        //Compute Qs and Rs for mesh moving problem
        // for (int i = 0; i < numNodesArlequinFine; i++){
        //     double x = nodesArlequinFine_[i] -> getCoordinateValue(0);
        //     double y = nodesArlequinFine_[i] -> getCoordinateValue(1);
        //     double xp = nodesArlequinFine_[i] -> getPreviousCoordinateValue(0);
        //     double yp = nodesArlequinFine_[i] -> getPreviousCoordinateValue(1);
        //     double vx = nodesArlequinFine_[i] -> getMeshVelocity(0);
        //     double vy = nodesArlequinFine_[i] -> getMeshVelocity(1);
        //     double ax = nodesArlequinFine_[i] -> getMeshAcceleration(0);
        //     double ay = nodesArlequinFine_[i] -> getMeshAcceleration(1);

        //     double accelx = (x - xp) / (0.25 * dTime * dTime) - vx / (0.25 * dTime) - ax * (0.5/0.25 - 1.0);
        //     double accely = (y - yp) / (0.25 * dTime * dTime) - vy / (0.25 * dTime) - ay * (0.5/0.25 - 1.0);

        //     nodesArlequinFine_[i] -> setMeshAccelerationComponent(0,accelx);
        //     nodesArlequinFine_[i] -> setMeshAccelerationComponent(1,accely);

        //     double velx = 0.5 * dTime * accelx + vx + dTime * (1.0 - 0.5) * ax;
        //     double vely = 0.5 * dTime * accely + vy + dTime * (1.0 - 0.5) * ay;

        //     nodesArlequinFine_[i] -> setMeshVelocityComponent(0,velx);
        //     nodesArlequinFine_[i] -> setMeshVelocityComponent(1,vely);
        // }


        for (int i = 0; i < numNodesArlequinFine; i++){
            VecDouble x = nodesArlequinFine_[i] -> getCoordinates();
            nodesArlequinFine_[i] -> setPreviousCoordinates(0,x[0]);
            nodesArlequinFine_[i] -> setPreviousCoordinates(1,x[1]);
        };


        // std::cout << "AQUI4 " << rank << std::endl;
        // MPI_Barrier(PETSC_COMM_WORLD);

        for (int i = 0; i < numNodesSolid; i++){
            int dof = 3*i+1;
            if (rank == 0) getposition_(&dof,&Y_k[3*i  ]);
            dof++;
            if (rank == 0) getposition_(&dof,&Y_k[3*i+1]);
            MPI_Bcast(&Y_k[3*i  ],1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
            MPI_Bcast(&Y_k[3*i+1],1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
        };
      
        saidaOmega << std::endl << "Passo de tempo " << iTimeStep << std::endl;

        // std::cout << "AQUI5 " << rank << std::endl;
        // MPI_Barrier(PETSC_COMM_WORLD);
        //COMPUTING PREDICTOR
        for (int i = 0; i < numNodesSolid; i++){
            if (rank == 0) {
                double v_ = 0.;
                double v_prev = 0.;
                int dof = 3*i+1;
                getposition_(&dof,&Y_k[3*i  ]);
                getvelocity_(&dof,&v_);
                getpreviousvelocity_(&dof,&v_prev);

                Y_k[3*i  ] += dTime * (1.5 * v_ - 0.5 * v_prev);

                setposition_(&dof,&Y_k[3*i  ]);
                
                dof++;

                getposition_(&dof,&Y_k[3*i+1]);
                getvelocity_(&dof,&v_);
                getpreviousvelocity_(&dof,&v_prev);

                Y_k[3*i+1] += dTime * (1.5 * v_ - 0.5 * v_prev);

                setposition_(&dof,&Y_k[3*i+1]);
            };
            MPI_Bcast(&Y_k[3*i  ],1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
            MPI_Bcast(&Y_k[3*i+1],1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
        };
           
        // std::cout << "AQUI6 " << rank << std::endl;
        // MPI_Barrier(PETSC_COMM_WORLD);
        //X_k = Y_k;

        double residual = 1.e10;
        int iterations = 0;

        //Gauss-Seidel iterations
        while (((residual > 1.e-5) && (iterations < 5))){//} || (iterations < 3)) {
            if(rank == 0) {std::cout << 
                    "........................... GAUSS-SEIDEL ITERATION "
                                     << iterations + 1 << 
                    " ..........................." << std:: endl;};

            for (int i=0; i<sizeSolid; i++) X_k[i] = Y_k[i];
        
            // std::cout << "AQUI-2 " << rank << std::endl;            
            //     MPI_Barrier(PETSC_COMM_WORLD);
            updateArlequinMesh();
            
            // std::cout << "AQUI7 " << rank << std::endl;
            // MPI_Barrier(PETSC_COMM_WORLD); 
            // std::cout << "AQUI-1 " << rank << std::endl;            
            //     MPI_Barrier(PETSC_COMM_WORLD);
            transferSolidVelocityArlequin();

            // std::cout << "AQUI8 " << rank << std::endl;
            // MPI_Barrier(PETSC_COMM_WORLD);
                // std::cout << "AQUI0 " << rank << std::endl;            
                // MPI_Barrier(PETSC_COMM_WORLD);
            // if (iTimeStep < 50){
            arlequinModel->solveFSIArlequin(4, 1.e-3, 2, iTimeStep);

                // std::cout << "AQUI4 " << rank << std::endl;
                // MPI_Barrier(PETSC_COMM_WORLD);
            // }else{
            //     arlequinModel->solveFSIArlequin(5, 1.e-3, 2, iTimeStep);
            // }
            
            // std::cout << "AQUI9 " << rank << std::endl;
            // MPI_Barrier(PETSC_COMM_WORLD);

            if (rank == 0) transferArlequinLoad();
            // if ((rank == 0) && (iTimeStep < 250)) {
            //     int inode = 43;
            //     double load[2];
            //     load[0] = 0.;
            //     load[1] = 500.;
            //     setcouplingload_(&load[0],&load[1],&inode);
            //     inode = 39;
            //     load[1] = -500;
            //     setcouplingload_(&load[0],&load[1],&inode);
            // }
            // std::cout << "AQUI5 " << rank << std::endl;
            // std::cout << "AQUI10 " << rank << std::endl;
            // MPI_Barrier(PETSC_COMM_WORLD);

            //  if (iTimeStep > 3){
            // if (rank == 0) solveframestructure_(&iTimeStep);
            // };
            // std::cout << "AQUI6 " << rank << std::endl;
            // std::cout << "AQUI11 " << rank << std::endl;
            // MPI_Barrier(PETSC_COMM_WORLD);
            
            for (int i = 0; i < numNodesSolid; i++){
                int dof = 3*i+1;
                if (rank == 0) getposition_(&dof,&Y_k[3*i  ]);
                dof++;
                if (rank == 0) getposition_(&dof,&Y_k[3*i+1]);
                // std::cout << "SASDASD " << Y_k[3*i  ] << " " << Y_k[3*i+1] << " " << rank << std::endl;
                MPI_Bcast(&Y_k[3*i  ],1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
                MPI_Bcast(&Y_k[3*i+1],1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
            };
// std::cout << "AQUI7 " << rank << std::endl;
            // std::cout << "AQUI12 " << rank << std::endl;
            // MPI_Barrier(PETSC_COMM_WORLD);
            
            for (int i=0; i<sizeSolid; i++) deltaXii[i] = X_k[i] - Y_k[i];
            
            //Aitken Relaxation
             if (iterations > 0){
                double aux1 = 0.;
                double aux2 = 0.;

                for (int i=0; i<sizeSolid; i++){
                    aux1 += (deltaXi[i]-deltaXii[i]) * deltaXii[i];
                    aux2 += (deltaXi[i]-deltaXii[i]) * (deltaXi[i]-deltaXii[i]);
                }
                // std::cout << "AASDASD " << aux1 << " " << aux2 << std::endl;
                mu = mu + (mu - 1.0) * aux1 / sqrt(aux2);
            };        
            
            residual = 0.;

            for (int i=0; i<sizeSolid; i++) {
                deltaXi[i] = deltaXii[i];
                residual += deltaXii[i] * deltaXii[i];
            }
            
            residual = sqrt(residual);

            omega = 1. - mu;

            for (int i=0; i<sizeSolid; i++) Y_k[i] = (1. - omega) * X_k[i] + omega * Y_k[i];
            

            for (int i = 0; i < numNodesSolid; i++){
                if (rank == 0) {
                    int dof = 3*i+1;
                    setposition_(&dof,&Y_k[3*i  ]);
                    dof++;
                    setposition_(&dof,&Y_k[3*i+1]);
                };
                MPI_Bcast(&Y_k[3*i  ],1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
                MPI_Bcast(&Y_k[3*i+1],1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
            };
       
            saidaOmega << std::scientific << omega << " " << std::scientific << residual << " " << iterations << std::endl;

            iterations++;
// std::cout << "AQUI8 " << rank << std::endl;
            if(rank == 0) std::cout << "GAUSS-SEIDEL OMEGA = " 
                                    << omega <<  std::endl;

            if(rank == 0) std::cout << "GAUSS-SEIDEL RESIDUAL = " 
                                    << std::scientific << residual << std::endl;

        };
        
        if (rank == 0) updatesolid_(&iTimeStep);

        //Updates SUPG Parameter
        // for (int i = 0; i < numElemFluid; i++){
        //     elementsFluid_[i] -> getParameterSUPG();
        // };
                
        std::clock_t t2 = std::clock();

        if (rank == 0) {
            std::cout << "****************************************"
                      << "****************************************"
                      << std::endl;
            std::cout << "********************** PROCESSING TIME = " << 
                std::fixed << 1000.*(t2-t1)/CLOCKS_PER_SEC/1000.
                      << " seconds **********************" 
                       << std::endl;
            std::cout << "****************************************"
                      << "****************************************"
                      << std::endl;
        };

        // Printing Results
        if (rank == 0) {
            if(iTimeStep % 1 == 0){
                arlequinModel->printResultsCoarse(iTimeStep);
                arlequinModel->printResultsFine(iTimeStep);
                printstructure_();
            };
        };


    };//Time Steps

    delete [] X_k;
    delete [] Y_k;
    delete [] deltaXi;
    delete [] deltaXii;

    return;
};



template class FSInteraction<2,1>;
template class FSInteraction<2,2>;
template class FSInteraction<3,1>;
template class FSInteraction<3,2>;