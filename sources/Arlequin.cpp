#include "Arlequin.h"
#include "hdf5.h"

//------------------------------------------------------------------------------
//--------------------------------IMPLEMENTATION--------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//-------------------------COMPUTE ELEMENT REGIONS/BOXES------------------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Arlequin<DIM,DEG>::setElementBoxes() {
    
    VecInt connec;
    VecDouble xk(DIM), Xk(DIM);

    //Compute element boxes for coarse model
    //Only function for straight elements
    for (int jel = 0; jel < numElemCoarse; jel++){
        connec = elementsCoarse_[jel] -> getConnectivity();
        VecDouble x1 = (*nodesCoarse_)[connec[0]] -> getCoordinates();
        VecDouble x2 = (*nodesCoarse_)[connec[1]] -> getCoordinates();
        VecDouble x3 = (*nodesCoarse_)[connec[2]] -> getCoordinates();      

        xk[0] = std::min(x1[0],std::min(x2[0], x3[0]));
        xk[1] = std::min(x1[1],std::min(x2[1], x3[1]));

        Xk[0] = std::max(x1[0],std::max(x2[0], x3[0]));
        Xk[1] = std::max(x1[1],std::max(x2[1], x3[1]));        
        
        elementsCoarse_[jel] -> setIntersectionParameters(xk, Xk);
    };

    //Compute element boxes for fine model
    //Only function for straight elements
    for (int jel = 0; jel < numElemFine; jel++){
        connec = elementsFine_[jel] -> getConnectivity();
        VecDouble x1 = (*nodesFine_)[connec[0]] -> getCoordinates();
        VecDouble x2 = (*nodesFine_)[connec[1]] -> getCoordinates();
        VecDouble x3 = (*nodesFine_)[connec[2]] -> getCoordinates();      

        xk[0] = std::min(x1[0],std::min(x2[0], x3[0]));
        xk[1] = std::min(x1[1],std::min(x2[1], x3[1]));

        Xk[0] = std::max(x1[0],std::max(x2[0], x3[0]));
        Xk[1] = std::max(x1[1],std::max(x2[1], x3[1]));        
        
        elementsFine_[jel] -> setIntersectionParameters(xk, Xk);
    };

    return;
};


//------------------------------------------------------------------------------
//-------------------COMPUTE NODAL CORRESPONDECE WITH ELEMENTS------------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Arlequin<DIM,DEG>::searchNodeCorrespondence(VecDouble &x,std::vector<Nodes *> nodes, 
                                           std::vector<Elements *> elements, 
                                           int numElem, int &elCorr, VecDouble &xsiCorr, int elSearch){
    
    ShapeFunction<DIM,DEG> shapeQuad;
    VecDouble phi_(nElNodes);

    MatrixDouble ainv(DIM,DIM);
    
    VecDouble xsiCC(3);
    std::pair<VecDouble,VecDouble> XK;

    elCorr = 150000;
    VecDouble xsi(DIM);
    VecDouble x_(DIM);
    VecDouble deltaX(DIM);
    VecDouble deltaXsi(DIM);
    bool flg = true;
    VecInt connec = elements[elSearch] -> getConnectivity();
    
    xsiCC.fill(1.e10);
    xsiCorr.fill(1.e50);
    xsi.fill(1./3.);
    x_.fill(0.);
    
    shapeQuad.evaluate(xsi,phi_);

    for (int i = 0; i < nElNodes; i++){
        VecDouble xint = nodes[connec[i]] -> getCoordinates();
        x_ += xint * phi_[i];
    };

    double error = 1.e6;
    int iterations = 0;

    while ((error > 1.e-8) && (iterations < 4)) {
        
        iterations++;
        
        deltaX = x - x_;
        deltaXsi.setZero();
        
        double djac_ = 0.;
        elements[elSearch] -> getJacobianMatrix(xsi,ainv,djac_);

        // for (int i = 0; i < DIM; i++)
        //     for (int j = 0; j < DIM; j++)
        //         deltaXsi[i] += ainv(j,i) * deltaX[j];
        deltaXsi = ainv.transpose()*deltaX;    

        xsi += deltaXsi;
        x_.setZero();
        
        shapeQuad.evaluate(xsi,phi_);
        
        for (int i=0; i<nElNodes; i++){
            VecDouble xint = nodes[connec[i]] -> getCoordinates();
            x_ += xint * phi_[i];
        };

        error = std::sqrt(deltaXsi[0]*deltaXsi[0] + deltaXsi[1]*deltaXsi[1]);
    };
    
    double t1 = -1.e-2;
    double t2 =  1. - t1;
    
    xsiCC[0] = xsi[0];
    xsiCC[1] = xsi[1];       
    xsiCC[2] = 1. - xsiCC[0] - xsiCC[1];

    if ((xsiCC[0] >= t1) && (xsiCC[1] >= t1) && (xsiCC[2] >= t1) &&
        (xsiCC[0] <= t2) && (xsiCC[1] <= t2) && (xsiCC[2] <= t2)){

        xsiCorr[0] = xsi[0]; xsiCorr[1] = xsi[1];
        elCorr = elSearch;
        // return;
    } else {

        int nEl;
        // if (iTimeStep < 2){
            nEl = numElem;
        // } else {
        //     nEl = elements[elSearch] -> getNumberOfNeighborElements(); 
        // }   

        for (int jel = 0; jel < nEl; jel++){

            // if (iTimeStep < 2){
                connec = elements[jel] -> getConnectivity();
            // } else {
            //     connec = elements[elements[elSearch] -> getNeighborElement(jel)] -> getConnectivity();
            // }   

            //get boxes information        
            XK = elements[jel] -> getXIntersectionParameter();

            //Chech if the node is inside the element box
            if ((x[0] < XK.first[0]) || (x[0] > XK.second[0]) ||
                (x[1] < XK.first[1]) || (x[1] > XK.second[1])) continue;
            
            //Compute nodal correspondence
            for (int i = DIM+1; i--; ) xsiCC[i] = 1.e10;
    
            for (int i = DIM; i--; ){
                xsi[i] = 1. / 3.;
                x_[i] = 0.;
            }

            shapeQuad.evaluate(xsi,phi_);

            for (int i = 0; i < nElNodes; i++){
                VecDouble xint = nodes[connec[i]] -> getCoordinates();
                for (int k = DIM; k--; )
                    x_[k] += xint[k] * phi_[i];
            };
            
            double error = 1.e6;
            int iterations = 0;

            while ((error > 1.e-8) && (iterations < 4)) {
                
                iterations++;

                for (int k = DIM; k--; ){
                    deltaX[k] = x[k] - x_[k];
                    deltaXsi[k] = 0.;
                }
                
                double djac_ = 0.;
                elements[jel] -> getJacobianMatrix(xsi,ainv,djac_);
            
                for (int i = 0; i < DIM; i++)
                    for (int j = 0; j < DIM; j++)
                        deltaXsi[i] += ainv(j,i) * deltaX[j];
                
                for (int k = DIM; k--; ){
                    xsi[k] += deltaXsi[k];
                    x_[k] = 0.;
                }
            
                shapeQuad.evaluate(xsi,phi_);
                
                for (int i=0; i<nElNodes; i++){
                    VecDouble xint = nodes[connec[i]] -> getCoordinates();
                    for (int k = DIM; k--; ) x_[k] += xint[k] * phi_[i];
                };
                        
                error = std::sqrt(deltaXsi[0]*deltaXsi[0] + deltaXsi[1]*deltaXsi[1]);
            };
            
            double t1 = -1.e-2;
            double t2 =  1. - t1;
            
            xsiCC[0] = xsi[0];
            xsiCC[1] = xsi[1];       
            xsiCC[2] = 1. - xsiCC[0] - xsiCC[1];

            if ((xsiCC[0] >= t1) && (xsiCC[1] >= t1) && (xsiCC[2] >= t1) &&
                (xsiCC[0] <= t2) && (xsiCC[1] <= t2) && (xsiCC[2] <= t2)){

                xsiCorr[0] = xsi[0]; xsiCorr[1] = xsi[1];
                elCorr = jel;
                break;
            }
        }
    };

    if (fabs(xsi[0]) > 2.) std::cout << "PROBEM SEARCHING NODE CORRESPONDENCE " 
                                     << std::endl;

    return;
};

//------------------------------------------------------------------------------
//--------COMPUTE NODAL CORRESPONDECE OF FINE NODES WITH COARSE ELEMENTS--------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Arlequin<DIM,DEG>::setNodalCorrespondenceFine() {

    //FINE MESH
    double &alpha_f = parametersFine -> getAlphaF();

    for (int inode = 0; inode < numNodesGlueZoneFine; inode++) {
        
        VecDouble x = (*nodesFine_)[nodesGlueZoneFine_[inode]] -> getCoordinates();

        int elCorr = 0;
        VecDouble xsiCorr(DIM);
        // searchNodeCorrespondence(x, nodesCoarse_, elementsCoarse_[(*nodesFine_)[nodesGlueZoneFine_[inode]]->getNodalElemCorrespondence()], 
        //                          1,elCorr,xsiCorr, );

        // double t1 = -1.e-2;
        // double t2 =  1. - t1;

        // if ((xsiCorr[0] >= t1) && (xsiCorr[1] >= t1) && (1. - xsiCorr[0] - xsiCorr[1] >= t1) &&
        //     (xsiCorr[0] <= t2) && (xsiCorr[1] <= t2) && (1. - xsiCorr[0] - xsiCorr[1] <= t2)){
            
        //     (*nodesFine_)[nodesGlueZoneFine_[inode]] -> setNodalCorrespondence(elCorr,xsiCorr);
        // } else {
        searchNodeCorrespondence(x, (*nodesCoarse_), elementsCoarse_, 
                                 elementsCoarse_.size(),elCorr,xsiCorr,(*nodesFine_)[nodesGlueZoneFine_[inode]] -> getNodalElemCorrespondence());
        (*nodesFine_)[nodesGlueZoneFine_[inode]] -> setNodalCorrespondence(elCorr,xsiCorr);             
        // }

        // std::cout << "INODE " << nodesGlueZoneFine_[inode] << " " << elCorr << " " << xsiCorr[0] << " " << xsiCorr[1] << std::endl;
   
        // std::cout << "CORRESP " << elCorr << " " << corresp.first << std::endl 
                  // << xsiCorr[0] << " " << xsiCorr[1] << " " << corresp.second[0] << " " << corresp.second[1] << std::endl;

        if (elCorr == 150000) {
            std::cout << "PROBLEM IN NODAL CORRESPONDENCE OF FINE NODE 1" 
                      << nodesGlueZoneFine_[inode] << std::endl;
            exit(0);
        };
            
        // std::cout << "corresp " << corresp.first 
        //           << " " << corresp.second(0) << " " << corresp.second(1)
        //           << std::endl;

    };

    // for (int i=0; i<numNodesFine; i++){
    //     int elem = (*nodesFine_)[i] -> getNodalElemCorrespondence();
    //     xsi = (*nodesFine_)[i] -> getNodalXsiCorrespondence();
    //     std::cout << "node " << i << " elem " << elem 
    //               << " " << xsi(0) << " " << xsi(1) << std::endl;
    // };

    //Compute correspondence of integration points
    //int numberIntPoints = elementsFine_[0] -> getNumberOfIntegrationPoints();
    //if (rank == 0) std::cout << "Int Points " << numberIntPoints << std::endl;

    for (int ielem = 0; ielem < numElemGlueZoneFine; ielem++) {
        
        VecDouble x1(nElNodes), x2(nElNodes);
        VecInt connec = elementsFine_[elementsGlueZoneFine_[ielem]] -> getConnectivity();
        
        for (int i = 0; i < nElNodes; i++){
            VecDouble x = (*nodesFine_)[connec[i]] -> getCoordinates();
            VecDouble xp = (*nodesFine_)[connec[i]] -> getPreviousCoordinates();
            
            x1[i] = alpha_f * x[0] + (1. - alpha_f) * xp[0];
            x2[i] = alpha_f * x[1] + (1. - alpha_f) * xp[1];
        };

        // std::cout << "XX1 " << x2 << " " << x22 << " " << x222 << std::endl;

        int numberIntPoints = elementsFine_[elementsGlueZoneFine_[ielem]] -> 
            getNumberOfIntegrationPoints();

        for (int i = 0; i < numberIntPoints; i++){

            VecDouble x(DIM);

            x[0] = quad.interpolateQuadraticVariable(x1,i);
            x[1] = quad.interpolateQuadraticVariable(x2,i);

            int elCorr = 0;
            VecDouble xsiCorr(DIM);
        
            searchNodeCorrespondence(x,(*nodesCoarse_),elementsCoarse_,
                                     elementsCoarse_.size(),elCorr,xsiCorr,elementsFine_[elementsGlueZoneFine_[ielem]] -> getIntegPointCorrespondenceElement(i));
                
            elementsFine_[elementsGlueZoneFine_[ielem]] -> setIntegrationPointCorrespondence(i,elCorr,xsiCorr);

            if (elCorr == 150000) {
                std::cout << "PROBLEM IN NODAL CORRESPONDENCE OF FINE NODE " 
                       << std::endl;
                std::cout << "INODE " << nodesGlueZoneFine_[ielem] << std::endl;
                std::cout << "elsglue " << ielem << " " <<  elementsGlueZoneFine_[ielem] << std::endl;
                       std::cout << elementsFine_[elementsGlueZoneFine_[ielem]] -> getIntegPointCorrespondenceElement(i) << std::endl;
                       std::cout << x1 << " " << x2 << std::endl;
                       std::cout << x[0] << " " << x[1] << std::endl;
                       std::cout << "CORRESP " << elCorr << std::endl
                                 << xsiCorr[0] << " " << xsiCorr[1]  << std::endl;

                exit(0);
            };
        };

    };
};

//------------------------------------------------------------------------------
//--------------------SETS THE COUPLING ZONE IN COARSE MODEL--------------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Arlequin<DIM,DEG>::setSignaledDistance(){

    int bconnec[nBdNodes];
    double dist;

    for (int inode = 0 ; inode < numNodesFine; inode++){ 
        (*nodesFine_)[inode] -> clearInnerNormal();
        (*nodesFine_)[inode] -> setDistFunction(0.0);
    };
    for (int inode = 0 ; inode < numNodesCoarse; inode++){ 
        (*nodesCoarse_)[inode] -> setDistFunction(0.0);
    };
    //approximate normal calculation

    for (int i = 0; i < numBoundElemFine; i++){
        if (boundaryFine_[i]->getConstrain(0) == 2){
            VecInt connec = elementsFine_[boundaryFine_[i] -> getElement()] -> getConnectivity();
            
            elementsFine_[boundaryFine_[i] -> getElement()] -> getBoundaryNodes(bconnec);
            
            //Loop over the 1D element segments
            for (int iSeg = 0; iSeg < DEG; iSeg++){
                int no1,no2;
                if (iSeg == 0){
                    no1 = bconnec[0];
                    if (DEG == 1) {
                        no2 = bconnec[1];
                    } else {
                        no2 = bconnec[2];
                    }
                } else {
                    no1 = bconnec[iSeg+1];
                    if (DEG == 2 || iSeg == 2){
                        no2 = bconnec[1];
                    } else {
                        no2 = bconnec[3];
                    }
                }
                
                VecDouble x1 = (*nodesFine_)[no1] -> getCoordinates();
                VecDouble x2 = (*nodesFine_)[no2] -> getCoordinates();

                double sLength = sqrt((x2[1] - x1[1]) * (x2[1] - x1[1]) +
                                      (x1[0] - x2[0]) * (x1[0] - x2[0]));

                VecDouble n(2);
                n[0] = (x2[1] - x1[1]) / sLength;
                n[1] = (x1[0] - x2[0]) / sLength;

                (*nodesFine_)[no1] -> setInnerNormal(n);
                (*nodesFine_)[no2] -> setInnerNormal(n);
            }
        };
    };


    //Fine mesh nodes
    for (int ino = 0; ino < numNodesFine; ino++){
        VecDouble x = (*nodesFine_)[ino] -> getCoordinates();
        dist=10000000000000000000000000000.;
        
        for (int i = 0; i < numBoundElemFine; i++){
            if (boundaryFine_[i] -> getConstrain(0) == 2){
                VecInt connec = elementsFine_[boundaryFine_[i] -> getElement()] -> getConnectivity();
                
                elementsFine_[boundaryFine_[i] -> getElement()] -> getBoundaryNodes(bconnec);

                

                for (int iSeg = 0; iSeg < DEG ; ++iSeg){
                    
                
                    //first segment
                    int no1,no2;
                    if (iSeg == 0){
                        no1 = bconnec[0];
                        if (DEG == 1) {
                            no2 = bconnec[1];
                        } else {
                            no2 = bconnec[2];
                        }
                    } else {
                        no1 = bconnec[iSeg+1];
                        if (DEG == 2 || iSeg == 2){
                            no2 = bconnec[1];
                        } else {
                            no2 = bconnec[3];
                        }
                    }
                  
                    VecDouble x1 = (*nodesFine_)[no1] -> getCoordinates();
                    VecDouble x2 = (*nodesFine_)[no2] -> getCoordinates();
                    
                    double aux0 =  sqrt((x2[1] - x1[1]) * (x2[1] - x1[1]) +
                                        (x2[0] - x1[0]) * (x2[0] - x1[0]));
                    double aux1 = ((x[0] - x1[0]) * (x2[0] - x1[0])+
                                   (x[1] - x1[1]) * (x2[1] - x1[1])) / aux0;
                    double dist2 =-((x2[1] - x1[1]) * x[0] - 
                                    (x2[0] - x1[0]) * x[1] +
                                    x2[0] * x1[1] - x2[1] * x1[0]) / aux0;
                    
                    if (aux1 > aux0){
                        dist2 = sqrt((x2[1] - x[1]) * (x2[1] - x[1]) +
                                     (x2[0] - x[0]) * (x2[0] - x[0]));
                        //find signal
                        //side normal vector
                        VecDouble n = (*nodesFine_)[no2] -> getInnerNormal();
                        double test[2];

                        test[0] = x[0] - x2[0];
                        test[1] = x[1] - x2[1];
                        double signaltest = n[0]*test[0] + n[1]*test[1];
                        double signal = -1.;
                        
                        if (signaltest <= -0.001)signal = 1.;
                        
                        dist2 *= signal;
                    };

                    if (aux1 < 0.){
                        dist2 = sqrt((x[1] - x1[1]) * (x[1] - x1[1]) +
                                     (x[0] - x1[0]) * (x[0] - x1[0]));
                        //find signal
                        //side normal vector
                        VecDouble n = (*nodesFine_)[no1] -> getInnerNormal();
                        double test[2];

                        test[0] = x[0] - x1[0];
                        test[1] = x[1] - x1[1];
                        double signaltest = n[0]*test[0] + n[1]*test[1];
                        double signal = -1.;
                        
                        if (signaltest <= -0.001) signal = 1.;
                        
                        dist2 *= signal;
                    };
                    
                    if (fabs(dist2) < fabs(dist)) dist = dist2;
                }
            }; //if bf is the glue boundary
        }; //
        // if(dist < 0) dist = 0;
        // dist = x[0];
        (*nodesFine_)[ino] -> setDistFunction(dist);
     };

    //Coarse mesh
     for (int ino = 0; ino < numNodesCoarse; ino++){
        VecDouble x = (*nodesCoarse_)[ino] -> getCoordinates();
        dist=10000000000000000000000000000.;
        
        for (int i = 0; i < numBoundElemFine; i++){
            if (boundaryFine_[i] -> getConstrain(0) == 2){

                VecInt connec = elementsFine_[boundaryFine_[i] -> getElement()] -> getConnectivity();
                elementsFine_[boundaryFine_[i] -> getElement()] -> getBoundaryNodes(bconnec);

                for (int iSeg = 0; iSeg < DEG ; ++iSeg){
                    //first segment
                    int no1,no2;
                    if (iSeg == 0){
                        no1 = bconnec[0];
                        if (DEG == 1) {
                            no2 = bconnec[1];
                        } else {
                            no2 = bconnec[2];
                        }
                    } else {
                        no1 = bconnec[iSeg+1];
                        if (DEG == 2 || iSeg == 2){
                            no2 = bconnec[1];
                        } else {
                            no2 = bconnec[3];
                        }
                    }
                    // std::cout<<no1<<" nos "<<no2<<std::endl;
                  
                    VecDouble x1 = (*nodesFine_)[no1] -> getCoordinates();
                    VecDouble x2 = (*nodesFine_)[no2] -> getCoordinates();
                    
                    double aux0 =  sqrt((x2[1] - x1[1]) * (x2[1] - x1[1]) +
                                        (x2[0] - x1[0]) * (x2[0] - x1[0]));
                    double aux1 = ((x[0] - x1[0]) * (x2[0] - x1[0])+
                                   (x[1] - x1[1]) * (x2[1] - x1[1])) / aux0;
                    double dist2 =-((x2[1] - x1[1]) * x[0] - 
                                    (x2[0] - x1[0]) * x[1] +
                                    x2[0] * x1[1] - x2[1] * x1[0]) / aux0;
                    
                    if (aux1 > aux0){
                        dist2 = sqrt((x2[1] - x[1]) * (x2[1] - x[1]) +
                                     (x2[0] - x[0]) * (x2[0] - x[0]));
                        //find signal
                        //side normal vector
                        VecDouble n = (*nodesFine_)[no2] -> getInnerNormal();
                        
                        double test[2];
                        test[0] = x[0] - x2[0];
                        test[1] = x[1] - x2[1];
                        double signaltest = n[0]*test[0] + n[1]*test[1];
                        double signal = -1.;
                        
                        if (signaltest <= -0.001)signal = 1.;
                        
                        dist2 *= signal;
                    };

                    if (aux1 < 0.){
                        dist2 = sqrt((x[1] - x1[1]) * (x[1] - x1[1]) +
                                     (x[0] - x1[0]) * (x[0] - x1[0]));
                        //find signal
                        //side normal vector
                        VecDouble n = (*nodesFine_)[no1] -> getInnerNormal();
                        
                        double test[2];
                        test[0] = x[0] - x1[0];
                        test[1] = x[1] - x1[1];
                        double signaltest = n[0]*test[0] + n[1]*test[1];
                        double signal = -1.;
                        
                        if (signaltest <= -0.001) signal = 1.;
                        
                        dist2 *= signal;
                    };
                    
                    if (fabs(dist2) < fabs(dist)) dist = dist2;
                }
            }; //if bf is the glue boundary
        }; //
    
        if (fabs((*nodesCoarse_)[ino] -> getDistFunction()) < 1.e-2){
            (*nodesCoarse_)[ino] -> setDistFunction(dist); 
        };
        // dist = x[0] - 0.5;
        (*nodesCoarse_)[ino] -> setDistFunction(dist); 
     };


};

//------------------------------------------------------------------------------
//--------------------SETS THE COUPLING ZONE IN COARSE MODEL--------------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Arlequin<DIM,DEG>::setCouplingZone(){

    // double dist;
    int flag;
    int nodesCZ[numNodesFine];
    int nodesCZ2[numNodesCoarse];

    double lim1 = 0.06251;
    double lim2 = 0.93749;
    double tick = 0.01249;

    for (int i = 0; i < numNodesFine; i++) nodesCZ[i] = 0;    

    elementsGlueZoneFine_.reserve(numElemFine / 3);
    nodesGlueZoneFine_.reserve(numNodesFine / 3);
    glueZoneFine_.reserve(numNodesFine / 3);
    int index = 0; 

    //Defines a criterion to select the elements that are in the glue zone
    for (int jel = 0; jel < numElemFine; jel++){
        
        VecInt connec = elementsFine_[jel] -> getConnectivity();
        flag = 0;

        for (int ino = 0; ino < nElNodes; ino++){
            VecDouble x = (*nodesFine_)[connec[ino]] -> getCoordinates();
            double dist = (*nodesFine_)[connec[ino]] -> getDistFunction();
            //  std::cout << "DIST " << dist << std::endl;
            if (dist <= fineModel.glueZoneThickness + 0.001) flag++;
            
        };

        if (flag == nElNodes) {
            elementsGlueZoneFine_.push_back(jel);
            elementsFine_[jel] -> setGlueZone();

            GlueZone *el = new GlueZone(index++,jel);
            glueZoneFine_.push_back(el);
            
        };        
    };


    //Defines which nodes are in the glue zone
    numElemGlueZoneFine = elementsGlueZoneFine_.size();
    for (int i = 0; i < numElemGlueZoneFine; i++){
        VecInt connec = elementsFine_[elementsGlueZoneFine_[i]] -> getConnectivity();

        for (int ino = 0; ino < nElNodes; ino++) nodesCZ[connec[ino]]++;
    };

    //Compute number of nodes in the glue zone
    numNodesGlueZoneFine = 0;
    for (int i = 0; i < numNodesFine; i++){
        if(nodesCZ[i] > 0) {
            numNodesGlueZoneFine++;
            nodesGlueZoneFine_.push_back(i);
        };
    };

    if (rank == 0) std::cout << "GLUE ZONE - Number of Nodes = " 
                             << numNodesGlueZoneFine 
                             << " - Number of Elements = " 
                             << elementsGlueZoneFine_.size() << std::endl;
    
    for (int i = 0; i < numNodesGlueZoneFine; i++){
        VecDouble x = (*nodesFine_)[nodesGlueZoneFine_[i]] -> getCoordinates();
        
        Nodes *no = new Nodes(x,i);
        nodesLagrangeFine_.push_back(no);
    };

    for (int i = 0; i < numElemGlueZoneFine; i++){
        VecInt connecAux(nElNodes);

        VecInt connec = elementsFine_[elementsGlueZoneFine_[i]] -> getConnectivity();
        
        for (int ino = 0; ino < numNodesGlueZoneFine; ino++)
            for (int k = 0; k < nElNodes; k++)
                if (nodesGlueZoneFine_[ino] == connec[k]) connecAux[k] = ino;
        
        
        glueZoneFine_[i] -> setConnectivity(connecAux);
        // glueZoneFine_[i] -> setNodes(nodesLagrangeFine_);

    };

    // Glue Zone in coarse mesh
    for (int i = 0; i < numNodesCoarse; i++) nodesCZ2[i] = 0;    

    elementsGlueZoneCoarse_.reserve(numElemCoarse / 3);
    nodesGlueZoneCoarse_.reserve(numNodesCoarse / 3);

    //Defines a criterion to select the elements that are in the glue zone
    for (int jel = 0; jel < numElemCoarse; jel++){
        
        VecInt connec = elementsCoarse_[jel] -> getConnectivity();
        flag = 0;

        for (int ino = 0; ino < nElNodes; ino++){
            VecDouble x = (*nodesCoarse_)[connec[ino]] -> getCoordinates();

            if ((x[0] < lim1) || (x[0] > lim2) || 
                (x[1] < lim1) || (x[1] > lim2)){

            }else{
                if ((x[0] > lim1 + tick) && (x[0] < lim2 - tick) &&
                    (x[1] > lim1 + tick) && (x[1] < lim2 - tick)){

                }else{
                    flag = 1;
                    break;
                };
            };
        };
        if (flag > 0) {
            elementsGlueZoneCoarse_.push_back(jel);
            elementsCoarse_[jel] -> setGlueZone();

            GlueZone *el = new GlueZone(index++,jel);
            glueZoneCoarse_.push_back(el);
            
            for (int i=0; i < elementsCoarse_[jel] -> getNumberOfIntegrationPoints(); i++){
                VecDouble x = elementsCoarse_[jel] -> getIntegPointCoordinatesValue(i);

                if ((x[0] < lim1) || (x[0] > lim2) || 
                    (x[1] < lim1) || (x[1] > lim2)){

                }else{
                    if ((x[0] > lim1 + tick) && (x[0] < lim2 - tick) &&
                        (x[1] > lim1 + tick) && (x[1] < lim2 - tick)){

                    }else{
                    elementsCoarse_[jel] -> setIntegPointInGlueZone(i);
                    };
                };
            };
        };        
    };

    //Defines which nodes are in the glue zone
    numElemGlueZoneCoarse = elementsGlueZoneCoarse_.size();
    for (int i = 0; i < numElemGlueZoneCoarse; i++){
        VecInt connec = elementsCoarse_[elementsGlueZoneCoarse_[i]] -> getConnectivity();

        for (int ino = 0; ino < nElNodes; ino++){
            nodesCZ2[connec[ino]] += 1;
        };
    };

    //Compute number of nodes in the glue zone
    numNodesGlueZoneCoarse = 0;
    for (int i = 0; i < numNodesCoarse; i++){
        if(nodesCZ2[i] > 0) {
            numNodesGlueZoneCoarse += 1;
            nodesGlueZoneCoarse_.push_back(i);
        };
    };

    for (int i = 0; i < numNodesGlueZoneCoarse; i++){
        VecDouble x = (*nodesCoarse_)[nodesGlueZoneCoarse_[i]] -> getCoordinates();
        
        Nodes *no = new Nodes(x,i);
        nodesLagrangeCoarse_.push_back(no);
    };

    for (int i = 0; i < numElemGlueZoneCoarse; i++){
        VecInt connecAux(nElNodes);

        VecInt connec = elementsCoarse_[elementsGlueZoneCoarse_[i]] -> getConnectivity();
        
        for (int ino = 0; ino < numNodesGlueZoneCoarse; ino++)
            for (int k = 0; k < nElNodes; k++)
                if (nodesGlueZoneCoarse_[ino] == connec[k]) connecAux[k] = ino;
            
        
        glueZoneCoarse_[i] -> setConnectivity(connecAux);
    };
     
};


//------------------------------------------------------------------------------
//---------------------------SETS THE WEIGHT FUNCTION---------------------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Arlequin<DIM,DEG>::setWeightFunction(double val){
    
    double wFuncValue;

    double epsilon = coarseModel.arlequinEpsilon;
    double lambda = fineModel.glueZoneThickness*1.01;
 
    for (int i = 0; i < numNodesCoarse; i++){
        
        double r = (*nodesCoarse_)[i] -> getDistFunction();
            
        if (r < 0){
            wFuncValue = 1.;
            //if (fabs(r) < 1.e-5) wFuncValue = 0.5;
        } else {
            if (r >= lambda){
                wFuncValue = epsilon;
            } else {
                wFuncValue = 1. - (1. - epsilon) / lambda * r;
                //wFuncValue = 0.5;
                //wFuncValue = epsilon;
                // wFuncValue = 1+3.*(epsilon-1.)/(lambda*lambda) * r * r
                //     -2.*(epsilon-1.)/(lambda*lambda*lambda) * r * r * r;
                
                if (wFuncValue < epsilon) wFuncValue = epsilon;
            };
        };  
              
        (*nodesCoarse_)[i] -> setWeightFunction(wFuncValue);
    };         

    for (int jel = 0; jel < numElemCoarse; jel++){
        elementsCoarse_[jel] -> setIntegPointWeightFunction();        
    };


    for (int i=0; i<numNodesFine; i++){

        double r = (*nodesFine_)[i] -> getDistFunction();
        
        if (r >= lambda){
            wFuncValue = 1. - epsilon;
        } else {
            wFuncValue = (1. - epsilon) / lambda * r;
            //wFuncValue = 0.5;
            //wFuncValue = 1. - epsilon;
            // wFuncValue = -3.*(epsilon-1.)/(lambda*lambda) * r * r
            //         +2.*(epsilon-1.)/(lambda*lambda*lambda) * r * r * r;

            if (wFuncValue > (1. - epsilon)) wFuncValue = 1. - epsilon;

        };  
        
        // wFuncValue = weightFunctionFineValue(r,epsilon);

        (*nodesFine_)[i] -> setWeightFunction(wFuncValue);
    };
    for (int jel = 0; jel < numElemFine; jel++){
        elementsFine_[jel] -> setIntegPointWeightFunction();        
    };  
     
    return;
};


//------------------------------------------------------------------------------
//-----------------------------PRINT COARSE RESULTS-----------------------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Arlequin<DIM,DEG>::printResultsCoarse(int step) {

    // PRINT COARSE MODEL RESULTS
    std::string result;
    std::ostringstream convert;

    convert << step+100000;
    result = convert.str();

    std::string s = "saidaVelCoarse"+result+".xdmf";
    std::fstream output_v(s.c_str(), std::ios_base::out);

    std::string s1 = "resultCoarse"+result+".h5";
    std::fstream filename(s1.c_str(), std::ios_base::out);

    std::string s2 = "geometryCoarse.h5";
    if (step == 0){
        std::fstream filename(s2.c_str(), std::ios_base::out);
    }

    //Auxiliary vectors to write HDF5 file
    double *pointVector;
    int *connec2;
    double *pointScalar;
    int *intCellScalar;
    double *douCellScalar;

    pointVector = new double[3*numNodesCoarse]();
    connec2 = new int[nElNodes*numElemCoarse]();
    pointScalar = new double[numNodesCoarse]();
    intCellScalar = new int[numElemCoarse]();
    douCellScalar = new double[numElemCoarse]();

    for (int i = 0; i < numNodesCoarse; i++) {
        VecDouble x = (*nodesCoarse_)[i] -> getCoordinates();
        pointVector[3*i  ] = x[0];
        pointVector[3*i+1] = x[1];
        pointVector[3*i+2] = 0.0;
    }
    for (int iElem = 0; iElem < numElemCoarse; iElem++){
        VecInt con = elementsCoarse_[iElem] -> getConnectivity();

        if (DIM == 2){
            switch (DEG)
            {
            case 1:
                for (int i = 0; i < nElNodes; ++i) connec2[nElNodes*iElem+i] = con[i];
                break;
            case 2:
                {
                    connec2[nElNodes*iElem+0] = con[4];
                    connec2[nElNodes*iElem+1] = con[1];
                    connec2[nElNodes*iElem+2] = con[3];
                    connec2[nElNodes*iElem+3] = con[0];
                    connec2[nElNodes*iElem+4] = con[5];
                    connec2[nElNodes*iElem+5] = con[2];
                    break;
                }
            case 3:
                {
                    connec2[nElNodes*iElem+0] = con[6];
                    connec2[nElNodes*iElem+1] = con[9];
                    connec2[nElNodes*iElem+2] = con[5];
                    connec2[nElNodes*iElem+3] = con[1];
                    connec2[nElNodes*iElem+4] = con[4];
                    connec2[nElNodes*iElem+5] = con[3];
                    connec2[nElNodes*iElem+6] = con[0];
                    connec2[nElNodes*iElem+7] = con[8];
                    connec2[nElNodes*iElem+8] = con[7];
                    connec2[nElNodes*iElem+9] = con[2];
                    break;
                }
            
            default:
                PanicButton();
                break;
            }
        } else {
            for (int i = 0; i < nElNodes; ++i) connec2[nElNodes*iElem+i] = con[i];
        }
        
        
    }
   
    hid_t file, file2; 
    hid_t dataset; 
    hid_t dataspace;

    herr_t status;
    hsize_t xdim = numNodesCoarse;
    hsize_t eldim = numElemCoarse;
    hsize_t elnod = nElNodes;
    hsize_t pointVectorDims[2] = { xdim, 3 };
    hsize_t pointScalarDims[2] = { xdim, 1 };
    hsize_t pointCellScalarDims[2] = { eldim, 1 };
    hsize_t connec2Dims[2] = { eldim, elnod };

    std::string elType;
    if (DIM == 2) elType = "TRIANGLE";
    if (DIM == 3) elType = "TETRAHEDRON";

    //Create HDF5 file
    file = H5Fcreate(s1.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (step == 0) file2 = H5Fcreate(s2.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    //Write xdmf file
    output_v << "<?xml version=\"1.0\"?>" << std::endl
             << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>" << std::endl
             << "<Xdmf Version=\"2.0\" xmlns:xi=\"http://www.w3.org/2001/XInclude\" >" << std::endl
             << "<Domain>" << std::endl
             << "  <Grid>"  << std::endl
             << "    <Topology TopologyType=\"" << elType << "\" NumberOfElements=\"" << numElemCoarse << "\" >" << std::endl
             << "      <DataItem Format=    \"HDF\" NumberType=\"int\" Dimensions=\"" << numElemCoarse << " " << nElNodes << "\" >" << std::endl;
    
    //Connectivity
    output_v << "        " << s2 << ":/connec" << std::endl;
    //Start connectivity
    if (step == 0){
        dataspace = H5Screate_simple(2, connec2Dims, NULL);
        dataset = H5Dcreate2(file2, "/connec", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &connec2[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
    }
    //End connectivity

    output_v << "      </DataItem>" << std::endl
             << "    </Topology>" << std::endl
             << "    <Geometry GeometryType=\"XYZ\">" << std::endl
             << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< numNodesCoarse <<" 3\">" << std::endl;

    //Coordinates   
    output_v << "        " << s2 << ":/coords" << std::endl;
    //Start coordinates
    if (step == 0){
        dataspace = H5Screate_simple(2, pointVectorDims, NULL);
        dataset = H5Dcreate2(file2, "/coords", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointVector[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
    }
    //End coordinates
    output_v << "      </DataItem>" << std::endl
             << "    </Geometry>" << std::endl;

    //LISTS
    //Velocity
    if (coarseModel.printVelocity){
        output_v << "    <Attribute Name=\"Velocity\" Center=\"Node\" AttributeType=\"Vector\" >" << std::endl
                 << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< numNodesCoarse <<" 3\">" << std::endl;  
        output_v << "        " << s1 << ":/velocity" << std::endl;     
        //Start Velocity
        for (int i = 0; i < numNodesCoarse; i++) {
            pointVector[3*i  ] = (*nodesCoarse_)[i] -> getVelocity(0);
            pointVector[3*i+1] = (*nodesCoarse_)[i] -> getVelocity(1);
            pointVector[3*i+2] = 0.0;
        };
        dataspace = H5Screate_simple(2, pointVectorDims, NULL);
        dataset = H5Dcreate2(file, "/velocity", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointVector[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
        //End Velocity
        output_v << "      </DataItem>" << std::endl
                 << "    </Attribute>" << std::endl;
    }

    //Acceleration
    if (coarseModel.printAcceleration){
        output_v << "    <Attribute Name=\"Acceleration\" Center=\"Node\" AttributeType=\"Vector\" >" << std::endl
                 << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< numNodesCoarse <<" 3\">" << std::endl;  
        output_v << "        " << s1 << ":/acceleration" << std::endl;     
        //Start Acceleration
        for (int i = 0; i < numNodesCoarse; i++) {
            pointVector[3*i  ] = (*nodesCoarse_)[i] -> getAcceleration(0);
            pointVector[3*i+1] = (*nodesCoarse_)[i] -> getAcceleration(1);
            pointVector[3*i+2] = 0.0;
        };
        dataspace = H5Screate_simple(2, pointVectorDims, NULL);
        dataset = H5Dcreate2(file, "/acceleration", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointVector[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
        //End Acceleration
        output_v << "      </DataItem>" << std::endl
                 << "    </Attribute>" << std::endl;
    }

    //Real Velocity
    if (coarseModel.printRealVelocity){
        output_v << "    <Attribute Name=\"Real Velocity\" Center=\"Node\" AttributeType=\"Vector\" >" << std::endl
                 << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< numNodesCoarse <<" 3\">" << std::endl;  
        output_v << "        " << s1 << ":/realVelocity" << std::endl;     
        //Start Real Velocity
        for (int i = 0; i < numNodesCoarse; i++) {
            pointVector[3*i  ] = (*nodesCoarse_)[i] -> getVelocityArlequin(0);
            pointVector[3*i+1] = (*nodesCoarse_)[i] -> getVelocityArlequin(1);
            pointVector[3*i+2] = 0.0;
        };
        dataspace = H5Screate_simple(2, pointVectorDims, NULL);
        dataset = H5Dcreate2(file, "/realVelocity", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointVector[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
        //End Real Velocity
        output_v << "      </DataItem>" << std::endl
                 << "    </Attribute>" << std::endl;
    }

    //Lagrange Multipliers
    if (coarseModel.printLagrangeMultipliers){
        output_v << "    <Attribute Name=\"Lagrange Multipliers\" Center=\"Node\" AttributeType=\"Vector\" >" << std::endl
                 << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< numNodesCoarse <<" 3\">" << std::endl;  
        output_v << "        " << s1 << ":/lagrangeMultiplers" << std::endl;     
        //Start Lagrange Multipliers
        for (int i = 0; i < numNodesCoarse; i++) {
            pointVector[3*i  ] = (*nodesCoarse_)[i] -> getLagrangeMultiplier(0);
            pointVector[3*i+1] = (*nodesCoarse_)[i] -> getLagrangeMultiplier(1);
            pointVector[3*i+2] = 0.0;
        };
        dataspace = H5Screate_simple(2, pointVectorDims, NULL);
        dataset = H5Dcreate2(file, "/lagrangeMultiplers", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointVector[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
        //End Real Velocity
        output_v << "      </DataItem>" << std::endl
                 << "    </Attribute>" << std::endl;
    }

    //Dist Function
    if (coarseModel.printDistFunction){
        output_v << "    <Attribute Name=\"Dist Function\" Center=\"Node\" AttributeType=\"Scalar\" >" << std::endl
                 << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\"1 "<< numNodesCoarse <<"\">" << std::endl;  
        output_v << "        " << s1 << ":/distfunction" << std::endl;     
        //Start Dist Function
        for (int i = 0; i < numNodesCoarse; i++) pointScalar[i] = (*nodesCoarse_)[i] -> getDistFunction();
        dataspace = H5Screate_simple(2, pointScalarDims, NULL);
        dataset = H5Dcreate2(file, "/distfunction", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointScalar[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
        //End Dist Function
        output_v << "      </DataItem>" << std::endl
                 << "    </Attribute>" << std::endl;
    }

    //Weight Function
    if (coarseModel.printEnergyWeightFunction){
        output_v << "    <Attribute Name=\"Weight Function\" Center=\"Node\" AttributeType=\"Scalar\" >" << std::endl
                 << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\"1 "<< numNodesCoarse <<"\">" << std::endl;  
        output_v << "        " << s1 << ":/weightFunction" << std::endl;     
        //Start Weigth Function
        for (int i = 0; i < numNodesCoarse; i++) pointScalar[i] = (*nodesCoarse_)[i] -> getWeightFunction();
        dataspace = H5Screate_simple(2, pointScalarDims, NULL);
        dataset = H5Dcreate2(file, "/weightFunction", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointScalar[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
        //End Weight Function
        output_v << "      </DataItem>" << std::endl
                 << "    </Attribute>" << std::endl;
    }

    //Pressure
    if (coarseModel.printPressure){
        output_v << "    <Attribute Name=\"Pressure\" Center=\"Node\" AttributeType=\"Scalar\" >" << std::endl
                 << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< numNodesCoarse <<" 1\">" << std::endl;  
        output_v << "        " << s1 << ":/pressure" << std::endl;     
        //Start Pressure
        for (int i = 0; i < numNodesCoarse; i++) pointScalar[i] = (*nodesCoarse_)[i] -> getPressure();
        dataspace = H5Screate_simple(2, pointScalarDims, NULL);
        dataset = H5Dcreate2(file, "/pressure", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointScalar[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
        //End Pressure
        output_v << "      </DataItem>" << std::endl
                 << "    </Attribute>" << std::endl;
    }

    //Real Pressure
    if (coarseModel.printRealPressure){
        output_v << "    <Attribute Name=\"Real Pressure\" Center=\"Node\" AttributeType=\"Scalar\" >" << std::endl
                 << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< numNodesCoarse <<" 1\">" << std::endl;  
        output_v << "        " << s1 << ":/realPressure" << std::endl;     
        //Start Real Pressure
        for (int i = 0; i < numNodesCoarse; i++) pointScalar[i] = (*nodesCoarse_)[i] -> getPressureArlequin();
        dataspace = H5Screate_simple(2, pointScalarDims, NULL);
        dataset = H5Dcreate2(file, "/realPressure", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointScalar[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
        //End Real Pressure
        output_v << "      </DataItem>" << std::endl
                 << "    </Attribute>" << std::endl;
    }

    //Process
    if (coarseModel.printProcess){
        output_v << "    <Attribute Name=\"Process\" Center=\"Cell\" AttributeType=\"Scalar\" >" << std::endl
                 << "      <DataItem Format=\"HDF\" NumberType=\"int\" Dimensions=\""<< numElemCoarse <<" 1\">" << std::endl;  
        output_v << "        " << s1 << ":/process" << std::endl;
        //Start Process
        dataspace = H5Screate_simple(2, pointCellScalarDims, NULL);
        dataset = H5Dcreate2(file, "/process", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &domDecompCoarse.first[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
        //End Process
        output_v << "      </DataItem>" << std::endl
                 << "    </Attribute>" << std::endl;
    }

    //WeightFunction2
    if (coarseModel.printEnergyWeightFunction){
        output_v << "    <Attribute Name=\"Weight Function\" Center=\"Cell\" AttributeType=\"Scalar\" >" << std::endl
                 << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< numElemCoarse <<" 1\">" << std::endl;  
        output_v << "        " << s1 << ":/weightFunction2" << std::endl;     
        //Start Weight Function
        for (int i = 0; i < numElemCoarse; i++) douCellScalar[i] = elementsCoarse_[i] -> getIntegPointWeightFunction(0);
        dataspace = H5Screate_simple(2, pointCellScalarDims, NULL);
        dataset = H5Dcreate2(file, "/weightFunction2", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &douCellScalar[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
        //End Weight Function
        output_v << "      </DataItem>" << std::endl
                 << "    </Attribute>" << std::endl;
    }

    //Lines
    if (coarseModel.printLines){
        output_v << "    <Attribute Name=\"Lines\" Center=\"Cell\" AttributeType=\"Scalar\" >" << std::endl
                 << "      <DataItem Format=\"HDF\" NumberType=\"int\" Dimensions=\""<< numElemCoarse <<" 1\">" << std::endl;  
        output_v << "        " << s1 << ":/lines" << std::endl;     
        //Start Lines
        for (int i = 0; i < numElemCoarse; i++){
            int res = 0;
            for (int j = 0; j < numBoundElemCoarse; j++)
               if (boundaryCoarse_[j] -> getElement() == i) res = boundaryCoarse_[j] -> getBoundaryGroup();
            intCellScalar[i] = res;
        };
        dataspace = H5Screate_simple(2, pointCellScalarDims, NULL);
        dataset = H5Dcreate2(file, "/lines", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &intCellScalar[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
        //End Lines
        output_v << "      </DataItem>" << std::endl
                 << "    </Attribute>" << std::endl;
    }

    // Jacobian
    if (coarseModel.printJacobian){
        output_v << "    <Attribute Name=\"Jacobian\" Center=\"Cell\" AttributeType=\"Scalar\" >" << std::endl
                 << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< numElemCoarse <<" 1\">" << std::endl;  
        output_v << "        " << s1 << ":/jacobian" << std::endl;     
        //Start Weight Function
        for (int i = 0; i < numElemCoarse; i++) douCellScalar[i] = elementsCoarse_[i] -> getJacobian();
        dataspace = H5Screate_simple(2, pointCellScalarDims, NULL);
        dataset = H5Dcreate2(file, "/jacobian", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &douCellScalar[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
        //End Weight Function
        output_v << "      </DataItem>" << std::endl
                 << "    </Attribute>" << std::endl;
    }

    // END FILE
    output_v << "  </Grid>" << std::endl
             << "</Domain>" << std::endl
             << "</Xdmf>" << std::endl;

    //Delete auxiliary vectors
    delete [] pointVector;
    delete [] connec2;
    delete [] pointScalar;
    delete [] intCellScalar;
    delete [] douCellScalar;

    //End HDF5 file
    status = H5Fclose(file);
    if (step == 0) status = H5Fclose(file2);

    return;
};

//------------------------------------------------------------------------------
//------------------------------PRINT FINE RESULTS------------------------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Arlequin<DIM,DEG>::printResultsFine(int step) {

    //PRINT FINE MODEL RESULTS
    std::string result;
    std::ostringstream convert;

    convert << step+100000;
    result = convert.str();

    std::string s = "saidaVelFine"+result+".xdmf";
    std::fstream output_v(s.c_str(), std::ios_base::out);

    std::string s1 = "resultFine"+result+".h5";
    std::fstream filename(s1.c_str(), std::ios_base::out);

    //Auxiliary vectors to write HDF5 file
    double *pointVector;
    int *connec2;
    double *pointScalar;
    int *intCellScalar;
    double *douCellScalar;

    pointVector = new double[3*numNodesFine]();
    connec2 = new int[nElNodes*numElemFine]();
    pointScalar = new double[numNodesFine]();
    intCellScalar = new int[numElemFine]();
    douCellScalar = new double[numElemFine]();

    for (int i = 0; i < numNodesFine; i++) {
        VecDouble x = (*nodesFine_)[i] -> getCoordinates();
        pointVector[3*i  ] = x[0];
        pointVector[3*i+1] = x[1];
        pointVector[3*i+2] = 0.0;
    }

    for (int iElem = 0; iElem < numElemFine; iElem++){
        VecInt con = elementsFine_[iElem] -> getConnectivity();
        if (DIM == 2){
            switch (DEG)
            {
            case 1:
                for (int i = 0; i < nElNodes; ++i) connec2[nElNodes*iElem+i] = con[i];
                break;
            case 2:
                {
                    connec2[nElNodes*iElem+0] = con[4];
                    connec2[nElNodes*iElem+1] = con[1];
                    connec2[nElNodes*iElem+2] = con[3];
                    connec2[nElNodes*iElem+3] = con[0];
                    connec2[nElNodes*iElem+4] = con[5];
                    connec2[nElNodes*iElem+5] = con[2];
                    break;
                }
            case 3:
                {
                    connec2[nElNodes*iElem+0] = con[6];
                    connec2[nElNodes*iElem+1] = con[9];
                    connec2[nElNodes*iElem+2] = con[5];
                    connec2[nElNodes*iElem+3] = con[1];
                    connec2[nElNodes*iElem+4] = con[4];
                    connec2[nElNodes*iElem+5] = con[3];
                    connec2[nElNodes*iElem+6] = con[0];
                    connec2[nElNodes*iElem+7] = con[8];
                    connec2[nElNodes*iElem+8] = con[7];
                    connec2[nElNodes*iElem+9] = con[2];
                    break;
                }      
            
            default:
                PanicButton();
                break;
            }
        } else {
            for (int i = 0; i < nElNodes; ++i) connec2[nElNodes*iElem+i] = con[i];
        }
    } 

    hid_t file; 
    hid_t dataset; 
    hid_t dataspace;

    herr_t status;
    hsize_t xdim = numNodesFine;
    hsize_t eldim = numElemFine;
    hsize_t elnod = nElNodes;
    hsize_t pointVectorDims[2] = { xdim, 3 };
    hsize_t pointScalarDims[2] = { xdim, 1 };
    hsize_t pointCellScalarDims[2] = { eldim, 1 };
    hsize_t connec2Dims[2] = { eldim, elnod };

    std::string elType;
    if (DIM == 2) elType = "TRIANGLE";
    if (DIM == 3) elType = "TETRAHEDRON";

    //Create HDF5 file
    file = H5Fcreate(s1.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    //Write xdmf file
    output_v << "<?xml version=\"1.0\"?>" << std::endl
             << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>" << std::endl
             << "<Xdmf Version=\"2.0\" xmlns:xi=\"http://www.w3.org/2001/XInclude\" >" << std::endl
             << "<Domain>" << std::endl
             << "  <Grid>"  << std::endl
             << "    <Topology TopologyType=\"" << elType << "\" NumberOfElements=\"" << numElemFine << "\" >" << std::endl
             << "      <DataItem Format=    \"HDF\" NumberType=\"int\" Dimensions=\"" << numElemFine << " " << nElNodes << "\" >" << std::endl;
    
    //Connectivity
    output_v << "        " << s1 << ":/connec" << std::endl;
    //Start connectivity
    dataspace = H5Screate_simple(2, connec2Dims, NULL);
    dataset = H5Dcreate2(file, "/connec", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &connec2[0]);
    status = H5Dclose(dataset);
    status = H5Sclose(dataspace);
    //End connectivity

    output_v << "      </DataItem>" << std::endl
             << "    </Topology>" << std::endl
             << "    <Geometry GeometryType=\"XYZ\">" << std::endl
             << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< numNodesFine <<" 3\">" << std::endl;

    //Coordinates   
    output_v << "        " << s1 << ":/coords" << std::endl;
    //Start coordinates
    dataspace = H5Screate_simple(2, pointVectorDims, NULL);
    dataset = H5Dcreate2(file, "/coords", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointVector[0]);
    status = H5Dclose(dataset);
    status = H5Sclose(dataspace);
    //End coordinates
    output_v << "      </DataItem>" << std::endl
             << "    </Geometry>" << std::endl;

    //LISTS
    //Velocity
    if (fineModel.printVelocity){
        output_v << "    <Attribute Name=\"Velocity\" Center=\"Node\" AttributeType=\"Vector\" >" << std::endl
                 << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< numNodesFine <<" 3\">" << std::endl;  
        output_v << "        " << s1 << ":/velocity" << std::endl;     
        //Start Velocity
        for (int i = 0; i < numNodesFine; i++) {
            pointVector[3*i  ] = (*nodesFine_)[i] -> getVelocity(0);
            pointVector[3*i+1] = (*nodesFine_)[i] -> getVelocity(1);
            pointVector[3*i+2] = 0.0;
        };
        dataspace = H5Screate_simple(2, pointVectorDims, NULL);
        dataset = H5Dcreate2(file, "/velocity", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointVector[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
        //End Velocity
        output_v << "      </DataItem>" << std::endl
                 << "    </Attribute>" << std::endl;
    }

    //Acceleration
    if (fineModel.printAcceleration){
        output_v << "    <Attribute Name=\"Acceleration\" Center=\"Node\" AttributeType=\"Vector\" >" << std::endl
                 << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< numNodesFine <<" 3\">" << std::endl;  
        output_v << "        " << s1 << ":/acceleration" << std::endl;     
        //Start Acceleration
        for (int i = 0; i < numNodesFine; i++) {
            pointVector[3*i  ] = (*nodesFine_)[i] -> getAcceleration(0);
            pointVector[3*i+1] = (*nodesFine_)[i] -> getAcceleration(1);
            pointVector[3*i+2] = 0.0;
        };
        dataspace = H5Screate_simple(2, pointVectorDims, NULL);
        dataset = H5Dcreate2(file, "/acceleration", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointVector[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
        //End Acceleration
        output_v << "      </DataItem>" << std::endl
                 << "    </Attribute>" << std::endl;
    }

    //Real Velocity
    if (fineModel.printRealVelocity){
        output_v << "    <Attribute Name=\"Real Velocity\" Center=\"Node\" AttributeType=\"Vector\" >" << std::endl
                 << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< numNodesFine <<" 3\">" << std::endl;  
        output_v << "        " << s1 << ":/realVelocity" << std::endl;     
        //Start Real Velocity
        for (int i = 0; i < numNodesFine; i++) {
            pointVector[3*i  ] = (*nodesFine_)[i] -> getVelocityArlequin(0);
            pointVector[3*i+1] = (*nodesFine_)[i] -> getVelocityArlequin(1);
            pointVector[3*i+2] = 0.0;
        };
        dataspace = H5Screate_simple(2, pointVectorDims, NULL);
        dataset = H5Dcreate2(file, "/realVelocity", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointVector[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
        //End Real Velocity
        output_v << "      </DataItem>" << std::endl
                 << "    </Attribute>" << std::endl;
    }

    //Lagrange Multipliers
    if (fineModel.printLagrangeMultipliers){
        output_v << "    <Attribute Name=\"Lagrange Multipliers\" Center=\"Node\" AttributeType=\"Vector\" >" << std::endl
                 << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< numNodesFine <<" 3\">" << std::endl;  
        output_v << "        " << s1 << ":/lagrangeMultiplers" << std::endl;     
        //Start Lagrange Multipliers
        for (int i = 0; i < numNodesFine; i++) {
            pointVector[3*i  ] = (*nodesFine_)[i] -> getLagrangeMultiplier(0);
            pointVector[3*i+1] = (*nodesFine_)[i] -> getLagrangeMultiplier(1);
            pointVector[3*i+2] = 0.0;
        };
        dataspace = H5Screate_simple(2, pointVectorDims, NULL);
        dataset = H5Dcreate2(file, "/lagrangeMultiplers", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointVector[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
        //End Lagrange Multipliers
        output_v << "      </DataItem>" << std::endl
                 << "    </Attribute>" << std::endl;
    }

    //Normal vector
    if (fineModel.printInnerNormal){
        output_v << "    <Attribute Name=\"Normal Vector\" Center=\"Node\" AttributeType=\"Vector\" >" << std::endl
                 << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< numNodesFine <<" 3\">" << std::endl;  
        output_v << "        " << s1 << ":/normalVector" << std::endl;     
        //Start Normal vector
        for (int i = 0; i < numNodesFine; i++) {
            VecDouble n = (*nodesFine_)[i] -> getInnerNormal();
            pointVector[3*i  ] = n[0];
            pointVector[3*i+1] = n[1];
            pointVector[3*i+2] = 0.0;
        };
        dataspace = H5Screate_simple(2, pointVectorDims, NULL);
        dataset = H5Dcreate2(file, "/normalVector", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointVector[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
        //End Normal Vector
        output_v << "      </DataItem>" << std::endl
                 << "    </Attribute>" << std::endl;
    }

    //Mesh Velocity
    if (fineModel.printMeshVelocity){
        output_v << "    <Attribute Name=\"Mesh Velocity\" Center=\"Node\" AttributeType=\"Vector\" >" << std::endl
                 << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< numNodesFine <<" 3\">" << std::endl;  
        output_v << "        " << s1 << ":/meshVelocity" << std::endl;     
        //Start mesh Velocity
        for (int i = 0; i < numNodesFine; i++) {
            pointVector[3*i  ] = (*nodesFine_)[i] -> getMeshVelocity(0);
            pointVector[3*i+1] = (*nodesFine_)[i] -> getMeshVelocity(1);
            pointVector[3*i+2] = 0.0;
        };
        dataspace = H5Screate_simple(2, pointVectorDims, NULL);
        dataset = H5Dcreate2(file, "/meshVelocity", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointVector[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
        //End Mesh Velocity
        output_v << "      </DataItem>" << std::endl
                 << "    </Attribute>" << std::endl;
    };

    //Mesh Displacement
    if (fineModel.printMeshDisplacement){
        output_v << "    <Attribute Name=\"Mesh Displacement\" Center=\"Node\" AttributeType=\"Vector\" >" << std::endl
                 << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< numNodesFine <<" 3\">" << std::endl;  
        output_v << "        " << s1 << ":/meshDisplacement" << std::endl;     
        //Start Mesh displacement
        for (int i = 0; i < numNodesFine; i++) {
            VecDouble x = (*nodesFine_)[i] -> getCoordinates();
            VecDouble xi = (*nodesFine_)[i] -> getInitialCoordinates();
            pointVector[3*i  ] = x[0] - xi[0];
            pointVector[3*i+1] = x[1] - xi[1];
            pointVector[3*i+2] = 0.0;
        };
        dataspace = H5Screate_simple(2, pointVectorDims, NULL);
        dataset = H5Dcreate2(file, "/meshDisplacement", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointVector[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
        //End Mesh displacement
        output_v << "      </DataItem>" << std::endl
                 << "    </Attribute>" << std::endl;
    };

    //Element correspondence
    if (fineModel.printElementCorrespondence){
        output_v << "    <Attribute Name=\"Element\" Center=\"Node\" AttributeType=\"Scalar\" >" << std::endl
                 << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\"1 "<< numNodesFine <<"\">" << std::endl;  
        output_v << "        " << s1 << ":/elementCorresp" << std::endl;     
        //Start Element correspondence
        for (int i = 0; i < numNodesFine; i++) pointScalar[i] = (*nodesFine_)[i] -> getNodalElemCorrespondence();
        dataspace = H5Screate_simple(2, pointScalarDims, NULL);
        dataset = H5Dcreate2(file, "/elementCorresp", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointScalar[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
        //End Element correspondence
        output_v << "      </DataItem>" << std::endl
                 << "    </Attribute>" << std::endl;
    }

    //Dist Function
    if (fineModel.printDistFunction){
        output_v << "    <Attribute Name=\"Dist Function\" Center=\"Node\" AttributeType=\"Scalar\" >" << std::endl
                 << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\"1 "<< numNodesFine <<"\">" << std::endl;  
        output_v << "        " << s1 << ":/distfunction" << std::endl;     
        //Start Dist Function
        for (int i = 0; i < numNodesFine; i++) pointScalar[i] = (*nodesFine_)[i] -> getDistFunction();
        dataspace = H5Screate_simple(2, pointScalarDims, NULL);
        dataset = H5Dcreate2(file, "/distfunction", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointScalar[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
        //End Dist Function
        output_v << "      </DataItem>" << std::endl
                 << "    </Attribute>" << std::endl;
    }

    //Weight Function
    if (fineModel.printEnergyWeightFunction){
        output_v << "    <Attribute Name=\"Weight Function\" Center=\"Node\" AttributeType=\"Scalar\" >" << std::endl
                 << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\"1 "<< numNodesFine <<"\">" << std::endl;  
        output_v << "        " << s1 << ":/weightFunction" << std::endl;     
        //Start Weigth Function
        for (int i = 0; i < numNodesFine; i++) pointScalar[i] = (*nodesFine_)[i] -> getWeightFunction();
        dataspace = H5Screate_simple(2, pointScalarDims, NULL);
        dataset = H5Dcreate2(file, "/weightFunction", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointScalar[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
        //End Weight Function
        output_v << "      </DataItem>" << std::endl
                 << "    </Attribute>" << std::endl;
    }

    //Pressure
    if (fineModel.printPressure){
        output_v << "    <Attribute Name=\"Pressure\" Center=\"Node\" AttributeType=\"Scalar\" >" << std::endl
                 << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< numNodesFine <<" 1\">" << std::endl;  
        output_v << "        " << s1 << ":/pressure" << std::endl;     
        //Start Pressure
        for (int i = 0; i < numNodesFine; i++) pointScalar[i] = (*nodesFine_)[i] -> getPressure();
        dataspace = H5Screate_simple(2, pointScalarDims, NULL);
        dataset = H5Dcreate2(file, "/pressure", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointScalar[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
        //End Pressure
        output_v << "      </DataItem>" << std::endl
                 << "    </Attribute>" << std::endl;
    }

    //Real Pressure
    if (fineModel.printRealPressure){
        output_v << "    <Attribute Name=\"Real Pressure\" Center=\"Node\" AttributeType=\"Scalar\" >" << std::endl
                 << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< numNodesFine <<" 1\">" << std::endl;  
        output_v << "        " << s1 << ":/realPressure" << std::endl;     
        //Start Real Pressure
        for (int i = 0; i < numNodesFine; i++) pointScalar[i] = (*nodesFine_)[i] -> getPressureArlequin();
        dataspace = H5Screate_simple(2, pointScalarDims, NULL);
        dataset = H5Dcreate2(file, "/realPressure", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointScalar[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
        //End Real Pressure
        output_v << "      </DataItem>" << std::endl
                 << "    </Attribute>" << std::endl;
    }

    if (fineModel.printVorticity){
        output_v << "    <Attribute Name=\"Vorticity\" Center=\"Node\" AttributeType=\"Scalar\" >" << std::endl
                 << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< numNodesFine <<" 1\">" << std::endl;  
        output_v << "        " << s1 << ":/vorticity" << std::endl;     
        //Start Real Pressure
        for (int i = 0; i < numNodesFine; i++) pointScalar[i] = (*nodesFine_)[i] -> getVorticity();
        dataspace = H5Screate_simple(2, pointScalarDims, NULL);
        dataset = H5Dcreate2(file, "/vorticity", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointScalar[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
        //End Real Pressure
        output_v << "      </DataItem>" << std::endl
                 << "    </Attribute>" << std::endl;
    };

    //Process
    if (fineModel.printProcess){
        output_v << "    <Attribute Name=\"Process\" Center=\"Cell\" AttributeType=\"Scalar\" >" << std::endl
                 << "      <DataItem Format=\"HDF\" NumberType=\"int\" Dimensions=\""<< numElemFine <<" 1\">" << std::endl;  
        output_v << "        " << s1 << ":/process" << std::endl;
        //Start Process
        dataspace = H5Screate_simple(2, pointCellScalarDims, NULL);
        dataset = H5Dcreate2(file, "/process", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &domDecompFine.first[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
        //End Process
        output_v << "      </DataItem>" << std::endl
                 << "    </Attribute>" << std::endl;
    }

    //WeightFunction2
    if (fineModel.printEnergyWeightFunction){
        output_v << "    <Attribute Name=\"Weight Function\" Center=\"Cell\" AttributeType=\"Scalar\" >" << std::endl
                 << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< numElemFine <<" 1\">" << std::endl;  
        output_v << "        " << s1 << ":/weightFunction2" << std::endl;     
        //Start Weight Function
        for (int i = 0; i < numElemFine; i++) douCellScalar[i] = elementsFine_[i] -> getIntegPointWeightFunction(0);
        dataspace = H5Screate_simple(2, pointCellScalarDims, NULL);
        dataset = H5Dcreate2(file, "/weightFunction2", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &douCellScalar[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
        //End Weight Function
        output_v << "      </DataItem>" << std::endl
                 << "    </Attribute>" << std::endl;
    }

    //Lines
    if (fineModel.printLines){
        output_v << "    <Attribute Name=\"Lines\" Center=\"Cell\" AttributeType=\"Scalar\" >" << std::endl
                 << "      <DataItem Format=\"HDF\" NumberType=\"int\" Dimensions=\""<< numElemFine <<" 1\">" << std::endl;  
        output_v << "        " << s1 << ":/lines" << std::endl;     
        //Start Lines
        for (int i = 0; i < numElemFine; i++){
            int res = 0;
            for (int j = 0; j < numBoundElemFine; j++)
               if (boundaryFine_[j] -> getElement() == i) res = boundaryFine_[j] -> getBoundaryGroup();
            intCellScalar[i] = res;
        };
        dataspace = H5Screate_simple(2, pointCellScalarDims, NULL);
        dataset = H5Dcreate2(file, "/lines", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &intCellScalar[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
        //End Lines
        output_v << "      </DataItem>" << std::endl
                 << "    </Attribute>" << std::endl;
    }

    // Jacobian
    if (fineModel.printJacobian){
        output_v << "    <Attribute Name=\"Jacobian\" Center=\"Cell\" AttributeType=\"Scalar\" >" << std::endl
                 << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< numElemFine <<" 1\">" << std::endl;  
        output_v << "        " << s1 << ":/jacobian" << std::endl;     
        //Start Weight Function
        for (int i = 0; i < numElemFine; i++) douCellScalar[i] = elementsFine_[i] -> getJacobian();
        dataspace = H5Screate_simple(2, pointCellScalarDims, NULL);
        dataset = H5Dcreate2(file, "/jacobian", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &douCellScalar[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
        //End Weight Function
        output_v << "      </DataItem>" << std::endl
                 << "    </Attribute>" << std::endl;
    }

    //Glue Zone
    if (fineModel.printGlueZone){
        output_v << "    <Attribute Name=\"Glue Zone\" Center=\"Cell\" AttributeType=\"Scalar\" >" << std::endl
                 << "      <DataItem Format=\"HDF\" NumberType=\"int\" Dimensions=\""<< numElemFine <<" 1\">" << std::endl;  
        output_v << "        " << s1 << ":/glueZone" << std::endl;     
        //Start Lines
        for (int i = 0; i < numElemFine; i++) intCellScalar[i] = elementsFine_[i] -> getGlueZoneInt();
        dataspace = H5Screate_simple(2, pointCellScalarDims, NULL);
        dataset = H5Dcreate2(file, "/glueZone", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &intCellScalar[0]);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);
        //End Lines
        output_v << "      </DataItem>" << std::endl
                 << "    </Attribute>" << std::endl;
    }

    // END FILE
    output_v << "  </Grid>" << std::endl
             << "</Domain>" << std::endl
             << "</Xdmf>" << std::endl;

    //Delete auxiliary vectors
    delete [] pointVector;
    delete [] connec2;
    delete [] pointScalar;
    delete [] intCellScalar;
    delete [] douCellScalar;

    //End HDF5 file
    status = H5Fclose(file);

    return;
};


//------------------------------------------------------------------------------
//------------SETS COARSE/FINE MESHES AND GETS ITS BASIC INFORMATIONS-----------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Arlequin<DIM,DEG>::setFluidModels(FluidMesh& coarse, FluidMesh& fine){

    coarseModel = coarse;
    fineModel = fine;

    //Gets Fine and Coarse models basic information
    numElemCoarse = coarseModel.elements_.size();
    numElemFine   = fineModel.elements_.size();
    numNodesCoarse = coarseModel.nodes_.size();
    numNodesFine   = fineModel.nodes_.size();
 
    nodesCoarse_  = &coarseModel.nodes_;
    nodesFine_    = &fineModel.nodes_;

    elementsCoarse_ = coarseModel.elements_;
    elementsFine_   = fineModel.elements_;
 
    boundaryCoarse_ = coarseModel.boundary_;
    boundaryFine_   = fineModel.boundary_;
    fProbType = fineModel.getProblemType();

    numBoundElemFine = boundaryFine_.size();
    numBoundElemCoarse = boundaryCoarse_.size();

    domDecompCoarse = coarseModel.getDomainDecomposition();
    domDecompFine = fineModel.getDomainDecomposition();
    
    numTimeSteps = fineModel.getNumberOfTimeSteps();
    dTime = fineModel.getTimeStep();

    for (int i=0; i < numElemFine; i++){
        elementsFine_[i] -> setModel(true);
    };
    for (int i=0; i < numElemCoarse; i++){
        elementsCoarse_[i] -> setModel(false);
    };

    //Sets the element boxes for all elements in both coarse and fine models
    setElementBoxes();
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    // // Rotating fine mesh
    // for (int i = 0; i < numNodesFine; ++i){
    //     typename Nodes::VecLocD x, xn;
    //     x = (*nodesFine_)[i] -> getCoordinates();       
    
    //     double a = -90 * pi / 180;

    //     xn(0) = 0.0 + (x(0)-0.0) * cos(a) - (x(1)-0.0) * sin(a);
    //     xn(1) = 0.0 + (x(0)-0.0) * sin(a) + (x(1)-0.0) * cos(a);

    //     (*nodesFine_)[i] -> setCoordinates(xn);
    //     (*nodesFine_)[i] -> setPreviousCoordinates(0,xn(0));
    //     (*nodesFine_)[i] -> setPreviousCoordinates(1,xn(1));
    // }

    parametersFine = &fineModel.fluidParameters;
    parametersCoarse = &coarseModel.fluidParameters;

    // std::cout << "AA1 " << rank << std::endl;
    // MPI_Barrier(PETSC_COMM_WORLD);

    setSignaledDistance();

    // std::cout << "AA2 " << rank << std::endl;
    // MPI_Barrier(PETSC_COMM_WORLD);

    //Construct the glue zone based on some defined criterion
    setCouplingZone();
    // std::cout << "AA3 " << rank << std::endl;
    // MPI_Barrier(PETSC_COMM_WORLD);
    //Computes the Weight function for all the finite elements
    setWeightFunction(16.); 
    // std::cout << "AA4 " << rank << std::endl;
    // MPI_Barrier(PETSC_COMM_WORLD);
    //Update domain decomposition - Start
    int size;
    MPI_Comm_size(PETSC_COMM_WORLD, &size);

    Vec  b;
    PetscErrorCode    ierr;
    PetscInt start[size], end[size];
    PetscInt numDOF = (DIM+1)*numNodesCoarse + (DIM+1)*numNodesFine + DIM*numNodesGlueZoneFine;

    ierr = VecCreate(PETSC_COMM_WORLD,&b);
    ierr = VecSetSizes(b,PETSC_DECIDE,numDOF);
    VecSetFromOptions(b);
    ierr = VecGetOwnershipRange(b,&start[rank],&end[rank]);

    VecDestroy(&b);
    
    for (int i = 0; i < size; ++i){
        MPI_Bcast(&start[i],1,MPI_INT,i,PETSC_COMM_WORLD);
        MPI_Bcast(&end[i],1,MPI_INT,i,PETSC_COMM_WORLD);
    }
    // std::cout << "AA5 " << rank << std::endl;
    // MPI_Barrier(PETSC_COMM_WORLD);
    for (int i = 0; i < numElemCoarse; i++){
        VecInt connec = elementsCoarse_[i] -> getConnectivity();
        for (int j = 0; j < size; j++){
            if (((DIM+1)*connec[0] >= start[j]) && ((DIM+1)*connec[0] <= end[j])) {
                domDecompCoarse.first[i] = j;
                break;
            }
        }
    }    
    for (int i = 0; i < numElemFine; i++){
        VecInt connec = elementsFine_[i] -> getConnectivity();
        for (int j = 0; j < size; j++){
            if (((DIM+1)*(numNodesCoarse + connec[0]) >= start[j]) && ((DIM+1)*(numNodesCoarse + connec[0]) <= end[j])) {
                domDecompFine.first[i] = j;
                break;
            }
        }
    }    
    // std::cout << "AA6 " << rank << std::endl;
    // MPI_Barrier(PETSC_COMM_WORLD);
    //Update domain decomposition - End

    if(rank == 0){
        std::cout << "---------------------ARLEQUIN DATA---------------------" 
                  << std::endl;
        std::cout << "Coarse Model: " << numNodesCoarse << " nodes, " << numElemCoarse << " elements." << std::endl;
        std::cout << "Fine Model: " << numNodesFine << " nodes, " << numElemFine << " elements." << std::endl;
        std::cout << "Lagrange Multipliers: " << numNodesGlueZoneFine << " nodes, " << numElemGlueZoneFine << " elements." << std::endl;
        std::cout << "Number of Degrees of Freedom: " << numDOF << std::endl; 
    }

    if(rank == 0){
        printResultsCoarse(10);
        printResultsFine(10);
    }


    // (*nodesFine_)[269] -> setWeightFunction(0.);
    // (*nodesFine_)[80] -> setWeightFunction(0.);





};

//------------------------------------------------------------------------------
//----------------------COMPUTES DRAG AND LIFT COEFFICIENTS---------------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Arlequin<DIM,DEG>::dragAndLiftCoefficients(std::ofstream& dragLift){

    double dragCoefficient = 0.;
    double liftCoefficient = 0.;
    double pressureDragCoefficient = 0.;
    double pressureLiftCoefficient = 0.;
    double frictionDragCoefficient = 0.;
    double frictionLiftCoefficient = 0.;
    double theta = 0.;
    double pitchingMomentCoefficient = 0.;
    double pMom = 0.;
    double per = 0.;

    double rhoInf = 1.0;
    double velocityInf[DIM];
    velocityInf[0] = 1;
    velocityInf[1] = 0.;

    for (int jel = 0; jel < numBoundElemFine; jel++){   
        

        
        double dForce = 0.;
        double lForce = 0.;
        double pDForce = 0.;
        double pLForce = 0.;
        double fDForce = 0.;
        double fLForce = 0.;
        double aux_Mom = 0.;
        double aux_Per = 0.;



        for (int i=0; i<fineModel.numberOfLines; i++){
            if (boundaryFine_[jel] -> getBoundaryGroup() == fineModel.dragAndLiftBoundary[i]){
                
                int iel = boundaryFine_[jel] -> getElement();
                elementsFine_[iel] -> computeDragAndLiftForces(pDForce, pLForce, fDForce, fLForce, dForce, lForce, aux_Mom, aux_Per);

                pMom += aux_Mom;
                per += aux_Per;
            };
        };
        
        pressureDragCoefficient += pDForce ;/// 
            //(0.5 * rhoInf * velocityInf[0] * velocityInf[0]);
        pressureLiftCoefficient += pLForce ;/// 
            //(0.5 * rhoInf * velocityInf[0] * velocityInf[0]);
        
        frictionDragCoefficient += fDForce / 
            (0.5 * rhoInf * velocityInf[0] * velocityInf[0] * 1);
        frictionLiftCoefficient += fLForce / 
            (0.5 * rhoInf * velocityInf[0] * velocityInf[0] * 1);
        
        dragCoefficient += dForce / 
            (0.5 * rhoInf * velocityInf[0] * velocityInf[0] * 1);
        liftCoefficient += lForce / 
            (0.5 * rhoInf * velocityInf[0] * velocityInf[0] * 1);
        
    };
    
    pitchingMomentCoefficient = pMom / (rhoInf * velocityInf[0] * velocityInf[0] * per);

    if (rank == 0) {
        const int timeWidth = 13;
        const int numWidth = 13;
        dragLift << std::setprecision(5) << std::scientific;
        dragLift << std::left << std::setw(timeWidth) << iTimeStep * dTime;
        dragLift << std::setw(numWidth) << pressureDragCoefficient;
        dragLift << std::setw(numWidth) << pressureLiftCoefficient;
        dragLift << std::setw(numWidth) << frictionDragCoefficient;
        dragLift << std::setw(numWidth) << frictionLiftCoefficient;
        dragLift << std::setw(numWidth) << dragCoefficient;
        dragLift << std::setw(numWidth) << liftCoefficient;
        dragLift << std::setw(numWidth) << pitchingMomentCoefficient;
        //dragLift << std::setw(numWidth) << theta;
        dragLift << std::endl;
    }

    return;
}

//------------------------------------------------------------------------------
//----------------COMPUTE ARLEQUIN COUPLED NAVIER-STOKES PROBLEM----------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Arlequin<DIM,DEG>::setMatVecValuesCoarseModel(MatrixDouble &matrix, VecDouble &rhs, VecInt &connec){

    //Disperse local contributions into the global matrix
    for (int i = 0; i < nElNodes; i++){
        for (int j = 0; j < nElNodes; j++){
            for (int k = 0; k < DIM; k++){
                for (int l = 0; l < DIM; l++){
                    //Matrix K and C            
                    PetscInt dof_i = (DIM+1) * connec[i] + k;
                    PetscInt dof_j = (DIM+1) * connec[j] + l;
                    MatSetValues(A,1,&dof_i,1,&dof_j,&matrix(DIM*i+k,DIM*j+l),ADD_VALUES);
                }
                //Matrix G and Gt
                PetscInt dof_i = (DIM+1) * connec[i] + k;
                PetscInt dof_j = (DIM+1) * connec[j] + DIM;
                MatSetValues(A,1,&dof_i,1,&dof_j,&matrix(DIM*i+k,DIM*nElNodes+j),ADD_VALUES);
                MatSetValues(A,1,&dof_j,1,&dof_i,&matrix(DIM*nElNodes+j,DIM*i+k),ADD_VALUES);
            }
            // Matrix Q
            PetscInt dof_i = (DIM+1) * connec[i] + DIM;
            PetscInt dof_j = (DIM+1) * connec[j] + DIM;
            MatSetValues(A,1,&dof_j,1,&dof_i,&matrix(DIM*nElNodes+i,DIM*nElNodes+j),ADD_VALUES);
        };
        for (int k = 0; k < DIM; k++){          
            //Rhs vector
            PetscInt dof_i = (DIM+1) * connec[i] + k;
            VecSetValues(b,1,&dof_i,&rhs[DIM*i+k],ADD_VALUES);
        }
        PetscInt dof_i = (DIM+1) * connec[i] + DIM;
        VecSetValues(b,1,&dof_i,&rhs[DIM*nElNodes+i],ADD_VALUES);
    };

    return;
}


template<int DIM, int DEG>
void Arlequin<DIM,DEG>::setVecValuesCoarseModel(VecDouble &rhs, VecInt &connec){

    //Disperse local contributions into the global matrix
    for (int i = 0; i < nElNodes; i++){
        for (int k = 0; k < DIM; k++){          
            //Rhs vector
            PetscInt dof_i = (DIM+1) * connec[i] + k;
            VecSetValues(b,1,&dof_i,&rhs[DIM*i+k],ADD_VALUES);
        }
        PetscInt dof_i = (DIM+1) * connec[i] + DIM;
        VecSetValues(b,1,&dof_i,&rhs[DIM*nElNodes+i],ADD_VALUES);
    };

    return;
}

template<int DIM, int DEG>
void Arlequin<DIM,DEG>::setMatValuesCoarseModel(MatrixDouble &matrix, VecInt &connec){

    //Disperse local contributions into the global matrix
    for (int i = 0; i < nElNodes; i++){
        for (int j = 0; j < nElNodes; j++){
            for (int k = 0; k < DIM; k++){
                for (int l = 0; l < DIM; l++){
                    //Matrix K and C            
                    PetscInt dof_i = (DIM+1) * connec[i] + k;
                    PetscInt dof_j = (DIM+1) * connec[j] + l;
                    MatSetValues(A,1,&dof_i,1,&dof_j,&matrix(DIM*i+k,DIM*j+l),ADD_VALUES);
                }
                //Matrix G and Gt
                PetscInt dof_i = (DIM+1) * connec[i] + k;
                PetscInt dof_j = (DIM+1) * connec[j] + DIM;
                MatSetValues(A,1,&dof_i,1,&dof_j,&matrix(DIM*i+k,DIM*nElNodes+j),ADD_VALUES);
                MatSetValues(A,1,&dof_j,1,&dof_i,&matrix(DIM*nElNodes+j,DIM*i+k),ADD_VALUES);
            }
            // Matrix Q
            PetscInt dof_i = (DIM+1) * connec[i] + DIM;
            PetscInt dof_j = (DIM+1) * connec[j] + DIM;
            MatSetValues(A,1,&dof_j,1,&dof_i,&matrix(DIM*nElNodes+i,DIM*nElNodes+j),ADD_VALUES);
        };
    };

    return;
}




//------------------------------------------------------------------------------
//----------------COMPUTE ARLEQUIN COUPLED NAVIER-STOKES PROBLEM----------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Arlequin<DIM,DEG>::setMatVecValuesCoarseModelPoisson(MatrixDouble &matrix, VecDouble &rhs, VecInt &connec){

    //Disperse local contributions into the global matrix
    for (int i = 0; i < nElNodes; i++){
        PetscInt dof_i = connec[i];
        for (int j = 0; j < nElNodes; j++){
            //Matrix K
            PetscInt dof_j = connec[j];
            MatSetValues(A,1,&dof_i,1,&dof_j,&matrix(i,j),ADD_VALUES);
        };
        //Rhs vector
        VecSetValues(b,1,&dof_i,&rhs[i],ADD_VALUES);
    };

    return;
}

//------------------------------------------------------------------------------
//----------------COMPUTE ARLEQUIN COUPLED NAVIER-STOKES PROBLEM----------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Arlequin<DIM,DEG>::setMatVecValuesFineModel(MatrixDouble &matrix, VecDouble &rhs, VecInt &connec){

    //Disperse local contributions into the global matrix
    for (int i = 0; i < nElNodes; i++){
        for (int j = 0; j < nElNodes; j++){
            for (int k = 0; k < DIM; k++){
                for (int l = 0; l < DIM; l++){
                    //Matrix K and C
                    PetscInt dof_i = (DIM+1) * numNodesCoarse + (DIM+1) * connec[i] + k;
                    PetscInt dof_j = (DIM+1) * numNodesCoarse + (DIM+1) * connec[j] + l;
                    MatSetValues(A,1,&dof_i,1,&dof_j,&matrix(DIM*i+k,DIM*j+l),ADD_VALUES);
                }          
                //Matrix G and Gt
                PetscInt dof_i = (DIM+1) * numNodesCoarse + (DIM+1) * connec[i] + k;
                PetscInt dof_j = (DIM+1) * numNodesCoarse + (DIM+1) * connec[j] + DIM;
                MatSetValues(A,1,&dof_i,1,&dof_j,&matrix(DIM*i+k,DIM*nElNodes+j),ADD_VALUES);
                MatSetValues(A,1,&dof_j,1,&dof_i,&matrix(DIM*nElNodes+j,DIM*i+k),ADD_VALUES);
            }
            //Matrix Q
            PetscInt dof_i = (DIM+1) * numNodesCoarse + (DIM+1) * connec[i] + DIM;
            PetscInt dof_j = (DIM+1) * numNodesCoarse + (DIM+1) * connec[j] + DIM;
            MatSetValues(A,1,&dof_j,1,&dof_i,&matrix(DIM*nElNodes+i,DIM*nElNodes+j),ADD_VALUES);
        };
        for (int k = 0; k < DIM; k++){  
            ///Rhs vector
            PetscInt dof_i = (DIM+1) * numNodesCoarse + (DIM+1) * connec[i] + k;
            VecSetValues(b,1,&dof_i,&rhs[DIM*i+k],ADD_VALUES);
        }
        PetscInt dof_i = (DIM+1) * numNodesCoarse + (DIM+1) * connec[i] + DIM;
        VecSetValues(b,1,&dof_i,&rhs[DIM*nElNodes+i],ADD_VALUES);
    }; 

    return;
};

template<int DIM, int DEG>
void Arlequin<DIM,DEG>::setVecValuesFineModel(VecDouble &rhs, VecInt &connec){

    //Disperse local contributions into the global matrix
    for (int i = 0; i < nElNodes; i++){
        for (int k = 0; k < DIM; k++){  
            ///Rhs vector
            PetscInt dof_i = (DIM+1) * numNodesCoarse + (DIM+1) * connec[i] + k;
            VecSetValues(b,1,&dof_i,&rhs[DIM*i+k],ADD_VALUES);
        }
        PetscInt dof_i = (DIM+1) * numNodesCoarse + (DIM+1) * connec[i] + DIM;
        VecSetValues(b,1,&dof_i,&rhs[DIM*nElNodes+i],ADD_VALUES);
    }; 

    return;
};

template<int DIM, int DEG>
void Arlequin<DIM,DEG>::setMatValuesFineModel(MatrixDouble &matrix, VecInt &connec){

    //Disperse local contributions into the global matrix
    for (int i = 0; i < nElNodes; i++){
        for (int j = 0; j < nElNodes; j++){
            for (int k = 0; k < DIM; k++){
                for (int l = 0; l < DIM; l++){
                    //Matrix K and C
                    PetscInt dof_i = (DIM+1) * numNodesCoarse + (DIM+1) * connec[i] + k;
                    PetscInt dof_j = (DIM+1) * numNodesCoarse + (DIM+1) * connec[j] + l;
                    MatSetValues(A,1,&dof_i,1,&dof_j,&matrix(DIM*i+k,DIM*j+l),ADD_VALUES);
                }          
                //Matrix G and Gt
                PetscInt dof_i = (DIM+1) * numNodesCoarse + (DIM+1) * connec[i] + k;
                PetscInt dof_j = (DIM+1) * numNodesCoarse + (DIM+1) * connec[j] + DIM;
                MatSetValues(A,1,&dof_i,1,&dof_j,&matrix(DIM*i+k,DIM*nElNodes+j),ADD_VALUES);
                MatSetValues(A,1,&dof_j,1,&dof_i,&matrix(DIM*nElNodes+j,DIM*i+k),ADD_VALUES);
            }
            //Matrix Q
            PetscInt dof_i = (DIM+1) * numNodesCoarse + (DIM+1) * connec[i] + DIM;
            PetscInt dof_j = (DIM+1) * numNodesCoarse + (DIM+1) * connec[j] + DIM;
            MatSetValues(A,1,&dof_j,1,&dof_i,&matrix(DIM*nElNodes+i,DIM*nElNodes+j),ADD_VALUES);
        };
    }; 

    return;
};

//------------------------------------------------------------------------------
//----------------COMPUTE ARLEQUIN COUPLED NAVIER-STOKES PROBLEM----------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Arlequin<DIM,DEG>::setMatVecValuesFineModelPoisson(MatrixDouble &matrix, VecDouble &rhs, VecInt &connec){

    //Disperse local contributions into the global matrix
    for (int i = 0; i < nElNodes; i++){
        PetscInt dof_i = numNodesCoarse + connec[i];
        for (int j = 0; j < nElNodes; j++){
            //Matrix K
            PetscInt dof_j = numNodesCoarse + connec[j];
            MatSetValues(A,1,&dof_i,1,&dof_j,&matrix(i,j),ADD_VALUES);
        };
        ///Rhs vector
        VecSetValues(b,1,&dof_i,&rhs[i],ADD_VALUES);
    }; 

    return;
};

//------------------------------------------------------------------------------
//----------------COMPUTE ARLEQUIN COUPLED NAVIER-STOKES PROBLEM----------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Arlequin<DIM,DEG>::setMatVecValuesLagMultFineFine(MatrixDouble &Ajac2, MatrixDouble &localMV_mat, 
                                                       MatrixDouble &ArlequinA1, MatrixDouble &ArlequinA2, 
                                                       VecDouble &Rhs2, VecDouble &rhsLagMult2,
                                                       VecDouble &localMV_vec, VecDouble &RhsArlequin2,
                                                       VecInt &connec, VecInt &connecL){

    double &alpha_f = parametersFine -> getAlphaF();
    double &alpha_m = parametersFine -> getAlphaM();
    double &gamma = parametersFine -> getGamma();
    
    double integ = alpha_f * gamma * dTime;
    //Disperse local contributions into the global matrix
    for (int i = 0; i < nElNodes; i++){
        for (int j = 0; j < nElNodes; j++){
            for (int k = 0; k < DIM; k++){
                for (int l = 0; l < DIM; l++){
                    //COUPLING OPERATOR
                    if (fabs(Ajac2(DIM*i+k,DIM*j+l)) >= 1.e-15){
                        PetscInt d_i = (DIM+1) * numNodesCoarse + (DIM+1) * numNodesFine + DIM * connecL[i] + k;
                        PetscInt d_j = (DIM+1) * numNodesCoarse + (DIM+1) * connec[j] + l;
                        double value = Ajac2(DIM*i+k,DIM*j+l)*integ;
                        MatSetValues(A,1,&d_i,1,&d_j,&value,ADD_VALUES);
                        MatSetValues(A,1,&d_j,1,&d_i,&Ajac2(DIM*i+k,DIM*j+l),ADD_VALUES);
                    };
                }
                //SUPG STABILIZATION
                if (fabs(localMV_mat(DIM*i+k,DIM*j+k)) >= 1.e-15){
                    PetscInt d_i = (DIM+1) * numNodesCoarse + (DIM+1) * numNodesFine + DIM*connecL[i] + k;
                    PetscInt d_j = (DIM+1) * numNodesCoarse + (DIM+1) * connec[j] + k;
                    MatSetValues(A,1,&d_j,1,&d_i,&localMV_mat(2*i  ,2*j  ),ADD_VALUES);
                };
                //PSPG STABILIZATION
                if (fabs(localMV_mat(DIM*i+k,DIM*nElNodes+j)) >= 1.e-15){
                    PetscInt dof_i = (DIM+1) * numNodesCoarse + (DIM+1) * connec[i] + DIM;
                    PetscInt dof_j = (DIM+1) * numNodesCoarse + (DIM+1) * numNodesFine + connecL[j] + k;
                    MatSetValues(A,1,&dof_i,1,&dof_j,&localMV_mat(DIM*i+k,DIM*nElNodes+j),ADD_VALUES);
                };
                //ARLEQUIN STABILIZATION
                if (fabs(ArlequinA2(DIM*i+k,DIM*j+k)) >= 1.e-15){
                    PetscInt d_i = (DIM+1) * numNodesCoarse + (DIM+1) * numNodesFine + DIM * connecL[i] + k;
                    PetscInt d_j = (DIM+1) * numNodesCoarse + (DIM+1) * connec[j] + k;
                    MatSetValues(A,1,&d_i,1,&d_j,&ArlequinA2(DIM*i+k,DIM*j+k),ADD_VALUES);
                };
                if (fabs(ArlequinA2(DIM*nElNodes+i,DIM*j+k)) >= 1.e-15){
                    PetscInt dof_i = (DIM+1) * numNodesCoarse + (DIM+1) * numNodesFine + DIM * connecL[i] + k;
                    PetscInt dof_j = (DIM+1) * numNodesCoarse + (DIM+1) * connec[j] + DIM;
                    MatSetValues(A,1,&dof_i,1,&dof_j,&ArlequinA2(DIM*nElNodes+i,DIM*j+k),ADD_VALUES);
                };
                if (fabs(ArlequinA1(DIM*i+k,DIM*j+k)) >= 1.e-15){
                    PetscInt dof_i = (DIM+1) * numNodesCoarse + (DIM+1) * numNodesFine + DIM * connecL[i] + k;
                    PetscInt dof_j = (DIM+1) * numNodesCoarse + (DIM+1) * numNodesFine + DIM * connecL[j] + k;
                    MatSetValues(A,1,&dof_i,1,&dof_j,&ArlequinA1(DIM*i+k,DIM*j+k),ADD_VALUES);
                };      
            }                      
        };
        //RHS VECTOR
        for (int k = 0; k < DIM; k++){
            //COUPLING OPERATOR
            PetscInt dof_i = (DIM+1) * numNodesCoarse + (DIM+1) * numNodesFine + DIM * connecL[i] + k;
            VecSetValues(b,1,&dof_i,&Rhs2[DIM*i+k],ADD_VALUES);

            dof_i = (DIM+1) * numNodesCoarse + (DIM+1) * connec[i] + k;
            VecSetValues(b,1,&dof_i,&rhsLagMult2[DIM*i+k],ADD_VALUES);

            //SUPG STABILIZATION
            dof_i = (DIM+1) * numNodesCoarse + (DIM+1) * connec[i] + k;
            VecSetValues(b,1,&dof_i,&localMV_vec[DIM*i+k],ADD_VALUES);

            //ARLEQUIN STABILIZATION
            dof_i = (DIM+1) * numNodesCoarse + (DIM+1) * numNodesFine + DIM*connecL[i] + k;
            VecSetValues(b,1,&dof_i,&RhsArlequin2[DIM*i+k],ADD_VALUES);
        }
        //PSPG STABILIZATION
        PetscInt dof_i = (DIM+1) * numNodesCoarse + (DIM+1) * connec[i] + DIM;
        VecSetValues(b,1,&dof_i,&localMV_vec[DIM*nElNodes+i],ADD_VALUES);
    };


    return;
};

template<int DIM, int DEG>
void Arlequin<DIM,DEG>::setVecValuesLagMultFineFine(VecDouble &Rhs2, VecDouble &rhsLagMult2,
                                                    VecDouble &localMV_vec, VecDouble &RhsArlequin2,
                                                    VecInt &connec, VecInt &connecL){

    double &alpha_f = parametersFine -> getAlphaF();
    double &alpha_m = parametersFine -> getAlphaM();
    double &gamma = parametersFine -> getGamma();
    
    double integ = alpha_f * gamma * dTime;
    //Disperse local contributions into the global matrix
    for (int i = 0; i < nElNodes; i++){
        //RHS VECTOR
        for (int k = 0; k < DIM; k++){
            //COUPLING OPERATOR
            PetscInt dof_i = (DIM+1) * numNodesCoarse + (DIM+1) * numNodesFine + DIM * connecL[i] + k;
            VecSetValues(b,1,&dof_i,&Rhs2[DIM*i+k],ADD_VALUES);

            dof_i = (DIM+1) * numNodesCoarse + (DIM+1) * connec[i] + k;
            VecSetValues(b,1,&dof_i,&rhsLagMult2[DIM*i+k],ADD_VALUES);

            //SUPG STABILIZATION
            dof_i = (DIM+1) * numNodesCoarse + (DIM+1) * connec[i] + k;
            VecSetValues(b,1,&dof_i,&localMV_vec[DIM*i+k],ADD_VALUES);

            //ARLEQUIN STABILIZATION
            dof_i = (DIM+1) * numNodesCoarse + (DIM+1) * numNodesFine + DIM*connecL[i] + k;
            VecSetValues(b,1,&dof_i,&RhsArlequin2[DIM*i+k],ADD_VALUES);
        }
        //PSPG STABILIZATION
        PetscInt dof_i = (DIM+1) * numNodesCoarse + (DIM+1) * connec[i] + DIM;
        VecSetValues(b,1,&dof_i,&localMV_vec[DIM*nElNodes+i],ADD_VALUES);
    };


    return;
};

template<int DIM, int DEG>
void Arlequin<DIM,DEG>::setMatValuesLagMultFineFine(MatrixDouble &Ajac2, MatrixDouble &localMV_mat, 
                                                    MatrixDouble &ArlequinA1, MatrixDouble &ArlequinA2,
                                                    VecInt &connec, VecInt &connecL){

    double &alpha_f = parametersFine -> getAlphaF();
    double &alpha_m = parametersFine -> getAlphaM();
    double &gamma = parametersFine -> getGamma();
    
    double integ = alpha_f * gamma * dTime;
    //Disperse local contributions into the global matrix
    for (int i = 0; i < nElNodes; i++){
        for (int j = 0; j < nElNodes; j++){
            for (int k = 0; k < DIM; k++){
                for (int l = 0; l < DIM; l++){
                    //COUPLING OPERATOR
                    if (fabs(Ajac2(DIM*i+k,DIM*j+l)) >= 1.e-15){
                        PetscInt d_i = (DIM+1) * numNodesCoarse + (DIM+1) * numNodesFine + DIM * connecL[i] + k;
                        PetscInt d_j = (DIM+1) * numNodesCoarse + (DIM+1) * connec[j] + l;
                        double value = Ajac2(DIM*i+k,DIM*j+l)*integ;
                        MatSetValues(A,1,&d_i,1,&d_j,&value,ADD_VALUES);
                        MatSetValues(A,1,&d_j,1,&d_i,&Ajac2(DIM*i+k,DIM*j+l),ADD_VALUES);
                    };
                }
                //SUPG STABILIZATION
                if (fabs(localMV_mat(DIM*i+k,DIM*j+k)) >= 1.e-15){
                    PetscInt d_i = (DIM+1) * numNodesCoarse + (DIM+1) * numNodesFine + DIM*connecL[i] + k;
                    PetscInt d_j = (DIM+1) * numNodesCoarse + (DIM+1) * connec[j] + k;
                    MatSetValues(A,1,&d_j,1,&d_i,&localMV_mat(2*i  ,2*j  ),ADD_VALUES);
                };
                //PSPG STABILIZATION
                if (fabs(localMV_mat(DIM*i+k,DIM*nElNodes+j)) >= 1.e-15){
                    PetscInt dof_i = (DIM+1) * numNodesCoarse + (DIM+1) * connec[i] + DIM;
                    PetscInt dof_j = (DIM+1) * numNodesCoarse + (DIM+1) * numNodesFine + connecL[j] + k;
                    MatSetValues(A,1,&dof_i,1,&dof_j,&localMV_mat(DIM*i+k,DIM*nElNodes+j),ADD_VALUES);
                };
                //ARLEQUIN STABILIZATION
                if (fabs(ArlequinA2(DIM*i+k,DIM*j+k)) >= 1.e-15){
                    PetscInt d_i = (DIM+1) * numNodesCoarse + (DIM+1) * numNodesFine + DIM * connecL[i] + k;
                    PetscInt d_j = (DIM+1) * numNodesCoarse + (DIM+1) * connec[j] + k;
                    MatSetValues(A,1,&d_i,1,&d_j,&ArlequinA2(DIM*i+k,DIM*j+k),ADD_VALUES);
                };
                if (fabs(ArlequinA2(DIM*nElNodes+i,DIM*j+k)) >= 1.e-15){
                    PetscInt dof_i = (DIM+1) * numNodesCoarse + (DIM+1) * numNodesFine + DIM * connecL[i] + k;
                    PetscInt dof_j = (DIM+1) * numNodesCoarse + (DIM+1) * connec[j] + DIM;
                    MatSetValues(A,1,&dof_i,1,&dof_j,&ArlequinA2(DIM*nElNodes+i,DIM*j+k),ADD_VALUES);
                };
                if (fabs(ArlequinA1(DIM*i+k,DIM*j+k)) >= 1.e-15){
                    PetscInt dof_i = (DIM+1) * numNodesCoarse + (DIM+1) * numNodesFine + DIM * connecL[i] + k;
                    PetscInt dof_j = (DIM+1) * numNodesCoarse + (DIM+1) * numNodesFine + DIM * connecL[j] + k;
                    MatSetValues(A,1,&dof_i,1,&dof_j,&ArlequinA1(DIM*i+k,DIM*j+k),ADD_VALUES);
                };      
            }                      
        };
    };
    return;
};

//------------------------------------------------------------------------------
//----------------COMPUTE ARLEQUIN COUPLED NAVIER-STOKES PROBLEM----------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Arlequin<DIM,DEG>::setMatVecValuesLagMultFineFinePoisson(MatrixDouble &Ajac2, MatrixDouble &localMV_mat, 
                                                              MatrixDouble &ArlequinA1, MatrixDouble &ArlequinA2, 
                                                              VecDouble &Rhs2, VecDouble &rhsLagMult2,
                                                              VecDouble &localMV_vec, VecDouble &RhsArlequin2,
                                                              VecInt &connec, VecInt &connecL){

    //Disperse local contributions into the global matrix
    for (int i = 0; i < nElNodes; i++){
        for (int j = 0; j < nElNodes; j++){
            //COUPLING OPERATOR
            if (fabs(Ajac2(i,j)) >= 1.e-15){
                PetscInt d_i = numNodesCoarse + numNodesFine + connecL[i];
                PetscInt d_j = numNodesCoarse + connec[j];
                MatSetValues(A,1,&d_i,1,&d_j,&Ajac2(i,j),ADD_VALUES);
                MatSetValues(A,1,&d_j,1,&d_i,&Ajac2(i,j),ADD_VALUES);
            };
            // ARLEQUIN STABILIZATION
            if (fabs(ArlequinA2(i,j)) >= 1.e-15){
                PetscInt dof_i = numNodesCoarse + numNodesFine + connecL[i];
                PetscInt dof_j = numNodesCoarse + connec[j];
                MatSetValues(A,1,&dof_i,1,&dof_j,&ArlequinA2(i,j),ADD_VALUES);
            };
            if (fabs(ArlequinA1(i,j)) >= 1.e-15){
                PetscInt dof_i = numNodesCoarse + numNodesFine + connecL[i];
                PetscInt dof_j = numNodesCoarse + numNodesFine + connecL[j];
                MatSetValues(A,1,&dof_i,1,&dof_j,&ArlequinA1(i,j),ADD_VALUES);
            };                   
        };
        //RHS VECTOR
        //COUPLING OPERATOR
        PetscInt dof_i = numNodesCoarse + numNodesFine + connecL[i];
        VecSetValues(b,1,&dof_i,&Rhs2[i],ADD_VALUES);

        dof_i = numNodesCoarse + connec[i];
        VecSetValues(b,1,&dof_i,&rhsLagMult2[i],ADD_VALUES);

        //ARLEQUIN STABILIZATION
        dof_i = numNodesCoarse + numNodesFine + connecL[i];
        VecSetValues(b,1,&dof_i,&RhsArlequin2[i],ADD_VALUES);
    };


    return;
};

//------------------------------------------------------------------------------
//----------------COMPUTE ARLEQUIN COUPLED NAVIER-STOKES PROBLEM----------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Arlequin<DIM,DEG>::setMatVecValuesLagMultFineCoarse(MatrixDouble &Ajac2, MatrixDouble &localMV_mat, 
                                                   MatrixDouble &ArlequinA1, MatrixDouble &ArlequinA2, 
                                                   VecDouble &Rhs2, VecDouble &rhsLagMult2,
                                                   VecDouble &localMV_vec, VecDouble &RhsArlequin2,
                                                   VecInt &connecC, VecInt &connecL){

    double &alpha_f = parametersFine -> getAlphaF();
    double &alpha_m = parametersFine -> getAlphaM();
    double &gamma = parametersFine -> getGamma();
    double integ = alpha_f * gamma * dTime;

    //Disperse local contribution into the global matrix
    for (int i = 0; i < nElNodes; i++){
        for (int j = 0; j < nElNodes; j++){
            for (int k = 0; k < DIM; k++){
                for (int l = 0; l < DIM; l++){
                    //COUPLING OPERATOR
                    if (fabs(Ajac2(DIM*i+k,DIM*j+l)) >= 1.e-15){
                        PetscInt dof_i = (DIM+1) * numNodesCoarse + (DIM+1) * numNodesFine + DIM * connecL[i] + k;
                        PetscInt dof_j = (DIM+1) * connecC[j] + l;
                        double value = Ajac2(DIM*i+k,DIM*j+l) * integ;
                        MatSetValues(A,1,&dof_i,1,&dof_j,&value,ADD_VALUES);
                        MatSetValues(A,1,&dof_j,1,&dof_i,&Ajac2(DIM*i+k,DIM*j+l),ADD_VALUES);
                    }
                }
                //SUPG STABILIZATION
                if (fabs(localMV_mat(DIM*i+k,DIM*j+k)) >= 1.e-15){
                    PetscInt dof_i = (DIM+1) * numNodesCoarse + (DIM+1) * numNodesFine + DIM * connecL[i] + k;
                    PetscInt dof_j = (DIM+1) * connecC[j] + k;
                    MatSetValues(A,1,&dof_j,1,&dof_i,&localMV_mat(DIM*i+k,DIM*j+k),ADD_VALUES);
                };
                //PSPG STABILIZATION
                if (fabs(localMV_mat(DIM*i+k,DIM*nElNodes+j)) >= 1.e-15){
                    PetscInt dof_i = (DIM+1) * connecC[i] + DIM;
                    PetscInt dof_j = (DIM+1) * numNodesCoarse + DIM * numNodesFine + connecL[j] + k;
                    MatSetValues(A,1,&dof_i,1,&dof_j,&localMV_mat(DIM*i+k,DIM*nElNodes+j),ADD_VALUES);
                };
                //ARLEQUIN STABILIZATION
                if (fabs(ArlequinA2(DIM*i+k,DIM*j+k)) >= 1.e-15){
                    PetscInt d_i = (DIM+1) * numNodesCoarse + (DIM+1) * numNodesFine + DIM * connecL[i] + k;
                    PetscInt d_j = (DIM+1) * connecC[j] + k;
                    MatSetValues(A,1,&d_i,1,&d_j,&ArlequinA2(2*i  ,2*j  ),ADD_VALUES);
                };
                if (fabs(ArlequinA2(DIM*nElNodes+i,DIM*j+k)) >= 1.e-15){
                    PetscInt dof_i = (DIM+1) * numNodesCoarse + (DIM+1) * numNodesFine + DIM * connecL[i] + 1;
                    PetscInt dof_j = (DIM+1) * connecC[j] + DIM;
                    MatSetValues(A,1,&dof_i,1,&dof_j,&ArlequinA2(DIM*nElNodes+i,DIM*j+k),ADD_VALUES);
                };
                if (fabs(ArlequinA1(DIM*i+k,DIM*j+k)) >= 1.e-15){
                    PetscInt dof_i = (DIM+1) * numNodesCoarse + (DIM+1)*numNodesFine + DIM * connecL[i] + k;
                    PetscInt dof_j = (DIM+1) * numNodesCoarse + (DIM+1)*numNodesFine + DIM * connecL[j] + k;
                    MatSetValues(A,1,&dof_i,1,&dof_j,&ArlequinA1(DIM*i+k,DIM*j+k),ADD_VALUES);
                };
            }
            // if (fabs(Ajac2(12+i,12+j)) >= 1.e-15){
            //     int dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL[i];
            //     int dof_j = 3*numNodesCoarse + 3*numNodesFine + 2*connecL[j];
            //     dof_i++; dof_j++;
            //     ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&Ajac2(12+i,12+j),ADD_VALUES);
            // };
        };
        for (int k = 0; k < DIM; k++){
            //RHS VECTOR
            //COUPLING OPERATOR
            PetscInt d_i = (DIM+1) * numNodesCoarse + (DIM+1) * numNodesFine + DIM * connecL[i] + k;
            VecSetValues(b,1,&d_i,&Rhs2[DIM*i+k],ADD_VALUES);

            PetscInt dof_i = (DIM+1) * connecC[i] + k;
            VecSetValues(b,1,&dof_i,&rhsLagMult2[DIM*i+k],ADD_VALUES);

            //SUPG STABILIZATION
            dof_i = (DIM+1) * connecC[i] + k;
            VecSetValues(b,1,&dof_i,&localMV_vec[DIM*i+k],ADD_VALUES);

            dof_i = (DIM+1) * numNodesCoarse + (DIM+1) * numNodesFine + DIM * connecL[i] + k;
            VecSetValues(b,1,&dof_i,&RhsArlequin2[DIM*i+k],ADD_VALUES);
        }
        //PSPG STABILIZATION
        PetscInt dof_i = (DIM+1) * connecC[i] + DIM;
        VecSetValues(b,1,&dof_i,&localMV_vec[DIM*nElNodes+i],ADD_VALUES);
    };
    return;
}


template<int DIM, int DEG>
void Arlequin<DIM,DEG>::setVecValuesLagMultFineCoarse(VecDouble &Rhs2, VecDouble &rhsLagMult2,
                                                      VecDouble &localMV_vec, VecDouble &RhsArlequin2,
                                                      VecInt &connecC, VecInt &connecL){

    double &alpha_f = parametersFine -> getAlphaF();
    double &alpha_m = parametersFine -> getAlphaM();
    double &gamma = parametersFine -> getGamma();
    double integ = alpha_f * gamma * dTime;

    //Disperse local contribution into the global matrix
    for (int i = 0; i < nElNodes; i++){
        for (int k = 0; k < DIM; k++){
            //RHS VECTOR
            //COUPLING OPERATOR
            PetscInt d_i = (DIM+1) * numNodesCoarse + (DIM+1) * numNodesFine + DIM * connecL[i] + k;
            VecSetValues(b,1,&d_i,&Rhs2[DIM*i+k],ADD_VALUES);

            PetscInt dof_i = (DIM+1) * connecC[i] + k;
            VecSetValues(b,1,&dof_i,&rhsLagMult2[DIM*i+k],ADD_VALUES);

            //SUPG STABILIZATION
            dof_i = (DIM+1) * connecC[i] + k;
            VecSetValues(b,1,&dof_i,&localMV_vec[DIM*i+k],ADD_VALUES);

            dof_i = (DIM+1) * numNodesCoarse + (DIM+1) * numNodesFine + DIM * connecL[i] + k;
            VecSetValues(b,1,&dof_i,&RhsArlequin2[DIM*i+k],ADD_VALUES);
        }
        //PSPG STABILIZATION
        PetscInt dof_i = (DIM+1) * connecC[i] + DIM;
        VecSetValues(b,1,&dof_i,&localMV_vec[DIM*nElNodes+i],ADD_VALUES);
    };
    return;
}


template<int DIM, int DEG>
void Arlequin<DIM,DEG>::setMatValuesLagMultFineCoarse(MatrixDouble &Ajac2, MatrixDouble &localMV_mat, 
                                                      MatrixDouble &ArlequinA1, MatrixDouble &ArlequinA2, 
                                                      VecInt &connecC, VecInt &connecL){

    double &alpha_f = parametersFine -> getAlphaF();
    double &alpha_m = parametersFine -> getAlphaM();
    double &gamma = parametersFine -> getGamma();
    double integ = alpha_f * gamma * dTime;

    //Disperse local contribution into the global matrix
    for (int i = 0; i < nElNodes; i++){
        for (int j = 0; j < nElNodes; j++){
            for (int k = 0; k < DIM; k++){
                for (int l = 0; l < DIM; l++){
                    //COUPLING OPERATOR
                    if (fabs(Ajac2(DIM*i+k,DIM*j+l)) >= 1.e-15){
                        PetscInt dof_i = (DIM+1) * numNodesCoarse + (DIM+1) * numNodesFine + DIM * connecL[i] + k;
                        PetscInt dof_j = (DIM+1) * connecC[j] + l;
                        double value = Ajac2(DIM*i+k,DIM*j+l) * integ;
                        MatSetValues(A,1,&dof_i,1,&dof_j,&value,ADD_VALUES);
                        MatSetValues(A,1,&dof_j,1,&dof_i,&Ajac2(DIM*i+k,DIM*j+l),ADD_VALUES);
                    }
                }
                //SUPG STABILIZATION
                if (fabs(localMV_mat(DIM*i+k,DIM*j+k)) >= 1.e-15){
                    PetscInt dof_i = (DIM+1) * numNodesCoarse + (DIM+1) * numNodesFine + DIM * connecL[i] + k;
                    PetscInt dof_j = (DIM+1) * connecC[j] + k;
                    MatSetValues(A,1,&dof_j,1,&dof_i,&localMV_mat(DIM*i+k,DIM*j+k),ADD_VALUES);
                };
                //PSPG STABILIZATION
                if (fabs(localMV_mat(DIM*i+k,DIM*nElNodes+j)) >= 1.e-15){
                    PetscInt dof_i = (DIM+1) * connecC[i] + DIM;
                    PetscInt dof_j = (DIM+1) * numNodesCoarse + DIM * numNodesFine + connecL[j] + k;
                    MatSetValues(A,1,&dof_i,1,&dof_j,&localMV_mat(DIM*i+k,DIM*nElNodes+j),ADD_VALUES);
                };
                //ARLEQUIN STABILIZATION
                if (fabs(ArlequinA2(DIM*i+k,DIM*j+k)) >= 1.e-15){
                    PetscInt d_i = (DIM+1) * numNodesCoarse + (DIM+1) * numNodesFine + DIM * connecL[i] + k;
                    PetscInt d_j = (DIM+1) * connecC[j] + k;
                    MatSetValues(A,1,&d_i,1,&d_j,&ArlequinA2(2*i  ,2*j  ),ADD_VALUES);
                };
                if (fabs(ArlequinA2(DIM*nElNodes+i,DIM*j+k)) >= 1.e-15){
                    PetscInt dof_i = (DIM+1) * numNodesCoarse + (DIM+1) * numNodesFine + DIM * connecL[i] + 1;
                    PetscInt dof_j = (DIM+1) * connecC[j] + DIM;
                    MatSetValues(A,1,&dof_i,1,&dof_j,&ArlequinA2(DIM*nElNodes+i,DIM*j+k),ADD_VALUES);
                };
                if (fabs(ArlequinA1(DIM*i+k,DIM*j+k)) >= 1.e-15){
                    PetscInt dof_i = (DIM+1) * numNodesCoarse + (DIM+1)*numNodesFine + DIM * connecL[i] + k;
                    PetscInt dof_j = (DIM+1) * numNodesCoarse + (DIM+1)*numNodesFine + DIM * connecL[j] + k;
                    MatSetValues(A,1,&dof_i,1,&dof_j,&ArlequinA1(DIM*i+k,DIM*j+k),ADD_VALUES);
                };
            }
            // if (fabs(Ajac2(12+i,12+j)) >= 1.e-15){
            //     int dof_i = 3*numNodesCoarse + 3*numNodesFine + 2*connecL[i];
            //     int dof_j = 3*numNodesCoarse + 3*numNodesFine + 2*connecL[j];
            //     dof_i++; dof_j++;
            //     ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&Ajac2(12+i,12+j),ADD_VALUES);
            // };
        };
    };
    return;
}


//------------------------------------------------------------------------------
//----------------COMPUTE ARLEQUIN COUPLED NAVIER-STOKES PROBLEM----------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Arlequin<DIM,DEG>::setMatVecValuesLagMultFineCoarsePoisson(MatrixDouble &Ajac2, MatrixDouble &localMV_mat, 
                                                                MatrixDouble &ArlequinA1, MatrixDouble &ArlequinA2, 
                                                                VecDouble &Rhs2, VecDouble &rhsLagMult2,
                                                                VecDouble &localMV_vec, VecDouble &RhsArlequin2,
                                                                VecInt &connecC, VecInt &connecL){

    //Disperse local contribution into the global matrix
    for (int i = 0; i < nElNodes; i++){
        for (int j = 0; j < nElNodes; j++){
            //COUPLING OPERATOR
            if (fabs(Ajac2(i,j)) >= 1.e-15){
                PetscInt dof_i = numNodesCoarse + numNodesFine + connecL[i];
                PetscInt dof_j = connecC[j];
                MatSetValues(A,1,&dof_i,1,&dof_j,&Ajac2(i,j),ADD_VALUES);
                MatSetValues(A,1,&dof_j,1,&dof_i,&Ajac2(i,j),ADD_VALUES);
            }
             
            //ARLEQUIN STABILIZATION
            if (fabs(ArlequinA2(i,j)) >= 1.e-15){
                PetscInt dof_i = numNodesCoarse + numNodesFine + connecL[i];
                PetscInt dof_j = connecC[j];
                MatSetValues(A,1,&dof_i,1,&dof_j,&ArlequinA2(i,j),ADD_VALUES);
            };
            if (fabs(ArlequinA1(i,j)) >= 1.e-15){
                PetscInt dof_i = numNodesCoarse + numNodesFine + connecL[i];
                PetscInt dof_j = numNodesCoarse + numNodesFine + connecL[j];
                MatSetValues(A,1,&dof_i,1,&dof_j,&ArlequinA1(i,j),ADD_VALUES);
            };
        };
        
        //RHS VECTOR
        //COUPLING OPERATOR
        PetscInt d_i = numNodesCoarse + numNodesFine + connecL[i];
        VecSetValues(b,1,&d_i,&Rhs2[i],ADD_VALUES);

        PetscInt dof_i = connecC[i];
        VecSetValues(b,1,&dof_i,&rhsLagMult2[i],ADD_VALUES);


        dof_i = numNodesCoarse + numNodesFine + connecL[i];
        VecSetValues(b,1,&dof_i,&RhsArlequin2[i],ADD_VALUES);
    };
    return;
}

//------------------------------------------------------------------------------
//----------------COMPUTE ARLEQUIN COUPLED NAVIER-STOKES PROBLEM----------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Arlequin<DIM,DEG>::assembleArlequinSystem(){

    //Coarse mesh
    for (int jel = 0; jel < numElemCoarse; jel++){   
        if (domDecompCoarse.first[jel] == rank) {
            //Compute Element matrix
            MatrixDouble matrix(nLocDOF,nLocDOF);
            matrix.setZero();
            VecDouble rhs(nLocDOF);
            rhs.setZero();

            elementsCoarse_[jel] -> getTransientNavierStokes(matrix,rhs);

            setMatVecValuesCoarseModel(matrix,rhs,elementsCoarse_[jel] -> getConnectivity());

        };
    };

    //Fine mesh
    for (int jel = 0; jel < numElemFine; jel++){           
        if (domDecompFine.first[jel] == rank) {
            //Compute Element matrix                    
            MatrixDouble matrix(nLocDOF,nLocDOF);
            matrix.setZero();
            VecDouble rhs(nLocDOF);
            rhs.setZero();

            elementsFine_[jel] -> getTransientNavierStokes(matrix,rhs);
            
            setMatVecValuesFineModel(matrix,rhs,elementsFine_[jel] -> getConnectivity());

        };                
    };
      

    //Lagrange Multipliers
    for (int l=0; l< numElemGlueZoneFine; l++){
        int jel = elementsGlueZoneFine_[l];
        if (domDecompFine.first[jel] == rank) {
            
            VecInt connecC;
            VecInt connec = elementsFine_[jel] -> getConnectivity();
            VecInt connecL = glueZoneFine_[l] -> getConnectivity();
            //FINE MESH
            //Matrices
            MatrixDouble Ajac2(nLocDOF,nLocDOF);
            MatrixDouble localMV_mat(nLocDOF,nLocDOF);
            MatrixDouble ArlequinA1(nLocDOF,nLocDOF);
            MatrixDouble ArlequinA2(nLocDOF,nLocDOF);
            Ajac2.setZero();
            localMV_mat.setZero();
            ArlequinA1.setZero();
            ArlequinA2.setZero();

            //Vectors
            VecDouble rhsLagMult2(nLocDOF);
            VecDouble Rhs2(nLocDOF);
            VecDouble localMV_vec(nLocDOF);
            VecDouble RhsArlequin2(nLocDOF);
            rhsLagMult2.setZero();
            Rhs2.setZero();
            localMV_vec.setZero();
            RhsArlequin2.setZero();

            //Stabilization 
            double tArlequin;
            VecDouble Ml1(nLocDOF), t1(nLocDOF), j1(nLocDOF), 
                      k1(nLocDOF), p1(nLocDOF);
            Ml1.setZero(); t1.setZero(); j1.setZero(); k1.setZero(); p1.setZero();

            // FINE MESH
            //Computes element matrix

            elementsFine_[jel] -> getLagrangeMultipliersSameMesh(Ajac2, rhsLagMult2, Rhs2);
            
            if (fArlequinStab != ArlequinStabType::ENoStab){
                //PSPG and SUPG stabilizations
                elementsFine_[jel] -> getLagrangeMultipliersSUPG_PSPG_SameMesh(localMV_mat,localMV_vec);
            
                //Arlequin Stabilization
                elementsFine_[jel] -> getLagrangeMultipliersArlequinSameMesh(ArlequinA1, ArlequinA2, RhsArlequin2, Ml1, t1, j1, 
                                                                           k1, p1);

            }
            // stabilizeArlequinNew(Ml1, t1, j1, k1, p1, tArlequin);
            // ArlequinA1 *= tArlequin;
            // ArlequinA2 *= tArlequin;
            // RhsArlequin2 *= tArlequin;
            
            setMatVecValuesLagMultFineFine(Ajac2,localMV_mat,ArlequinA1,ArlequinA2, 
                                           Rhs2,rhsLagMult2,localMV_vec,RhsArlequin2,
                                           elementsFine_[jel] -> getConnectivity(),
                                           glueZoneFine_[l] -> getConnectivity());
            
            //COAESE MESH
            //Counts number of coarse mesh intersecting the fine element
            int numberIntPoints = elementsFine_[jel] -> 
                getNumberOfIntegrationPoints();
            int aux;
            
            std::vector<int> ele, diffElem;
            ele.clear();
            diffElem.clear();

            ele.reserve(3);
            for (int i=0; i<numberIntPoints; i++){
                aux = elementsFine_[jel] -> 
                    getIntegPointCorrespondenceElement(i);
                ele.push_back(aux);
                //std::cout << "Num elem inters " << jel << " " << aux << std::endl;
            };
            
            int numElemIntersect = 1;
            int flag = 0;
            diffElem.push_back(ele[0]);
            
            for (int i = 1; i<numberIntPoints; i++){
                flag = 0;
                for (int j = 0; j<numElemIntersect; j++){
                    if (ele[i] == diffElem[j]) {
                        break;
                    }else{
                        flag++;
                    };
                    if(flag == numElemIntersect){
                        numElemIntersect++;
                        diffElem.push_back(ele[i]);
                    };
                };
            };
            //Compute the Lagrange Multiplier element matrix
            for (int ielem = 0; ielem < numElemIntersect; ielem++){
                
                int iElemCoarse = diffElem[ielem];
                double pspg = 0;//elementsCoarse_[iElemCoarse] -> getPSPG();
                VecDouble press_(nElNodes), velX_(nElNodes), velY_(nElNodes), velXPrev_(nElNodes), velYPrev_(nElNodes), acelX_(nElNodes), acelY_(nElNodes), acelXPrev_(nElNodes), acelYPrev_(nElNodes);

                connecC = elementsCoarse_[iElemCoarse] -> getConnectivity();

                for (int k = 0; k < nElNodes; k++){
                    press_[k] = (*nodesCoarse_)[connecC[k]] -> getPressure();
                    velX_[k] = (*nodesCoarse_)[connecC[k]] -> getVelocity(0);
                    velY_[k] = (*nodesCoarse_)[connecC[k]] -> getVelocity(1);
                    velXPrev_[k] = (*nodesCoarse_)[connecC[k]] -> getPreviousVelocity(0);
                    velYPrev_[k] = (*nodesCoarse_)[connecC[k]] -> getPreviousVelocity(1);
                    acelX_[k] = (*nodesCoarse_)[connecC[k]] -> getAcceleration(0);
                    acelY_[k] = (*nodesCoarse_)[connecC[k]] -> getAcceleration(1);
                    acelXPrev_[k] = (*nodesCoarse_)[connecC[k]] -> getPreviousAcceleration(0);
                    acelYPrev_[k] = (*nodesCoarse_)[connecC[k]] -> getPreviousAcceleration(1);
                }
                
                Ajac2.setZero();
                localMV_mat.setZero();
                ArlequinA1.setZero();
                ArlequinA2.setZero();
                
                //Vectors
                rhsLagMult2.setZero();
                Rhs2.setZero();
                localMV_vec.setZero();
                RhsArlequin2.setZero();
                
                elementsFine_[jel] -> getLagrangeMultipliersDifferentMesh(iElemCoarse,pspg,press_,velX_,velY_,velXPrev_,velYPrev_,Ajac2,rhsLagMult2,Rhs2);

                if (fArlequinStab != ArlequinStabType::ENoStab){
                    elementsFine_[jel] -> getLagrangeMultipliersSUPG_PSPG_DifferentMesh(iElemCoarse,pspg,press_,velX_,velY_,localMV_mat,localMV_vec);

                    elementsFine_[jel] -> getLagrangeMultipliersArlequinDifferentMesh(iElemCoarse,pspg,press_,velX_,velY_,acelX_,acelY_,
                                                                                      acelXPrev_,acelYPrev_,ArlequinA1,ArlequinA2,RhsArlequin2);
                }

                // ArlequinA1 *= tArlequin;
                // ArlequinA2 *= tArlequin;
                // RhsArlequin2 *= tArlequin;

                setMatVecValuesLagMultFineCoarse(Ajac2, localMV_mat, ArlequinA1, ArlequinA2, 
                                                 Rhs2, rhsLagMult2, localMV_vec, RhsArlequin2,
                                                 elementsCoarse_[iElemCoarse] -> getConnectivity(), 
                                                 glueZoneFine_[l] -> getConnectivity());
                
            }; //Number of intersections
        }; // if element belongs to the glue zone
    }; // Glue zone

    return;
}






//------------------------------------------------------------------------------
//----------------COMPUTE ARLEQUIN COUPLED NAVIER-STOKES PROBLEM----------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Arlequin<DIM,DEG>::assembleArlequinSystemPoisson(){

    //Coarse mesh
    for (int jel = 0; jel < numElemCoarse; jel++){   
        if (domDecompCoarse.first[jel] == rank) {
            //Compute Element matrix
            MatrixDouble matrix(nElNodes,nElNodes);
            matrix.setZero();
            VecDouble rhs(nElNodes);
            rhs.setZero();

            elementsCoarse_[jel] -> getPoisson(matrix,rhs);

            setMatVecValuesCoarseModelPoisson(matrix,rhs,elementsCoarse_[jel] -> getConnectivity());

        };
    };

    //Fine mesh
    for (int jel = 0; jel < numElemFine; jel++){           
        if (domDecompFine.first[jel] == rank) {
            //Compute Element matrix                    
            MatrixDouble matrix(nElNodes,nElNodes);
            matrix.setZero();
            VecDouble rhs(nElNodes);
            rhs.setZero();

            elementsFine_[jel] -> getPoisson(matrix,rhs);
            
            setMatVecValuesFineModelPoisson(matrix,rhs,elementsFine_[jel] -> getConnectivity());

        };                
    };
      

    //Lagrange Multipliers
    for (int l=0; l< numElemGlueZoneFine; l++){
        int jel = elementsGlueZoneFine_[l];
        if (domDecompFine.first[jel] == rank) {
            
            VecInt connecC;
            VecInt connec = elementsFine_[jel] -> getConnectivity();
            VecInt connecL = glueZoneFine_[l] -> getConnectivity();
            //FINE MESH
            //Matrices
            MatrixDouble matC0(nElNodes,nElNodes),matC1(nElNodes,nElNodes);
            MatrixDouble matA0(nElNodes,nElNodes),matA1(nElNodes,nElNodes);
            MatrixDouble localMV_mat(nElNodes,nElNodes);
            MatrixDouble matE0(nElNodes,nElNodes),matE1(nElNodes,nElNodes);
            matC1.setZero();
            matA1.setZero();
            localMV_mat.setZero();
            matE1.setZero();

            double tARLQ0_, tARLQ1_;

            //Vectors
            VecDouble vecC0(nElNodes),vecC1(nElNodes);
            VecDouble vecU0(nElNodes),vecU1(nElNodes);
            VecDouble RhsA0(nElNodes),RhsA1(nElNodes);
            VecDouble localMV_vec(nElNodes);
            VecDouble vecE0(nElNodes),vecE1(nElNodes);
            vecC1.setZero();
            vecU1.setZero();
            localMV_vec.setZero();
            vecE1.setZero();

            // FINE MESH
            //Computes element matrix
            elementsFine_[jel] -> getLagrangeMultipliersSameMeshPoisson(matC1, vecC1, vecU1);
            
            if (fArlequinStab != ArlequinStabType::ENoStab){
                //Arlequin Stabilization
                elementsFine_[jel] -> getLagrangeMultipliersArlequinSameMeshPoisson(matE1, matA1, vecE1);
            }

            
            
            //COAESE MESH
            //Counts number of coarse mesh intersecting the fine element
            int numberIntPoints = elementsFine_[jel] -> 
                getNumberOfIntegrationPoints();
            int aux;
            
            std::vector<int> ele, diffElem;
            ele.clear();
            diffElem.clear();

            ele.reserve(3);
            for (int i=0; i<numberIntPoints; i++){
                aux = elementsFine_[jel] -> 
                    getIntegPointCorrespondenceElement(i);
                ele.push_back(aux);
                //std::cout << "Num elem inters " << jel << " " << aux << std::endl;
            };
            
            int numElemIntersect = 1;
            int flag = 0;
            diffElem.push_back(ele[0]);
            
            for (int i = 1; i<numberIntPoints; i++){
                flag = 0;
                for (int j = 0; j<numElemIntersect; j++){
                    if (ele[i] == diffElem[j]) {
                        break;
                    }else{
                        flag++;
                    };
                    if(flag == numElemIntersect){
                        numElemIntersect++;
                        diffElem.push_back(ele[i]);
                    };
                };
            };
            //Compute the Lagrange Multiplier element matrix
            for (int ielem = 0; ielem < numElemIntersect; ielem++){
                
                int iElemCoarse = diffElem[ielem];
                double pspg = 0;//elementsCoarse_[iElemCoarse] -> getPSPG();
                VecDouble press_(nElNodes), velX_(nElNodes), velY_(nElNodes), velXPrev_(nElNodes), velYPrev_(nElNodes);

                connecC = elementsCoarse_[iElemCoarse] -> getConnectivity();

                for (int k = 0; k < nElNodes; k++){
                    press_[k] = (*nodesCoarse_)[connecC[k]] -> getPressure();
                    velX_[k] = (*nodesCoarse_)[connecC[k]] -> getVelocity(0);
                    velY_[k] = (*nodesCoarse_)[connecC[k]] -> getVelocity(1);
                    velXPrev_[k] = (*nodesCoarse_)[connecC[k]] -> getPreviousVelocity(0);
                    velYPrev_[k] = (*nodesCoarse_)[connecC[k]] -> getPreviousVelocity(1);
                }
                
                matC0.setZero();
                localMV_mat.setZero();
                matE0.setZero();
                
                //Vectors
                vecC0.setZero();
                vecU0.setZero();
                localMV_vec.setZero();
                vecE0.setZero();
                matA0.setZero();
                
                elementsFine_[jel] -> getLagrangeMultipliersDifferentMeshPoisson(iElemCoarse,pspg,press_,velX_,velY_,velXPrev_,velYPrev_,matC0,vecC0,vecU0);

                if (fArlequinStab != ArlequinStabType::ENoStab){
                    elementsFine_[jel] -> getLagrangeMultipliersArlequinDifferentMeshPoisson(iElemCoarse,pspg,press_,velX_,velY_,matE0,matA0,vecE0);
                } 
                
                MatrixDouble matE = matE1;
                stabilizeArlequin(matA0,matA1,matC0,matC1,matE0,vecE0,vecE1,tARLQ0_,tARLQ1_);

                matE0 *= tARLQ0_;
                matA0 *= tARLQ0_;
                vecE0 *= tARLQ0_;
                setMatVecValuesLagMultFineCoarsePoisson(matC0, localMV_mat, matE0, matA0, 
                                                        vecU0, vecC0, localMV_vec, vecE0,
                                                        elementsCoarse_[iElemCoarse] -> getConnectivity(), 
                                                        glueZoneFine_[l] -> getConnectivity());
                
            }; //Number of intersections
        
            matE1 *= tARLQ1_;
            matA1 *= tARLQ1_;
            vecE1 *= tARLQ1_;
            setMatVecValuesLagMultFineFinePoisson(matC1,localMV_mat,matE1,matA1, 
                                                  vecU1,vecC0,localMV_vec,vecE1,
                                                  elementsFine_[jel] -> getConnectivity(),
                                                  glueZoneFine_[l] -> getConnectivity());
        }; // if element belongs to the glue zone
    }; // Glue zone

    return;
}

//------------------------------------------------------------------------------
//----------------COMPUTE ARLEQUIN COUPLED NAVIER-STOKES PROBLEM----------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
int Arlequin<DIM,DEG>::solveArlequinProblem(int iterNumber, double tolerance,
                                            int problem_type, int time_dependency){


    std::ofstream dragLift;
    dragLift.open("dragLift.dat", std::ofstream::out | std::ofstream::app);
    if (rank == 0) {
        dragLift << "Time   Pressure Drag   Pressure Lift " 
                 << "Friction Drag  Friction Lift Drag    Lift " 
                 << std::endl;
    };

    if ((problem_type < 1) || (problem_type > 2)){
        std::cout << "WRONG PROBLEM TYPE." << std::endl;
        return 0;
    };

    if (time_dependency == 0) numTimeSteps = 1;
    
    iTimeStep = 0.;
   
    // //Computes the Weight function for all the finite elements
    setWeightFunction(16.);

    //Computes the Nodal correspondence between fine nodes and coarse elements
    setNodalCorrespondenceFine();
    if (rank == 0) {
        printResultsCoarse(100);
        printResultsFine(100);
    }
    // Computes the system size
    int sysSize;
    if (fProbType == ProblemType::ENavierStokes){
        sysSize = (DIM+1)*numNodesCoarse + (DIM+1)*numNodesFine + DIM*numNodesGlueZoneFine;
    } else if (fProbType == ProblemType::EPoisson){
        sysSize = numNodesCoarse + numNodesFine + numNodesGlueZoneFine;
    }
    double integScheme = fineModel.integScheme;
 
    alpha_f = 1. / (1. + integScheme);
    alpha_m = 0.5 * (3. - integScheme) / (1. + integScheme);
    gamma = 0.5 + alpha_m - alpha_f;

    std::cout << "Time integration parameters = " << alpha_f << ", " << alpha_m << ", " << gamma <<std::endl; 

    for (iTimeStep = 0; iTimeStep < numTimeSteps; iTimeStep++){

        // if(iTimeStep == 10){
        //     double dd = 0.75;
        //     parametersFine -> setSpectralRadius(dd);
        //     parametersCoarse -> setSpectralRadius(dd);
        //     std::cout << "AQUI " << rank << std::endl;
        // }
           
        parametersCoarse -> setTimeInstant(iTimeStep);
        parametersFine -> setTimeInstant(iTimeStep);

        if (rank == 0) {std::cout << "----------------------------" 
                                  << " TIME STEP = "
                                  << iTimeStep 
                                  << " ---------------------------"
                                  << std::endl;}
        PetscMemoryGetCurrentUsage(&bytes);
        PetscPrintf(PETSC_COMM_WORLD,"Memory used22 %g M\n",bytes/(1024*1024));

        //Updates velocity and acceleration
        for (int i = 0; i < numNodesCoarse; i++){
            VecDouble accel(DIM), u(DIM), uprev(DIM);
            
            //Compute acceleration
            u[0] = (*nodesCoarse_)[i] -> getVelocity(0);
            u[1] = (*nodesCoarse_)[i] -> getVelocity(1);

            (*nodesCoarse_)[i] -> setPreviousVelocity(u);
            
            accel[0] = (*nodesCoarse_)[i] -> getAcceleration(0);
            accel[1] = (*nodesCoarse_)[i] -> getAcceleration(1);
            
            (*nodesCoarse_)[i] -> setPreviousAcceleration(accel);

            accel[0] *= (gamma - 1.) / gamma;
            accel[1] *= (gamma - 1.) / gamma;
            
            (*nodesCoarse_)[i] -> setAcceleration(accel);            

        };

        for (int i = 0; i < numNodesFine; i++){
            VecDouble accel(DIM), u(DIM), uprev(DIM), lag(DIM);
            
            //Compute acceleration
            u[0] = (*nodesFine_)[i] -> getVelocity(0);
            u[1] = (*nodesFine_)[i] -> getVelocity(1);

            (*nodesFine_)[i] -> setPreviousVelocity(u);
            
            accel[0] = (*nodesFine_)[i] -> getAcceleration(0);
            accel[1] = (*nodesFine_)[i] -> getAcceleration(1);
            
            (*nodesFine_)[i] -> setPreviousAcceleration(accel);

            accel[0] *= (gamma - 1.) / gamma;
            accel[1] *= (gamma - 1.) / gamma;
            
            (*nodesFine_)[i] -> setAcceleration(accel);
        };

        //STARTS NEWTON-RAPHSON
        for (int inewton = 0; inewton < iterNumber; inewton++){
            
            std::clock_t t1 = std::clock();
            
            // Preallocates the matrix
            if (parametersFine->getSolverType() == SolverType::ESuiteSparse){
                ierr = MatCreateSeqAIJ(PETSC_COMM_WORLD, sysSize, sysSize, 100,NULL,&A);
            } else {
                ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                                sysSize, sysSize,600,NULL,1500,NULL,&A); 
            }

            for (PetscInt i=0; i<sysSize; i++){
                double val = 1.e-10;
                ierr = MatSetValues(A,1,&i,1,&i,&val,ADD_VALUES);   
            }  
            CHKERRQ(ierr);
            
            // Divides the matrix between the processes
            ierr = MatGetOwnershipRange(A, &Istart, &Iend);CHKERRQ(ierr);
            
            //Create PETSc vectors
            ierr = VecCreate(PETSC_COMM_WORLD, &b); CHKERRQ(ierr);
            ierr = VecSetSizes(b, PETSC_DECIDE, sysSize); CHKERRQ(ierr);
            ierr = VecSetFromOptions(b); CHKERRQ(ierr); 
            ierr = VecDuplicate(b, &u); CHKERRQ(ierr);
                        
            // for (int i=0; i<sysSize; i++){
            //     double val = 1.e-20;
            //     ierr = MatSetValues(A,1,&i,1,&i,&val,ADD_VALUES);
                
            // }
            

            std::clock_t t3 = std::clock();
            if (fProbType == ProblemType::ENavierStokes){
                assembleArlequinSystem();
            } else if (fProbType == ProblemType::EPoisson){
                assembleArlequinSystemPoisson();
            }
            
            std::clock_t t4 = std::clock();

            //std::cout << "Enter PETSc " << rank << std::endl;
            
            //Assemble matrices and vectors
            ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
            ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
            
            ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
            ierr = VecAssemblyEnd(b);CHKERRQ(ierr);
            std::clock_t t5 = std::clock();

          
            
 // PetscViewer    viewer;

 // PetscViewerDrawOpen(PETSC_COMM_WORLD,NULL,NULL,0,0,300,300,&viewer);
 // PetscObjectSetName((PetscObject)viewer,"Line graph Plot");
 //  PetscViewerPushFormat(viewer,PETSC_VIEWER_DRAW_LG);

            // ierr = MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
            //ierr = VecView(b,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
            
            //Create KSP context to solve the linear system
            ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
            
            ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);


            switch (parametersFine->getSolverType())
            {
            case SolverType::ESuiteSparse:
                KSPGetPC(ksp, &pc);
                PCSetType(pc, PCLU);
                PCFactorSetMatSolverType(pc, MATSOLVERUMFPACK);
                break;
            case SolverType::EMumps:
                KSPGetPC(ksp, &pc);
                PCSetType(pc, PCLU);
                PCFactorSetMatSolverType(pc, MATSOLVERMUMPS);
                break;

            case SolverType::EIterative:
                KSPSetType(ksp,KSPFGMRES);
                KSPGetPC(ksp, &pc);
                PCSetType(pc,PCBJACOBI);
                KSPSetTolerances(ksp,1.e-10,PETSC_DEFAULT,PETSC_DEFAULT,200);
                break;

            default:
                PanicButton();
                break;
            }


            
            ierr = KSPSolve(ksp,b,u);CHKERRQ(ierr);
            
            ierr = KSPGetTotalIterations(ksp, &iterations);
            
            std::clock_t t6 = std::clock();

            if (rank == 0) std::cout << "TIME " << 1000.*(t4-t3)/CLOCKS_PER_SEC/1000. << " " 
                                                << 1000.*(t5-t4)/CLOCKS_PER_SEC/1000. << " " 
                                                << 1000.*(t6-t5)/CLOCKS_PER_SEC/1000. << std::endl;
            //if (rank == 0)std::cout << "GMRES Iterations = " << iterations << std::endl;
        
            //ierr = VecView(u,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
            
            //Gathers the solution vector to the master process
            ierr = VecScatterCreateToAll(u, &ctx, &All);CHKERRQ(ierr);
            
            ierr = VecScatterBegin(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);
            CHKERRQ(ierr);
            
            ierr = VecScatterEnd(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);
            CHKERRQ(ierr);
            
            ierr = VecScatterDestroy(&ctx);CHKERRQ(ierr);
            
            //Updates nodal values
            double u_[DIM];
            double normU = 0.;
            double normP = 0.;
            double normL = 0.;
            double p_;
            Ione = 1;

            if (fProbType == ProblemType::ENavierStokes){
                for (int i = 0; i < numNodesCoarse; ++i){
                    double w_ = (*nodesCoarse_)[i] -> getWeightFunction();
                    for (int k = 0; k < DIM; k++){
                        Ii = (DIM+1) * i + k;
                        ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
                        // if (nodesCoarse_[i] -> getDistFunction() > -1.2) val *= 1000.e0;
                        u_[k] = val;
                        normU += val*w_*val*w_;
                        (*nodesCoarse_)[i] -> incrementAcceleration(k,u_[k]);
                        (*nodesCoarse_)[i] -> incrementVelocity(k,u_[k]*gamma*dTime);
                    }
                    Ii = (DIM+1) * i + DIM;
                    ierr = VecGetValues(All,Ione,&Ii,&val);CHKERRQ(ierr);
                    p_ = val;
                    normP += val*w_*val*w_;
                    (*nodesCoarse_)[i] -> incrementPressure(p_);
                };
                
                for (int i = 0; i < numNodesFine; ++i){
                    double w_ = (*nodesFine_)[i] -> getWeightFunction();
                    for (int k = 0; k < DIM; k++){
                        Ii = (DIM+1) * numNodesCoarse + (DIM+1) * i + k;
                        ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
                        u_[k] = val;
                        normU += val*w_*val*w_;
                        (*nodesFine_)[i] -> incrementAcceleration(k,u_[k]);
                        (*nodesFine_)[i] -> incrementVelocity(k,u_[k]*gamma*dTime);
                    }        
                    Ii = (DIM+1) * numNodesCoarse + (DIM+1) * i + DIM;
                    ierr = VecGetValues(All,Ione,&Ii,&val);CHKERRQ(ierr);
                    p_ = val;
                    normP += val*w_*val*w_;
                    (*nodesFine_)[i] -> incrementPressure(p_);
                };
                
                for (int i = 0; i < numNodesGlueZoneFine; ++i){
                    for (int k = 0; k < DIM; k++){
                        Ii = (DIM+1) * numNodesCoarse + (DIM+1) * numNodesFine + DIM * i + k;
                        ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
                        u_[k] = val;
                        (*nodesFine_)[nodesGlueZoneFine_[i]] -> incrementLagrangeMultiplier(k,u_[k]);
                        normL += val*val;
                    }
                };
            } else if (fProbType == ProblemType::EPoisson){
                for (int i = 0; i < numNodesCoarse; ++i){
                    double w_ = (*nodesCoarse_)[i] -> getWeightFunction();
                    Ii = i;
                    ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
                    // if (nodesCoarse_[i] -> getDistFunction() > -1.2) val *= 1000.e0;
                    u_[0] = val;
                    normU += val*w_*val*w_;
                    (*nodesCoarse_)[i] -> incrementVelocity(0,u_[0]);
                };
                
                for (int i = 0; i < numNodesFine; ++i){
                    double w_ = (*nodesFine_)[i] -> getWeightFunction();
                    Ii = numNodesCoarse + i;
                    ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
                    u_[0] = val;
                    normU += val*w_*val*w_;
                    (*nodesFine_)[i] -> incrementVelocity(0,u_[0]);
                };
                
                for (int i = 0; i < numNodesGlueZoneFine; ++i){
                    Ii = numNodesCoarse + numNodesFine + i;
                    ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
                    u_[0] = val;
                    (*nodesFine_)[nodesGlueZoneFine_[i]] -> incrementLagrangeMultiplier(0,u_[0]);
                    normL += val*val;

                };
            }
            
            //Computes the solution vector norm
            ierr = VecNorm(u,NORM_2,&val);CHKERRQ(ierr);
            
            std::clock_t t2 = std::clock();
            
            if(rank == 0){
                std::cout<<"Iteration = " << inewton << " (" << iterations <<  
                    ")  Du Norm = " << std::scientific << sqrt(normU) 
                         << " " << sqrt(normP) 
                         << " " << sqrt(normL) 
                         << " " << val << 
                    "  Time (s) = " << std::fixed << 
                    1000.*(t2-t1)/CLOCKS_PER_SEC/1000. << std::endl;
            };
            
            ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
            ierr = VecDestroy(&b); CHKERRQ(ierr);
            ierr = VecDestroy(&u); CHKERRQ(ierr);
            ierr = VecDestroy(&All); CHKERRQ(ierr);
            ierr = MatDestroy(&A); CHKERRQ(ierr);
           // ierr = MatDestroy(&F); CHKERRQ(ierr);
            //ierr = MatDestroy(&C); CHKERRQ(ierr);
            
            if(val <= tolerance){
                break;            
            }; 
        };
        
        if (rank == 0){
            double normUUprev, normU;
            normUUprev = 0.;
            normU = 0.;
            for (int i = 0; i < numNodesCoarse; ++i){
                for (int k = 0; k < DIM; k++){
                    double u, uPr, weight;
                    weight = (*nodesCoarse_)[i] -> getWeightFunction();
                    u = (*nodesCoarse_)[i] -> getVelocity(k) * weight;
                    uPr = (*nodesCoarse_)[i] -> getPreviousVelocity(k) * weight;
                    normUUprev += (u-uPr) * (u-uPr);
                    normU += u*u;
                }
            }
            for (int i = 0; i < numNodesFine; ++i){
                for (int k = 0; k < DIM; k++){
                    double u, uPr, weight;
                    weight = (*nodesFine_)[i] -> getWeightFunction();
                    u = (*nodesFine_)[i] -> getVelocity(k) * weight;
                    uPr = (*nodesFine_)[i] -> getPreviousVelocity(k) * weight;
                    normUUprev += (u-uPr) * (u-uPr);
                    normU += u*u;
                }
            }
            std::cout << "NORM U  " << std::scientific <<  sqrt(normUUprev/normU) << std::endl;
        }

        //Compute real velocity
        ShapeFunction<DIM,DEG>                       shapeQuad;
        VecDouble phi_(nElNodes);
        
        for (int i = 0; i<numNodesFine; i++){
            for (int k = 0; k < DIM; k++) 
                (*nodesFine_)[i] -> setVelocityArlequin(k,(*nodesFine_)[i] -> getVelocity(k));
            (*nodesFine_)[i] -> setPressureArlequin((*nodesFine_)[i] ->getPressure());
        };
        
        for (int i = 0; i<numNodesGlueZoneFine; i++){
            double u_coarse[nElNodes], v_coarse[nElNodes], p_coarse[nElNodes];
            
            double u = 0.;
            double v = 0.;
            double p = 0.;
            
            int elCoarse = (*nodesFine_)[nodesGlueZoneFine_[i]] -> getNodalElemCorrespondence();
            VecDouble xsi = (*nodesFine_)[nodesGlueZoneFine_[i]] -> getNodalXsiCorrespondence();
                        
            VecInt connecCoarse = elementsCoarse_[elCoarse] -> getConnectivity();
            
            for (int j=0; j<nElNodes; j++){
                u_coarse[j] = (*nodesCoarse_)[connecCoarse[j]] -> getVelocity(0);
                v_coarse[j] = (*nodesCoarse_)[connecCoarse[j]] -> getVelocity(1);
                p_coarse[j] = (*nodesCoarse_)[connecCoarse[j]] -> getPressure();
            };
            
            shapeQuad.evaluate(xsi,phi_);
            
            for (int j=0; j<nElNodes; j++){
                u += u_coarse[j] * phi_[j];
                v += v_coarse[j] * phi_[j];
                p += p_coarse[j] * phi_[j];
            };
            
            double wFunc = (*nodesFine_)[nodesGlueZoneFine_[i]] -> 
                getWeightFunction();
            
            double u_int = (*nodesFine_)[nodesGlueZoneFine_[i]] -> getVelocity(0) * wFunc + u * (1. - wFunc);
            double v_int = (*nodesFine_)[nodesGlueZoneFine_[i]] -> getVelocity(1) * wFunc + v * (1. - wFunc);
            double p_int = (*nodesFine_)[nodesGlueZoneFine_[i]] -> getPressure() * wFunc + p * (1. - wFunc);
            
            (*nodesFine_)[nodesGlueZoneFine_[i]] -> setVelocityArlequin(0,u_int);
            (*nodesFine_)[nodesGlueZoneFine_[i]] -> setVelocityArlequin(1,v_int);
            (*nodesFine_)[nodesGlueZoneFine_[i]] -> setPressureArlequin(p_int);
            
        };
        
        for (int i=0; i<numNodesCoarse; i++){
            (*nodesCoarse_)[i] -> setVelocityArlequin(0,(*nodesCoarse_)[i] -> getVelocity(0));
            (*nodesCoarse_)[i] -> setVelocityArlequin(1,(*nodesCoarse_)[i] -> getVelocity(1));
            (*nodesCoarse_)[i] -> setPressureArlequin((*nodesCoarse_)[i] -> getPressure());
        };

        if (parametersCoarse->getExactSolutionPoisson() && parametersFine->getExactSolutionPoisson()) computeErrorPoisson();
        
        // Compute and print drag and lift coefficients
        if (fineModel.getComputeDragAndLift()){
            dragAndLiftCoefficients(dragLift);
        };

        if (rank == 0) {
            //Printing results
            printResultsCoarse(iTimeStep);
            printResultsFine(iTimeStep);
        };        
     };
        
    return 0;
};

//------------------------------------------------------------------------------
//----------------COMPUTE ARLEQUIN COUPLED NAVIER-STOKES PROBLEM----------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
int Arlequin<DIM,DEG>::solveArlequinProblemMoving(int iterNumber, double tolerance,
                                                  int problem_type, 
                                                  int time_dependency){


    std::ofstream dragLift;
    dragLift.open("dragLift.dat", std::ofstream::out | std::ofstream::app);
    if (rank == 0) {
        dragLift << "Time   Pressure Drag   Pressure Lift " 
                 << "Friction Drag  Friction Lift Drag    Lift " 
                 << std::endl;
    };   

    if ((problem_type < 1) || (problem_type > 2)){
        std::cout << "WRONG PROBLEM TYPE." << std::endl;
        return 0;
    };

    if (time_dependency == 0) numTimeSteps = 1;
    
    iTimeStep = 0.;
    
    // //Construct the glue zone based on some defined criterion
    // setCouplingZone();

    // //Computes the Weight function for all the finite elements
    // setWeightFunction(16.);
   
    //Computes the Nodal correspondence between fine nodes and coarse elements
    // setNodalCorrespondenceFine();

    double integScheme = fineModel.integScheme;

  
    // Computes the system size
    int sysSize;
    if (fProbType == ProblemType::ENavierStokes){
        sysSize = (DIM+1)*numNodesCoarse + (DIM+1)*numNodesFine + DIM*numNodesGlueZoneFine;
    } else if (fProbType == ProblemType::EPoisson){
        sysSize = numNodesCoarse + numNodesFine + numNodesGlueZoneFine;
    }

    alpha_f = 1. / (1. + integScheme);
    alpha_m = 0.5 * (3. - integScheme) / (1. + integScheme);
    gamma = 0.5 + alpha_m - alpha_f;

    for (iTimeStep = 0; iTimeStep < numTimeSteps; iTimeStep++){
        
        if (rank == 0) {std::cout << "------------------------- TIME STEP = "
                                  << iTimeStep << " -------------------------"
                                  << std::endl;}
        PetscMemoryGetCurrentUsage(&bytes);
        PetscPrintf(PETSC_COMM_WORLD,"Memory used33 %g M\n",bytes/(1024*1024));
        
        parametersCoarse -> setTimeInstant(iTimeStep);
        parametersFine -> setTimeInstant(iTimeStep);

        //Updates velocity and acceleration
        for (int i = 0; i < numNodesCoarse; i++){
            VecDouble accel(DIM), u(DIM), uprev(DIM);
            
            //Compute acceleration
            u[0] = (*nodesCoarse_)[i] -> getVelocity(0);
            u[1] = (*nodesCoarse_)[i] -> getVelocity(1);

            (*nodesCoarse_)[i] -> setPreviousVelocity(u);
            
            accel[0] = (*nodesCoarse_)[i] -> getAcceleration(0);
            accel[1] = (*nodesCoarse_)[i] -> getAcceleration(1);
            
            (*nodesCoarse_)[i] -> setPreviousAcceleration(accel);

            accel[0] *= (gamma - 1.) / gamma;
            accel[1] *= (gamma - 1.) / gamma;
            
            (*nodesCoarse_)[i] -> setAcceleration(accel);            

        };

        // double f = .35;
        // double w = 2 * pi * f;

        for (int i = 0; i < numNodesFine; i++){
            VecDouble accel(DIM), u(DIM), uprev(DIM), lag(DIM);
            
            //Compute acceleration
            u[0] = (*nodesFine_)[i] -> getVelocity(0);
            u[1] = (*nodesFine_)[i] -> getVelocity(1);

            (*nodesFine_)[i] -> setPreviousVelocity(u);
            
            accel[0] = (*nodesFine_)[i] -> getAcceleration(0);
            accel[1] = (*nodesFine_)[i] -> getAcceleration(1);
            
            (*nodesFine_)[i] -> setPreviousAcceleration(accel);

            accel[0] *= (gamma - 1.) / gamma;
            accel[1] *= (gamma - 1.) / gamma;
            
            (*nodesFine_)[i] -> setAcceleration(accel);



            VecDouble xn(DIM);
            VecDouble xi = (*nodesFine_)[i] -> getInitialCoordinates();       
            VecDouble x = (*nodesFine_)[i] -> getCoordinates();       
    
            double a = -20 * pi / 180 + 10 * pi / 180 * std::cos(2.*pi*iTimeStep*dTime);// + 10 * pi / 180;

            // std::cout << " AAA " << a << std::endl;

            xn[0] = 0.5 + (xi[0]-0.5) * std::cos(a) - (xi[1]-0.0) * std::sin(a);
            xn[1] = 0.0 + (xi[0]-0.5) * std::sin(a) + (xi[1]-0.0) * std::cos(a);

            u[0] = (xn[0] - x[0]) / dTime;
            u[1] = (xn[1] - x[1]) / dTime;


            (*nodesFine_)[i] -> setMeshVelocity(u);
      
            (*nodesFine_)[i] -> setPreviousCoordinates(0,x[0]);
            (*nodesFine_)[i] -> setPreviousCoordinates(1,x[1]);

            (*nodesFine_)[i] -> setCoordinates(xn);
        };
        
        setNodalCorrespondenceFine();
        setSignaledDistance();
        setWeightFunction(1.);
        
        //STARTS NEWTON-RAPHSON
        for (int inewton = 0; inewton < iterNumber; inewton++){
            
            std::clock_t t1 = std::clock();
            
            // Preallocates the matrix
            if (parametersFine->getSolverType() == SolverType::ESuiteSparse){
                ierr = MatCreateSeqAIJ(PETSC_COMM_WORLD, sysSize, sysSize, 100,NULL,&A);
            } else {
                ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                                sysSize, sysSize,400,NULL,600,NULL,&A); 
            }
            for (PetscInt i=0; i<sysSize; i++){
                double val = 1.e-10;
                ierr = MatSetValues(A,1,&i,1,&i,&val,ADD_VALUES);   
            }  
            CHKERRQ(ierr);
            
            // Divides the matrix between the processes
            ierr = MatGetOwnershipRange(A, &Istart, &Iend);CHKERRQ(ierr);
            
            //Create PETSc vectors
            ierr = VecCreate(PETSC_COMM_WORLD, &b); CHKERRQ(ierr);
            ierr = VecSetSizes(b, PETSC_DECIDE, sysSize); CHKERRQ(ierr);
            ierr = VecSetFromOptions(b); CHKERRQ(ierr); 
            ierr = VecDuplicate(b, &u); CHKERRQ(ierr);
                        
            if (fProbType == ProblemType::ENavierStokes){
                assembleArlequinSystem();
            } else if (fProbType == ProblemType::EPoisson){
                assembleArlequinSystemPoisson();
            }
            
            //Assemble matrices and vectors
            ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
            ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
            
            ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
            ierr = VecAssemblyEnd(b);CHKERRQ(ierr);
            
            //ierr = MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
            // ierr = VecView(b,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
            
            //Create KSP context to solve the linear system
            ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
            
            ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
            

            switch (parametersFine->getSolverType())
            {
            case SolverType::ESuiteSparse:
                KSPGetPC(ksp, &pc);
                PCSetType(pc, PCLU);
                PCFactorSetMatSolverType(pc, MATSOLVERUMFPACK);
                break;
            case SolverType::EMumps:
                KSPGetPC(ksp, &pc);
                PCSetType(pc, PCLU);
                PCFactorSetMatSolverType(pc, MATSOLVERMUMPS);
                break;

            case SolverType::EIterative:
                KSPSetType(ksp,KSPFGMRES);
                KSPGetPC(ksp, &pc);
                PCSetType(pc,PCBJACOBI);
                KSPSetTolerances(ksp,1.e-10,PETSC_DEFAULT,PETSC_DEFAULT,200);
                break;

            default:
                PanicButton();
                break;
            }

           // ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);

            ierr = KSPSolve(ksp,b,u);CHKERRQ(ierr);
            
            ierr = KSPGetTotalIterations(ksp, &iterations);
            
            //if (rank == 0)std::cout << "GMRES Iterations = " << iterations << std::endl;
            
            //ierr = VecView(u,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
            
            //Gathers the solution vector to the master process
            ierr = VecScatterCreateToAll(u, &ctx, &All);CHKERRQ(ierr);
            
            ierr = VecScatterBegin(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);
            CHKERRQ(ierr);
            
            ierr = VecScatterEnd(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);
            CHKERRQ(ierr);
            
            ierr = VecScatterDestroy(&ctx);CHKERRQ(ierr);
            
            //Updates nodal values
            double u_[DIM];
            double normU = 0.;
            double normP = 0.;
            double normL = 0.;
            double normT = 0.;
            double p_;
            Ione = 1;

            if (fProbType == ProblemType::ENavierStokes){
                for (int i = 0; i < numNodesCoarse; ++i){
                    double w_ = (*nodesCoarse_)[i] -> getWeightFunction();
                    for (int k = 0; k < DIM; k++){
                        Ii = (DIM+1) * i + k;
                        ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
                        // if (nodesCoarse_[i] -> getDistFunction() > -1.2) val *= 1000.e0;
                        u_[k] = val;
                        normU += val*w_*val*w_;
                        (*nodesCoarse_)[i] -> incrementAcceleration(k,u_[k]);
                        (*nodesCoarse_)[i] -> incrementVelocity(k,u_[k]*gamma*dTime);
                    }
                    Ii = (DIM+1) * i + DIM;
                    ierr = VecGetValues(All,Ione,&Ii,&val);CHKERRQ(ierr);
                    p_ = val;
                    normP += val*w_*val*w_;
                    (*nodesCoarse_)[i] -> incrementPressure(p_);
                };
                
                for (int i = 0; i < numNodesFine; ++i){
                    double w_ = (*nodesFine_)[i] -> getWeightFunction();
                    for (int k = 0; k < DIM; k++){
                        Ii = (DIM+1) * numNodesCoarse + (DIM+1) * i + k;
                        ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
                        u_[k] = val;
                        normU += val*w_*val*w_;
                        (*nodesFine_)[i] -> incrementAcceleration(k,u_[k]);
                        (*nodesFine_)[i] -> incrementVelocity(k,u_[k]*gamma*dTime);
                    }        
                    Ii = (DIM+1) * numNodesCoarse + (DIM+1) * i + DIM;
                    ierr = VecGetValues(All,Ione,&Ii,&val);CHKERRQ(ierr);
                    p_ = val;
                    normP += val*w_*val*w_;
                    (*nodesFine_)[i] -> incrementPressure(p_);
                };
                
                for (int i = 0; i < numNodesGlueZoneFine; ++i){
                    for (int k = 0; k < DIM; k++){
                        Ii = (DIM+1) * numNodesCoarse + (DIM+1) * numNodesFine + DIM * i + k;
                        ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
                        u_[k] = val;
                        (*nodesFine_)[nodesGlueZoneFine_[i]] -> incrementLagrangeMultiplier(k,u_[k]);
                        normL += val*val;
                    }
                };
            } else if (fProbType == ProblemType::EPoisson){
                for (int i = 0; i < numNodesCoarse; ++i){
                    double w_ = (*nodesCoarse_)[i] -> getWeightFunction();
                    Ii = i;
                    ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
                    // if (nodesCoarse_[i] -> getDistFunction() > -1.2) val *= 1000.e0;
                    u_[0] = val;
                    normU += val*w_*val*w_;
                    (*nodesCoarse_)[i] -> incrementVelocity(0,u_[0]);
                };
                
                for (int i = 0; i < numNodesFine; ++i){
                    double w_ = (*nodesFine_)[i] -> getWeightFunction();
                    Ii = numNodesCoarse + i;
                    ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
                    u_[0] = val;
                    normU += val*w_*val*w_;
                    (*nodesFine_)[i] -> incrementVelocity(0,u_[0]);
                };
                
                for (int i = 0; i < numNodesGlueZoneFine; ++i){
                    Ii = numNodesCoarse + numNodesFine + i;
                    ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
                    u_[0] = val;
                    (*nodesFine_)[nodesGlueZoneFine_[i]] -> incrementLagrangeMultiplier(0,u_[0]);
                    normL += val*val;

                };
            }

            //Computes the solution vector norm
            ierr = VecNorm(u,NORM_2,&val);CHKERRQ(ierr);
            
            std::clock_t t2 = std::clock();
            
            if(rank == 0){                
                std::cout<<"Iteration = " << inewton << " (" << iterations <<  
                    ")  Du Norm = " << std::scientific << sqrt(normU) 
                         << " " << sqrt(normP) 
                         << " " << sqrt(normL) 
                         << " " << val << 
                    "  Time (s) = " << std::fixed << 
                    1000.*(t2-t1)/CLOCKS_PER_SEC/1000. << std::endl;
            };
            
            ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
            ierr = VecDestroy(&b); CHKERRQ(ierr);
            ierr = VecDestroy(&u); CHKERRQ(ierr);
            ierr = VecDestroy(&All); CHKERRQ(ierr);
            ierr = MatDestroy(&A); CHKERRQ(ierr);
            
            if(val <= tolerance){
                break;            
            };          

            //Updates SUPG Parameter
            // for (int i = 0; i < numElemFine; i++){
            //     elementsFine_[i] -> getParameterSUPG();
            // };   
            // for (int i = 0; i < numElemCoarse; i++){
            //     elementsCoarse_[i] -> getParameterSUPG();
            // };

        };

        //Compute real velocity
        ShapeFunction<DIM,DEG>                       shapeQuad;
        VecDouble phi_(nElNodes);
        
        for (int i = 0; i<numNodesFine; i++){
            for (int k = 0; k < DIM; k++) 
                (*nodesFine_)[i] -> setVelocityArlequin(k,(*nodesFine_)[i] -> getVelocity(k));
            (*nodesFine_)[i] -> setPressureArlequin((*nodesFine_)[i] ->getPressure());
        };
        
        for (int i = 0; i<numNodesGlueZoneFine; i++){
            double u_coarse[nElNodes], v_coarse[nElNodes], p_coarse[nElNodes];
            
            double u = 0.;
            double v = 0.;
            double p = 0.;
            
            int elCoarse = (*nodesFine_)[nodesGlueZoneFine_[i]] -> getNodalElemCorrespondence();
            VecDouble xsi = (*nodesFine_)[nodesGlueZoneFine_[i]] -> getNodalXsiCorrespondence();
                        
            VecInt connecCoarse = elementsCoarse_[elCoarse] -> getConnectivity();
            
            for (int j=0; j<nElNodes; j++){
                u_coarse[j] = (*nodesCoarse_)[connecCoarse[j]] -> getVelocity(0);
                v_coarse[j] = (*nodesCoarse_)[connecCoarse[j]] -> getVelocity(1);
                p_coarse[j] = (*nodesCoarse_)[connecCoarse[j]] -> getPressure();
            };
            
            shapeQuad.evaluate(xsi,phi_);
            
            for (int j=0; j<nElNodes; j++){
                u += u_coarse[j] * phi_[j];
                v += v_coarse[j] * phi_[j];
                p += p_coarse[j] * phi_[j];
            };
            
            double wFunc = (*nodesFine_)[nodesGlueZoneFine_[i]] -> 
                getWeightFunction();
            
            double u_int = (*nodesFine_)[nodesGlueZoneFine_[i]] -> getVelocity(0) * wFunc + u * (1. - wFunc);
            double v_int = (*nodesFine_)[nodesGlueZoneFine_[i]] -> getVelocity(1) * wFunc + v * (1. - wFunc);
            double p_int = (*nodesFine_)[nodesGlueZoneFine_[i]] -> getPressure() * wFunc + p * (1. - wFunc);
            
            (*nodesFine_)[nodesGlueZoneFine_[i]] -> setVelocityArlequin(0,u_int);
            (*nodesFine_)[nodesGlueZoneFine_[i]] -> setVelocityArlequin(1,v_int);
            (*nodesFine_)[nodesGlueZoneFine_[i]] -> setPressureArlequin(p_int);
            
        };
        
        for (int i=0; i<numNodesCoarse; i++){
            (*nodesCoarse_)[i] -> setVelocityArlequin(0,(*nodesCoarse_)[i] -> getVelocity(0));
            (*nodesCoarse_)[i] -> setVelocityArlequin(1,(*nodesCoarse_)[i] -> getVelocity(1));
            (*nodesCoarse_)[i] -> setPressureArlequin((*nodesCoarse_)[i] -> getPressure());
        };
        
        // Compute and print drag and lift coefficients
        if (fineModel.getComputeDragAndLift()){
            dragAndLiftCoefficients(dragLift);
        };

        if (rank == 0) {
                       
            //Printing results
            printResultsCoarse(iTimeStep);
            printResultsFine(iTimeStep);
        };
    };
        
    return 0;

};



//------------------------------------------------------------------------------
//----------------COMPUTE ARLEQUIN COUPLED NAVIER-STOKES PROBLEM----------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
int Arlequin<DIM,DEG>::solveFSIArlequin(int iterNumber, double tolerance,
                                        int problem_type, int iTimeStep){


    std::ofstream dragLift;
    dragLift.open("dragLift.dat", std::ofstream::out | std::ofstream::app);

    if ((problem_type < 1) || (problem_type > 2)){
        std::cout << "WRONG PROBLEM TYPE." << std::endl;
        return 0;
    };
    
    // Computes the system size
    int sysSize = 3 * numNodesCoarse + 3 * numNodesFine + 2 * numNodesGlueZoneFine;
    
    setNodalCorrespondenceFine();
    setSignaledDistance();
    setWeightFunction(1.);

    double &alpha_f = parametersFine -> getAlphaF();
    double &alpha_m = parametersFine -> getAlphaM();
    double &gamma = parametersFine -> getGamma();
    
    //STARTS NEWTON-RAPHSON
    for (int inewton = 0; inewton < iterNumber; inewton++){
        
        std::clock_t t1 = std::clock();
        
        // Preallocates the matrix
        if (parametersFine->getSolverType() == SolverType::ESuiteSparse){
            ierr = MatCreateSeqAIJ(PETSC_COMM_WORLD, sysSize, sysSize, 100,NULL,&A);
        } else {
            ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                            sysSize, sysSize,600,NULL,1500,NULL,&A); 
        }

        for (PetscInt i=0; i<sysSize; i++){
            double val = 1.e-10;
            ierr = MatSetValues(A,1,&i,1,&i,&val,ADD_VALUES);   
        }  
        
        CHKERRQ(ierr);
        
        // Divides the matrix between the processes
        ierr = MatGetOwnershipRange(A, &Istart, &Iend);CHKERRQ(ierr);
        
        //Create PETSc vectors
        ierr = VecCreate(PETSC_COMM_WORLD, &b); CHKERRQ(ierr);
        ierr = VecSetSizes(b, PETSC_DECIDE, sysSize); CHKERRQ(ierr);
        ierr = VecSetFromOptions(b); CHKERRQ(ierr); 
        ierr = VecDuplicate(b, &u); CHKERRQ(ierr);
                            
        assembleArlequinSystem();
        
        //Assemble matrices and vectors
        ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        
        ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
        ierr = VecAssemblyEnd(b);CHKERRQ(ierr);
        
        //Create KSP context to solve the linear system
        ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
        
        ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
        
        // // ierr = MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE,0, NULL, &nullsp);
        // // ierr = MatSetNullSpace(A, nullsp);
        // // ierr = MatNullSpaceDestroy(&nullsp);


        // // if (iTimeStep > 5){

        //     ierr = KSPSetTolerances(ksp,1.e-7,1.e-10,PETSC_DEFAULT,
        //                             15);CHKERRQ(ierr);
            
        //     //ierr = KSPGMRESSetRestart(ksp, 10); CHKERRQ(ierr);
            
        //     ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
            
        //     ierr = PCSetType(pc,PCJACOBI);CHKERRQ(ierr);
            
        //     //ierr = KSPSetPCSide(ksp, PC_RIGHT);
        //     ierr = KSPSetType(ksp,KSPGMRES); CHKERRQ(ierr);
            
        //     ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
        //     // ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
        // // }else{
       

        switch (parametersFine->getSolverType())
        {
        case SolverType::ESuiteSparse:
            KSPGetPC(ksp, &pc);
            PCSetType(pc, PCLU);
            PCFactorSetMatSolverType(pc, MATSOLVERUMFPACK);
            break;
        case SolverType::EMumps:
            KSPGetPC(ksp, &pc);
            PCSetType(pc, PCLU);
            PCFactorSetMatSolverType(pc, MATSOLVERMUMPS);
            break;

        case SolverType::EIterative:
            KSPSetType(ksp,KSPFGMRES);
            KSPGetPC(ksp, &pc);
            PCSetType(pc,PCBJACOBI);
            KSPSetTolerances(ksp,1.e-10,PETSC_DEFAULT,PETSC_DEFAULT,200);
            break;

        default:
            PanicButton();
            break;
        }

        ierr = KSPSolve(ksp,b,u);CHKERRQ(ierr);
        
        ierr = KSPGetTotalIterations(ksp, &iterations);
        
        //if (rank == 0)std::cout << "GMRES Iterations = " << iterations << std::endl;
    
        //ierr = VecView(u,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
        
        //Gathers the solution vector to the master process
        ierr = VecScatterCreateToAll(u, &ctx, &All);CHKERRQ(ierr);
        
        ierr = VecScatterBegin(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);
        CHKERRQ(ierr);
        
        ierr = VecScatterEnd(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);
        CHKERRQ(ierr);
        
        ierr = VecScatterDestroy(&ctx);CHKERRQ(ierr);
        
        //Updates nodal values
        double u_[DIM];
        double normU = 0.;
        double normP = 0.;
        double normL = 0.;
        double normT = 0.;
        double p_;
        Ione = 1;

        for (int i = 0; i < numNodesCoarse; ++i){
            for (int k = 0; k < DIM; k++){
                Ii = (DIM+1) * i + k;
                ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
                // if (nodesCoarse_[i] -> getDistFunction() > -1.2) val *= 1000.e0;
                u_[k] = val;
                normU += val*val;
                (*nodesCoarse_)[i] -> incrementAcceleration(k,u_[k]);
                (*nodesCoarse_)[i] -> incrementVelocity(k,u_[k]*gamma*dTime);
            }
            Ii = (DIM+1) * i + DIM;
            ierr = VecGetValues(All,Ione,&Ii,&val);CHKERRQ(ierr);
            p_ = val;
            normP += val*val;
            (*nodesCoarse_)[i] -> incrementPressure(p_);
        };
        for (int i = 0; i < numNodesFine; ++i){
            for (int k = 0; k < DIM; k++){
                Ii = (DIM+1) * numNodesCoarse + (DIM+1) * i + k;
                ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
                u_[k] = val;
                normU += val*val;
                (*nodesFine_)[i] -> incrementAcceleration(k,u_[k]);
                (*nodesFine_)[i] -> incrementVelocity(k,u_[k]*gamma*dTime);
            }        
            Ii = (DIM+1) * numNodesCoarse + (DIM+1) * i + DIM;
            ierr = VecGetValues(All,Ione,&Ii,&val);CHKERRQ(ierr);
            p_ = val;
            normP += val*val;
            (*nodesFine_)[i] -> incrementPressure(p_);
        };
        for (int i = 0; i < numNodesGlueZoneFine; ++i){
            for (int k = 0; k < DIM; k++){
                Ii = (DIM+1) * numNodesCoarse + (DIM+1) * numNodesFine + DIM * i + k;
                ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
                u_[k] = val;
                (*nodesFine_)[nodesGlueZoneFine_[i]] -> incrementLagrangeMultiplier(k,u_[k]);
                normL += val*val;
            }
        };
        
        //Computes the solution vector norm
        ierr = VecNorm(u,NORM_2,&val);CHKERRQ(ierr);
        
        std::clock_t t2 = std::clock();
        
        if(rank == 0){
            std::cout<<"Iteration = " << inewton << " (" << iterations <<  
                ")  Du Norm = " << std::scientific << sqrt(normU) 
                     << " " << sqrt(normP) 
                     << " " << sqrt(normL) 
                     << " " << val << 
                "  Time (s) = " << std::fixed << 
                1000.*(t2-t1)/CLOCKS_PER_SEC/1000. << std::endl;
        };
        
        ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
        ierr = VecDestroy(&b); CHKERRQ(ierr);
        ierr = VecDestroy(&u); CHKERRQ(ierr);
        ierr = VecDestroy(&All); CHKERRQ(ierr);
        ierr = MatDestroy(&A); CHKERRQ(ierr);
        // ierr = MatDestroy(&F); CHKERRQ(ierr);
        
        
        if(val <= tolerance){
            break;            
        };          
        
    };

    if (rank == 0){
        double normUUprev, normU;
        normUUprev = 0.;
        normU = 0.;
        for (int i = 0; i < numNodesCoarse; ++i){
            for (int k = 0; k < DIM; k++){
                double u, uPr, weight;
                weight = (*nodesCoarse_)[i] -> getWeightFunction();
                u = (*nodesCoarse_)[i] -> getVelocity(k) * weight;
                uPr = (*nodesCoarse_)[i] -> getPreviousVelocity(k) * weight;
                normUUprev += (u-uPr) * (u-uPr);
                normU += u*u;
            }
        }
        for (int i = 0; i < numNodesFine; ++i){
            for (int k = 0; k < DIM; k++){
                double u, uPr, weight;
                weight = (*nodesFine_)[i] -> getWeightFunction();
                u = (*nodesFine_)[i] -> getVelocity(k) * weight;
                uPr = (*nodesFine_)[i] -> getPreviousVelocity(k) * weight;
                normUUprev += (u-uPr) * (u-uPr);
                normU += u*u;
            }
        }
        std::cout << "NORM U  " << std::scientific <<  sqrt(normUUprev / normU) << std::endl;
    }


    // std::cout << "AQUI1 " << rank << std::endl;


    //Compute real velocity
    ShapeFunction<DIM,DEG>                       shapeQuad;
    VecDouble phi_(nElNodes);
    
    for (int i = 0; i<numNodesFine; i++){
        (*nodesFine_)[i] -> setVelocityArlequin(0,(*nodesFine_)[i] -> getVelocity(0));
        (*nodesFine_)[i] -> setVelocityArlequin(1,(*nodesFine_)[i] -> getVelocity(1));
        (*nodesFine_)[i] -> setPressureArlequin((*nodesFine_)[i] ->getPressure());
    };
    
    for (int i = 0; i<numNodesGlueZoneFine; i++){
        double u_coarse[nElNodes], v_coarse[nElNodes], p_coarse[nElNodes];
        
        double u = 0.;
        double v = 0.;
        double p = 0.;
        
        int elCoarse = (*nodesFine_)[nodesGlueZoneFine_[i]] -> getNodalElemCorrespondence();
        VecDouble xsi = (*nodesFine_)[nodesGlueZoneFine_[i]] -> getNodalXsiCorrespondence();
        
        VecInt connecCoarse = elementsCoarse_[elCoarse] -> getConnectivity();
        
        for (int j=0; j<nElNodes; j++){
            u_coarse[j] = (*nodesCoarse_)[connecCoarse[j]] -> getVelocity(0);
            v_coarse[j] = (*nodesCoarse_)[connecCoarse[j]] -> getVelocity(1);
            p_coarse[j] = (*nodesCoarse_)[connecCoarse[j]] -> getPressure();
        };
        
        shapeQuad.evaluate(xsi,phi_);
        
        for (int j=0; j<nElNodes; j++){
            u += u_coarse[j] * phi_[j];
            v += v_coarse[j] * phi_[j];
            p += p_coarse[j] * phi_[j];
        };
        
        double wFunc = (*nodesFine_)[nodesGlueZoneFine_[i]] -> 
            getWeightFunction();
        
        double u_int = (*nodesFine_)[nodesGlueZoneFine_[i]] -> getVelocity(0) * wFunc + u * (1. - wFunc);
        double v_int = (*nodesFine_)[nodesGlueZoneFine_[i]] -> getVelocity(1) * wFunc + v * (1. - wFunc);
        double p_int = (*nodesFine_)[nodesGlueZoneFine_[i]] -> getPressure() * wFunc + p * (1. - wFunc);
        
        (*nodesFine_)[nodesGlueZoneFine_[i]] -> setVelocityArlequin(0,u_int);
        (*nodesFine_)[nodesGlueZoneFine_[i]] -> setVelocityArlequin(1,v_int);
        (*nodesFine_)[nodesGlueZoneFine_[i]] -> setPressureArlequin(p_int);
        
    };
    // std::cout << "AQUI2 " << rank << std::endl;

    for (int i=0; i<numNodesCoarse; i++){
        (*nodesCoarse_)[i] -> setVelocityArlequin(0,(*nodesCoarse_)[i] -> getVelocity(0));
        (*nodesCoarse_)[i] -> setVelocityArlequin(1,(*nodesCoarse_)[i] -> getVelocity(1));
        (*nodesCoarse_)[i] -> setPressureArlequin((*nodesCoarse_)[i] -> getPressure());
    };
    // std::cout << "AQUI3 " << rank << std::endl;
    
    return 0;

};

template<int DIM, int DEG>
void Arlequin<DIM,DEG>::computeErrorPoisson() {

    VecDouble errorFine(3),errorCoarse(3),errorTotal(3);

    coarseModel.computeError(errorCoarse);
    fineModel.computeError(errorFine);

    std::ofstream rprint("errorsArlequin.txt",std::ios::app);
    errorTotal = errorFine + errorCoarse;

    if (rank == 0){
        std::cout << "\n\nERROR REPORT:\n" << std::scientific << std::setprecision(10)
            << "L2 state var = " << sqrt(errorTotal[0]) << "\n" 
            << "Semi H1 state var = " << sqrt(errorTotal[1]) << "\n" 
            << "H1 state var = " << sqrt(errorTotal[2]) << "\n"; 
        rprint << sqrt(errorTotal[0]) << " " << sqrt(errorTotal[1]) << " " << sqrt(errorTotal[2]) << std::endl;
    }


}

template<int DIM, int DEG>
void Arlequin<DIM,DEG>::stabilizeArlequin(MatrixDouble &A0, MatrixDouble &A1, 
                                          MatrixDouble &C0, MatrixDouble &C1,
                                          MatrixDouble &E, VecDouble &b0, 
                                          VecDouble &b1, double &tArlq0, double &tArlq1){
    
    //There are in general three main options for taking the norm of a matrix: the L2, L2 and Linfty norms.
    //In eigen they can be simply obtained by:
    // normL2 = mat.norm(); 
    // normL1 = mat.lpNorm<1>(); 
    // normInfty = mat.lpNorm<Infinity>();

    switch (fArlequinStab)
    {
    case ArlequinStabType::ENoStab:
        {
            tArlq0 = 0.;
            tArlq1 = 0.;
            break;
        }
    
    case ArlequinStabType::EOption1:
        {
            double normC0 = C0.norm();
            double normC1 = C1.norm();
            double normA0 = A0.norm();
            double normA1 = A1.norm();
            double normB0 = b0.norm();
            double normB1 = b1.norm();
            double normE = E.norm();
            tArlq0 = std::min({normC0/normA0, normC0/normE, normC0/normB0});
            tArlq1 = std::min({normC1/normA1, normC1/normE, normC1/normB1});
            break;
        }

    case ArlequinStabType::EOption2:
        {
            double normC0 = C0.norm();
            double normC1 = C1.norm();
            double normA0 = A0.norm();
            double normA1 = A1.norm();
            double normB0 = b0.norm();
            double normB1 = b1.norm();
            double normE = E.norm();
            double aux1 = std::min({normC0/normA0, normC0/normE, normC0/normB0, normC1/normA1, normC1/normE, normC1/normB1});
            tArlq0 = aux1;
            tArlq1 = aux1;
            break;
        }
    
    case ArlequinStabType::EOption3:
        {
            double normC0 = C0.norm();
            double normC1 = C1.norm();
            double normA0 = A0.norm();
            double normA1 = A1.norm();
            double normE = E.norm();
            double normB0 = b0.norm();
            double normB1 = b1.norm();
            
            tArlq0 = std::min({normC0/normA0, normC0/normE, normC0/normB0});
            tArlq1 = std::min({normC1/normA1, normC1/normE, normC1/normB1});

            double normC = std::min(normC0,normC1);
            tArlq0 *= normC/normC0;
            tArlq1 *= normC/normC1;
            break;
        }

    case ArlequinStabType::EOption4:
        {
            double normC1 = C1.norm();
            double normA1 = A1.norm();
            double normE = E.norm();
            double normB1 = b1.norm();
            tArlq0 = 0.;
            tArlq1 = std::min({normC1/normA1, normC1/normE, normC1/normB1});
            break;
        }

    case ArlequinStabType::EOption5:
        {
            double normC0 = C0.norm();
            double normC1 = C1.norm();
            double normA1 = A1.norm();
            double normE = E.norm();
            double normB1 = b1.norm();
            double normC = std::min(normC0,normC1);
            tArlq0 = 0.;
            tArlq1 = std::min({normC/normA1, normC/normE, normC/normB1});
            break;
        }
    default:
        PanicButton();
        break;
    }

}

template<int DIM, int DEG>
void Arlequin<DIM,DEG>::stabilizeArlequinNew(VecDouble &Ml1, VecDouble &t1, 
                             VecDouble &j1, VecDouble &k1, VecDouble &p1, double &tArleq){
    
    double mnorm = Ml1.norm();
    double tnorm = t1.norm();
    double jnorm = j1.norm();
    double knorm = k1.norm();
    double pnorm = p1.norm();

    if (fabs(mnorm) < 1.e-20) mnorm = 1.e3;
    if (fabs(tnorm) < 1.e-20) tnorm = 1.e3;
    if (fabs(jnorm) < 1.e-20) jnorm = 1.e3;
    if (fabs(knorm) < 1.e-20) knorm = 1.e3;
    if (fabs(pnorm) < 1.e-20) pnorm = 1.e3;

    double tauA = mnorm/tnorm;
    double tauB = mnorm/jnorm;
    double tauC = mnorm/knorm;
    double tauD = mnorm/pnorm;

    

    tArleq = 1.e-4*pow((1/(tauA*tauA) + 1/(tauB*tauB) + 1/(tauC*tauC) + 1/(tauD*tauD)), -0.5);
    // std::cout << "Ml1 " << Ml1 << std::endl;
    // std::cout << "t1 " << t1 << std::endl;
    // std::cout << "j1 " << j1 << std::endl;
    // std::cout << "k1 " << k1 << std::endl;
    // std::cout << "p1 " << p1 << std::endl;
    // std::cout << "tauA = " << tauA << " , tauB = " << tauB << " , tauC = " << tauC << " , tauD = " << tauD << std::endl;
    // std::cout << "tArleq = " << tArleq << std::endl;
};

template<int DIM, int DEG>
PetscErrorCode Arlequin<DIM,DEG>::FormJacobian(SNES snes,Vec vecU, Mat matA, Mat matB,void *ptr){

    return 0;
}

template<int DIM, int DEG>
PetscErrorCode Arlequin<DIM,DEG>::FormFunction(SNES snes, Vec vecU,Vec vecB, void *ptr){

    return 0;
}



template class Arlequin<2,1>;
template class Arlequin<2,2>;
template class Arlequin<2,3>;
template class Arlequin<3,1>;
template class Arlequin<3,2>;



