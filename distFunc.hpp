//------------------------------------------------------------------------------
//-------------------------DISTANCE FUNCTION FINE MESH--------------------------
//------------------------------------------------------------------------------
// defines elements in the gluing zone (where to enhance quadrature)
// computes signed distance function to the fine mesh boundary
template<>
void Glue<2>::distFunctionFineMesh(){
    //TO BE CALLED AFTER DISTFUNCTION COARSE MESH
    typename Elements::Connectivity connec;
    typename Nodes::VecLocD x, x1, x2, n, test;
    typename Boundary::BoundConnect bconnec;
    double dist;
       
    
    for (int ino = 0; ino < numNodesFine; ino++){
        x = nodesFine_[ino] -> getCoordinates();
        dist=10000000000000000000000000000.;
        
        for (int i = 0; i < numBoundariesFine; i++){
            if (boundaryFine_[i]->getConstrain(0)==2){
                bconnec = boundaryFine_[i]->getBoundaryConnectivity();
                //first segment
                int no1 = bconnec(0);
                int no2 = bconnec(2);
                // std::cout<<no1<<" nos "<<no2<<std::endl;
              
                x1 = nodesFine_[no1] -> getCoordinates();
                x2 = nodesFine_[no2] -> getCoordinates();
                
                double aux0 =  std::sqrt((x2(1)-x1(1))*(x2(1)-x1(1))+(x2(0)-x1(0))*(x2(0)-x1(0)));
                double aux1 = ((x(0)-x1(0))*(x2(0)-x1(0))+(x(1)-x1(1))*(x2(1)-x1(1)))/aux0;
                double dist2 =-((x2(1)-x1(1))*x(0)-(x2(0)-x1(0))*x(1)+x2(0)*x1(1)-x2(1)*x1(0))/aux0;
                
                if(aux1>aux0){
                    dist2 =  std::sqrt((x2(1)-x(1))*(x2(1)-x(1))+(x2(0)-x(0))*(x2(0)-x(0)));
                    //find signal
                    //side normal vector
                    n = nodesFine_[no2] -> getInnerNormal();
                    
                    test(0)=x(0)-x2(0);
                    test(1)=x(1)-x2(1);
                    double signaltest = inner_prod(n,test);
                    double signal=-1.;
                    
                    if (signaltest < 0.)signal = 1.;
                    
                    dist2 *= signal;
                    
                    
                };
                if(aux1<0.){
                    dist2 =std::sqrt((x(1)-x1(1))*(x(1)-x1(1))+(x(0)-x1(0))*(x(0)-x1(0)));
                    //find signal
                    //side normal vector
                    n = nodesFine_[no1] -> getInnerNormal();
                    
                    test(0)=x(0)-x1(0);
                    test(1)=x(1)-x1(1);
                    double signaltest = inner_prod(n,test);
                    double signal=-1.;
                    
                    if (signaltest < 0.)signal = 1.;
                    
                    dist2 *= signal;
                };
                
                if (fabs(dist2) < fabs(dist)) dist = dist2;
                
                //second segment
                no1 = bconnec(2);
                no2 = bconnec(1);
                x1 = nodesFine_[no1] -> getCoordinates();
                x2 = nodesFine_[no2] -> getCoordinates();
                
                aux0 =  std::sqrt((x2(1)-x1(1))*(x2(1)-x1(1))+(x2(0)-x1(0))*(x2(0)-x1(0)));
                aux1 = ((x(0)-x1(0))*(x2(0)-x1(0))+(x(1)-x1(1))*(x2(1)-x1(1)))/aux0;
                dist2 =-((x2(1)-x1(1))*x(0)-(x2(0)-x1(0))*x(1)+x2(0)*x1(1)-x2(1)*x1(0))/aux0;
                if(aux1>aux0){
                    dist2 =  std::sqrt((x2(1)-x(1))*(x2(1)-x(1))+(x2(0)-x(0))*(x2(0)-x(0)));
                    n = nodesFine_[no2] -> getInnerNormal();
                    
                    test(0)=x(0)-x2(0);
                    test(1)=x(1)-x2(1);
                    double signaltest = inner_prod(n,test);
                    double signal=-1.;
                    
                    if (signaltest < 0.)signal = 1.;
                    
                    dist2 *= signal;
                    
                    
                };
                if(aux1<0.){
                    dist2 =std::sqrt((x(1)-x1(1))*(x(1)-x1(1))+(x(0)-x1(0))*(x(0)-x1(0)));
        //find signal
                    //side normal vector
                    n = nodesFine_[no1] -> getInnerNormal();
                    
                    test(0)=x(0)-x1(0);
                    test(1)=x(1)-x1(1);
                    double signaltest = inner_prod(n,test);
                    double signal=-1.;
                    
                    if (signaltest < 0.)signal = 1.;
                    
                    dist2 *= signal;
                    
                };
                if (fabs(dist2) < fabs(dist)) dist = dist2;
            }; //if bf is the glue boundary
        }; //
    
        nodesFine_[ino] -> setDistFunction(dist);
    };
};

//------------------------------------------------------------------------------
//-----------------------DISTANCE FUNCTION COARSE MESH--------------------------
//------------------------------------------------------------------------------
// computes coarse mesh nodes distance function to the fine mesh boundary
template<>
void Glue<2>::distFunctionCoarseMesh(){

    typename Elements::Connectivity connec;
    typename Nodes::VecLocD x, x1, x2, n, test;
    typename Boundary::BoundConnect bconnec;
    double dist;
    
    for (int inode = 0 ; inode < numNodesFine; inode++){ 
        nodesFine_[inode] -> clearInnerNormal();
    };

    
    //approximate normal calculation
    for (int i = 0; i < numBoundariesFine; i++){
        if (boundaryFine_[i]->getConstrain(0)==2){
            bconnec = boundaryFine_[i]->getBoundaryConnectivity();
            //first segment
            int no1 = bconnec(0);
            int no2 = bconnec(2);
            x1 = nodesFine_[no1] -> getCoordinates();
            x2 = nodesFine_[no2] -> getCoordinates();

            double sLength = sqrt((x2(1)-x1(1))*(x2(1)-x1(1))+(x1(0)-x2(0))*(x1(0)-x2(0)));
            n(0)=(x2(1)-x1(1))/sLength;
            n(1)=(x1(0)-x2(0))/sLength;
            nodesFine_[no1] -> setInnerNormal(n);
            nodesFine_[no2] -> setInnerNormal(n);
            //second segment
            no1 = bconnec(2);
            no2 = bconnec(1);
            x1 = nodesFine_[no1] -> getCoordinates();
            x2 = nodesFine_[no2] -> getCoordinates();

            sLength = sqrt((x2(1)-x1(1))*(x2(1)-x1(1))+(x1(0)-x2(0))*(x1(0)-x2(0)));

            n(0)=(x2(1)-x1(1))/sLength;
            n(1)=(x1(0)-x2(0))/sLength;
            nodesFine_[no1] -> setInnerNormal(n);
            nodesFine_[no2] -> setInnerNormal(n);
            
        };
    };
    // std::cout<<"zerar normais apos mesh update"<<std::endl;
    
    for (int ino = 0; ino < numNodesCoarse; ino++){
        x = nodesCoarse_[ino] -> getCoordinates();
        dist=10000000000000000000000000000.;
        
        for (int i = 0; i < numBoundariesFine; i++){
            if (boundaryFine_[i]->getConstrain(0)==2){
                bconnec = boundaryFine_[i]->getBoundaryConnectivity();
                //first segment
                int no1 = bconnec(0);
                int no2 = bconnec(2);
                x1 = nodesFine_[no1] -> getCoordinates();
                x2 = nodesFine_[no2] -> getCoordinates();
                
                double aux0 =  std::sqrt((x2(1)-x1(1))*(x2(1)-x1(1))+(x2(0)-x1(0))*(x2(0)-x1(0)));
                double aux1 = ((x(0)-x1(0))*(x2(0)-x1(0))+(x(1)-x1(1))*(x2(1)-x1(1)))/aux0;
                double dist2 =-((x2(1)-x1(1))*x(0)-(x2(0)-x1(0))*x(1)+x2(0)*x1(1)-x2(1)*x1(0))/aux0;
                
                if(aux1>aux0){
                    dist2 =  std::sqrt((x2(1)-x(1))*(x2(1)-x(1))+(x2(0)-x(0))*(x2(0)-x(0)));
                    //find signal
                    //side normal vector
                    n = nodesFine_[no2] -> getInnerNormal();
                    
                    test(0)=x(0)-x2(0);
                    test(1)=x(1)-x2(1);
                    double signaltest = inner_prod(n,test);
                    double signal=-1.;
                    
                    if (signaltest < 0.)signal = 1.;
                    
                    dist2 *= signal;
                    
                    
                };
                if(aux1<0.){
                    dist2 =std::sqrt((x(1)-x1(1))*(x(1)-x1(1))+(x(0)-x1(0))*(x(0)-x1(0)));
                    //find signal
                    //side normal vector
                    n = nodesFine_[no1] -> getInnerNormal();
                    
                    test(0)=x(0)-x1(0);
                    test(1)=x(1)-x1(1);
                    double signaltest = inner_prod(n,test);
                    double signal=-1.;
                    
                    if (signaltest < 0.)signal = 1.;
                    
                    dist2 *= signal;
                };
                
                if (fabs(dist2) < fabs(dist)) dist = dist2;
                
                //second segment
                no1 = bconnec(2);
                no2 = bconnec(1);
                x1 = nodesFine_[no1] -> getCoordinates();
                x2 = nodesFine_[no2] -> getCoordinates();
                
                aux0 =  std::sqrt((x2(1)-x1(1))*(x2(1)-x1(1))+(x2(0)-x1(0))*(x2(0)-x1(0)));
                aux1 = ((x(0)-x1(0))*(x2(0)-x1(0))+(x(1)-x1(1))*(x2(1)-x1(1)))/aux0;
                dist2 =-((x2(1)-x1(1))*x(0)-(x2(0)-x1(0))*x(1)+x2(0)*x1(1)-x2(1)*x1(0))/aux0;
                if(aux1>aux0){
                    dist2 =  std::sqrt((x2(1)-x(1))*(x2(1)-x(1))+(x2(0)-x(0))*(x2(0)-x(0)));
                    n = nodesFine_[no2] -> getInnerNormal();
                    
                    test(0)=x(0)-x2(0);
                    test(1)=x(1)-x2(1);
                    double signaltest = inner_prod(n,test);
                    double signal=-1.;
                    
                    if (signaltest < 0.)signal = 1.;
                    
                    dist2 *= signal;
                    
                    
                };
                if(aux1<0.){
                    dist2 =std::sqrt((x(1)-x1(1))*(x(1)-x1(1))+(x(0)-x1(0))*(x(0)-x1(0)));
        //find signal
                    //side normal vector
                    n = nodesFine_[no1] -> getInnerNormal();
                    
                    test(0)=x(0)-x1(0);
                    test(1)=x(1)-x1(1);
                    double signaltest = inner_prod(n,test);
                    double signal=-1.;
                    
                    if (signaltest < 0.)signal = 1.;
                    
                    dist2 *= signal;
                    
                };
                if (fabs(dist2) < fabs(dist)) dist = dist2;
            }; //if bf is the glue boundary
        }; //
    
        nodesCoarse_[ino] -> setDistFunction(dist);
    };

    
};
