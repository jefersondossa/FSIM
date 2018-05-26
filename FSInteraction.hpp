//------------------------------------------------------------------------------
// 
//                   Jeferson W D Fernandes and Rodolfo A K Sanches
//                             University of Sao Paulo
//                           (C) 2017 All Rights Reserved
//
// <LicenseText>
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//------------------------------------SOLID-------------------------------------
//------------------------------------------------------------------------------

#ifndef FSINTERACTION_H
#define FSINTERACTION_H

#include <boost/timer.hpp> 
#include <boost/thread.hpp>

#include "Arlequin_cyl.hpp"

//Solid extern functions (from porticomb.for)
extern "C" {void preprocessing_(char *solid_reading);};
extern "C" {void solveframestructure_(int *ipt);};
extern "C" {void getnumberofnodessolid_(int *n);};
extern "C" {void getnumberofelementssolid_(int *n);};
extern "C" {void getsolidposition_(int *node, double *x_, double *y_);};
extern "C" {void searchcorrespondencefluid_(double *x_, double *y_, 
                                            int *elem, double *xsi);};
extern "C" {void getupdatedcoordinates_(double *x_, double *y_, 
                                        int *elem, double *xsi);};
extern "C" {void getinterpolatedvelocity_(double *x_, double *y_, 
                                          int *elem, double *xsi);};
extern "C" {void setcouplingload_(double *fx_, double *fy_, int *node);};
extern "C" {void clearcouplingloads_();};
extern "C" {void printstructure_();};
extern "C" {void updateqsrs_();};
extern "C" {void getposition_(int *dof, double *p);};
extern "C" {void setposition_(int *dof, double *p);};
extern "C" {void getvelocity_(int *dof, double *p);};
extern "C" {void getpreviousvelocity_(int *dof, double *p);};
extern "C" {void updatesolid_(int *ipt);};

/// Mounts and solve the Fluid-Structure interaction problem.

template<int DIM>
class FSInteraction{
public:
    /// Defines locally the class Fluid
    typedef Fluid<DIM>                     FluidMesh;
    /// Defines locally the class Element
    typedef typename FluidMesh::Elements   Elements;
    /// Defines locally the class Node
    typedef typename FluidMesh::Node       Nodes;
    /// Defines locally the class Boundary
    typedef typename FluidMesh::Boundaries Boundary;
   
private:
    FluidMesh          fluidModel;
    
    std::vector<Nodes *>     nodesFluid_;
    std::vector<Nodes *>     nodesPFluid_;
    //    std::vector<Nodes *>     nodesSolid_;
    std::vector<std::vector<Nodes *> > nodesSolid_;
    std::vector<Elements *>  elementsFluid_;
    std::vector<Boundary *>  boundaryFluid_;

    int numElemFluid;
    int numElemFluidBoundary;
    int numNodesFluid;
    int numNodesPFluid;
    int numNodesSolid;
    int numElemSolid;
    int numInterfaces;
    double dTime;
    static const int one=1;

    int rank, size;

    double pi = M_PI;

    std::vector<int>         groupInterfaces;

    std::pair<idx_t*,idx_t*> domDecompFluid; //Fluid Model Domain Decomposition

public:

    /// Sets the fluid and solid models and perform the preprocessing tasks
    /// @param Fluid fluid model @param char* solid input file
    void setFluidAndSolidModels(FluidMesh fluid, char *in_solid);
    
    /// Perform the preprocessing tasks
    void preProcess();

    /// Compute and store the element boxes for improving the element 
    /// correspondence searching process
    void setElementBoxes();

    /// Searchs the solid node correspondence into the fluid mesh
    /// @param int fluid boundary interface index @param int interface index
    void searchSolidNodeCorrespondence(int interface, int i);

    /// Searchs the solid node correspondence into the fluid mesh
    /// @param int fluid boundary interface index
    void searchFluidNodeCorrespondence(int interface);

    /// Updates the fluid mesh solving the Laplace problem
    void updateFluidMesh();

    /// Transfer solid velocity to the fluid nodes
    void transferSolidVelocity();

    /// Transfer fluid loads to the solid nodes
    void transferFluidLoad();

    /// Solves the partitioned weakly coupled Dirichlet-Neumann fluid-structure
    /// interaction problem 
    /// @param int number of time steps
    void solveFSIProblem(int numTimeSteps);

    /// Solves the partitioned strong coupled fixed-point block Gauss-Seidel
    /// with Aitken relaxation fluid-structure interaction problem 
    /// @param int number of time steps
    void solveFSIProblemGaussSeidel(int numTimeSteps);


};

//------------------------------------------------------------------------------
//--------------------------------IMPLEMENTATION--------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//-------------------COMPUTE NODAL CORRESPONDECE WITH ELEMENTS------------------
//------------------------------------------------------------------------------
template<>
void FSInteraction<2>::searchSolidNodeCorrespondence(int interface, int iSol){
    
    for (int isolid = 0; isolid < numNodesSolid; isolid++){

        typename Nodes::VecLocD xint, x_, deltaX, deltaXsi, x;
        typename Elements::Connectivity connec;
        ublas::bounded_vector<double,2>            xsi, xsiC;
        QuadShapeFunction<2>                       shapeQuad;
        typename QuadShapeFunction<2>::Values      phi_;
        typename Elements::DimMatrix               ainv;
        double xsiCC[3];
        std::pair<typename Elements::DimVector,typename Elements::DimVector> XK;
        int elemC;

        x = nodesSolid_[iSol][isolid] -> getCoordinates();
        
        elemC = 150000;
        xsiC(0) = 1.e50;
        xsiC(1) = 1.e50; 
        
        for (int ibound = 0; ibound < numElemFluidBoundary; ibound++){
            
            if (boundaryFluid_[ibound] -> getBoundaryGroup() == interface){
                
                int jel = boundaryFluid_[ibound] -> getElement();
                
                connec = elementsFluid_[jel] -> getConnectivity();
                
                //get boxes information        
                XK = elementsFluid_[jel] -> getXIntersectionParameter();
                
                //Chech if the node is inside the element box
                // if ((x(0) < XK.first(0)) || (x(0) > XK.second(0)) ||
                //     (x(1) < XK.first(1)) || (x(1) > XK.second(1))) continue;
                
                //Compute nodal correspondence
                xsiCC[0] = 1.e10;
                xsiCC[1] = 1.e10;
                xsiCC[2] = 1.e10;
                
                xsi(0) = 1. / 3.;
                xsi(1) = 1. / 3.;
                
                shapeQuad.evaluate(xsi,phi_);
                
                x_.clear();
                
                for (int i = 0; i < 6; i++){
                    xint = nodesFluid_[connec(i)] -> getCoordinates();
                    x_(0) += xint(0) * phi_(i);
                    x_(1) += xint(1) * phi_(i);                    
                };
                
                double error = 1.e6;
                
                xsi(0) = 1. / 3.;
                xsi(1) = 1. / 3.;
                
                int iterations = 0;
                
                while ((error > 1.e-8) && (iterations < 4)) {
                    
                    iterations++;
                    
                    deltaX = x - x_;
                    
                    deltaXsi.clear();
                    
                    ainv = elementsFluid_[jel] -> getJacobianMatrixValues(xsi);
                    
                    noalias(deltaXsi) = prod(trans(ainv),deltaX);
                    
                    xsi += deltaXsi;
                    
                    x_.clear();
                    
                    shapeQuad.evaluate(xsi,phi_);
                    
                    for (int i=0; i<6; i++){
                        xint = nodesFluid_[connec(i)] -> getCoordinates();
                        x_(0) += xint(0) * phi_(i);
                        x_(1) += xint(1) * phi_(i);   
                    };                   
                    error = norm_2(deltaXsi);
                };
                
                double t1 = -1.e-1;
                double t2 =  1. - t1;
                
                xsiCC[0] = xsi(0);
                xsiCC[1] = xsi(1);       
                xsiCC[2] = 1. - xsiCC[0] - xsiCC[1];
                
                if ((xsiCC[0] >= t1) && (xsiCC[1] >= t1) && (xsiCC[2] >= t1) &&
                    (xsiCC[0] <= t2) && (xsiCC[1] <= t2) && (xsiCC[2] <= t2)){
                    
                    xsiC = xsi;
                    elemC = jel;
                    break;
                };           
            };
        };

        nodesSolid_[iSol][isolid] -> setNodalCorrespondence(elemC,xsiC);
        
        //  std::cout << "isolid " << isolid << " " << interface << " " << elemC << " " << x(0) << " " << x(1) << std::endl;
    };
};

//------------------------------------------------------------------------------
//-------------------COMPUTE NODAL CORRESPONDECE WITH ELEMENTS------------------
//------------------------------------------------------------------------------
template<>
void FSInteraction<2>::searchFluidNodeCorrespondence(int interface){
    
    for (int ibound = 0; ibound < numElemFluidBoundary; ibound++){

        if (boundaryFluid_[ibound] -> getBoundaryGroup() == interface){    

            typename Nodes::VecLocD         x, xsi;
            typename Boundary::BoundConnect connec;
            int elemC;
            double xsiC;

            connec = boundaryFluid_[ibound] -> getBoundaryConnectivity();
            
            for (int inode = 0; inode < 3; inode++){
                x = nodesFluid_[connec(inode)] -> getCoordinates();
                
                searchcorrespondencefluid_(&x(0), &x(1), &elemC, &xsiC);
                xsi.clear();
                xsi(0) = xsiC;

                //std::cout << "asdasd " << elemC << " " << xsiC << std::endl; 

                nodesFluid_[connec(inode)] -> setNodalCorrespondence(elemC,xsi);
                
                //std::cout << "isolid " << connec(inode) << " " << elemC << " " << xsi(0) << std::endl;
            };
        };        
    };
};

//------------------------------------------------------------------------------
//-------------------------COMPUTE ELEMENT REGIONS/BOXES------------------------
//------------------------------------------------------------------------------
template<>
void FSInteraction<2>::setElementBoxes() {
    
    typename Elements::Connectivity connec;
    typename Nodes::VecLocD x1, x2, x3;
    ublas::bounded_vector<double,2> d1, d2, d3, xk, Xk;
    double dCk[3], dck[3];
    std::vector<ublas::bounded_vector<double,2> > di1, di2;

    di1.reserve(3);

    //Compute element boxes for coarse model
    //Only function for straight elements
    for (int jel = 0; jel < numElemFluid; jel++){
        connec = elementsFluid_[jel] -> getConnectivity();
        x1 = nodesFluid_[connec(0)] -> getCoordinates();
        x2 = nodesFluid_[connec(1)] -> getCoordinates();
        x3 = nodesFluid_[connec(2)] -> getCoordinates();      

        d1(0) = x2(1) - x1(1);
        d1(1) = x1(0) - x2(0);
        
        di1.push_back(d1);

        d2(0) = x3(1) - x2(1);
        d2(1) = x2(0) - x3(0);

        di1.push_back(d2);

        d3(0) = x1(1) - x3(1);
        d3(1) = x3(0) - x1(0);

        di1.push_back(d3);

        xk(0) = std::min(x1(0),std::min(x2(0), x3(0)));
        xk(1) = std::min(x1(1),std::min(x2(1), x3(1)));

        Xk(0) = std::max(x1(0),std::max(x2(0), x3(0)));
        Xk(1) = std::max(x1(1),std::max(x2(1), x3(1)));        

        dCk[0] = std::max(inner_prod(d1,x1),
                          std::max(inner_prod(d1,x2),inner_prod(d1,x3)));
        dck[0] = std::min(inner_prod(d1,x1),
                          std::min(inner_prod(d1,x2),inner_prod(d1,x3)));

        dCk[1] = std::max(inner_prod(d2,x1),
                          std::max(inner_prod(d2,x2),inner_prod(d2,x3)));
        dck[1] = std::min(inner_prod(d2,x1),
                          std::min(inner_prod(d2,x2),inner_prod(d2,x3)));

        dCk[2] = std::max(inner_prod(d3,x1),
                          std::max(inner_prod(d3,x2),inner_prod(d3,x3)));
        dck[2] = std::min(inner_prod(d3,x1),
                          std::min(inner_prod(d3,x2),inner_prod(d3,x3)));
        
        elementsFluid_[jel] -> setIntersectionParameters(xk, Xk, dCk, dck,di1);
    };

};

//------------------------------------------------------------------------------
//---------SETS FLUID AND SOLID MODELS AND GETS ITS BASIC INFORMATIONS----------
//------------------------------------------------------------------------------
template<>
void FSInteraction<2>::preProcess(){
    
    
    numInterfaces = fluidModel.getNumberofFSIInterfaces();
    
    std::cout << "Numero de interfaces " << numInterfaces << std::endl;

    groupInterfaces.reserve(numInterfaces);
    nodesSolid_.reserve(numInterfaces);
    //    int groupInterfaces2[numInterfaces];

    int fl=0;
    //Sets fluid elements and sides on interface boundaries
    for (int i=0; i<numElemFluidBoundary; i++){
        if ((boundaryFluid_[i] -> getConstrain(0) == 3) ||
            (boundaryFluid_[i] -> getConstrain(1) == 3)) {

            typename Boundary::BoundConnect connectB;
            connectB = boundaryFluid_[i] -> getBoundaryConnectivity();

            for (int j=0; j<numElemFluid; j++){
                typename Elements::Connectivity connect;
                connect = elementsFluid_[j] -> getConnectivity();
                
                int flag = 0;
                int side[3];
                for (int k=0; k<6; k++){
                    if ((connectB(0) == connect(k)) || 
                        (connectB(1) == connect(k)) ||
                        (connectB(2) == connect(k))){
                        side[flag] = k;
                        flag++;
                    };
                };
                
                if (flag == 3){
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
                            if (boundaryFluid_[i] -> getBoundaryGroup() == 
                                groupInterfaces[m]) fl2++;
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
        std::vector<Nodes *> j;

        nodesSolid_.push_back(j);
    
        nodesSolid_[k].reserve(numNodesSolid);
        int index = 0;
        
        for (int i=0; i<numNodesSolid; i++){
            typename Nodes::VecLocD x;
            int inode = i+1;
            getsolidposition_(&inode,&x(0),&x(1));
            
            Nodes *node = new Nodes(x,index++);
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
       
    domDecompFluid = fluidModel.getDomainDecomposition();

};

//------------------------------------------------------------------------------
//---------SETS FLUID AND SOLID MODELS AND GETS ITS BASIC INFORMATIONS----------
//------------------------------------------------------------------------------
template<>
void FSInteraction<2>::setFluidAndSolidModels(FluidMesh fluid, char *in_solid){
    
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &size);

    fluidModel = fluid;

    numElemFluid = fluidModel.elements_.size();
    numNodesFluid = fluidModel.nodes_.size();
    numElemFluidBoundary = fluidModel.boundary_.size();

    nodesFluid_  = fluidModel.nodes_;
    elementsFluid_ = fluidModel.elements_;
    boundaryFluid_ = fluidModel.boundary_;

    // Reads Solid input file
    preprocessing_(in_solid);
    
    getnumberofnodessolid_(&numNodesSolid);
    getnumberofelementssolid_(&numElemSolid);

    dTime = fluidModel.getTimeStep();

    std::cout << "NumElemSolid " << numElemSolid << std::endl;

    //Pre Processing data
    preProcess();


};

//------------------------------------------------------------------------------
//---------SETS FLUID AND SOLID MODELS AND GETS ITS BASIC INFORMATIONS----------
//------------------------------------------------------------------------------
template<>
void FSInteraction<2>::updateFluidMesh(){
    
                   


    for (int i = 0; i < numInterfaces; i++){

        int interf = groupInterfaces[i];       

        for (int ibound = 0; ibound < numElemFluidBoundary; ibound++){          
            if (boundaryFluid_[ibound] -> getBoundaryGroup() == interf){
                
                typename Boundary::BoundConnect connec;
                connec = boundaryFluid_[ibound] -> getBoundaryConnectivity();
              
                typename Nodes::VecLocD x, xsi, xant;

                for (int k=0; k<3; k++){

                    x.clear();
                    
                    int elem = nodesFluid_[connec(k)] ->
                        getNodalElemCorrespondence();
                    xsi = nodesFluid_[connec(k)] -> 
                        getNodalXsiCorrespondence();
                    
                    
                    if (rank == 0) 
                        getupdatedcoordinates_(&x(0),&x(1),&elem,&xsi(0));
              
                    MPI_Bcast(&x(0),1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
                    MPI_Bcast(&x(1),1,MPI_DOUBLE,0,PETSC_COMM_WORLD);

                    nodesFluid_[connec(k)] -> setUpdatedCoordinates(x);
                    // xant = nodesFluid_[connec(k)] -> getCoordinates();

                    //if(connec(k) == 29)std::cout << "Updated Coord " << x(0) << " " << x(1) << std::endl;
                };
                
            };//if interface
        };//ibound
    };//i


    //Solves Laplace Smoothing mesh moving scheme
    fluidModel.solveSteadyLaplaceProblem(1,1.e-6);
   
    for (int i = 0; i < numNodesFluid; i++){
        typename Nodes::VecLocD x, xp, up;
        double u[2];
            
        x = nodesFluid_[i] -> getCoordinates();
        xp = nodesFluid_[i] -> getPreviousCoordinates();
        up(0) = nodesFluid_[i] -> getPreviousMeshVelocity(0);
        up(1) = nodesFluid_[i] -> getPreviousMeshVelocity(1);
        
        u[0] = (x(0) - xp(0)) / dTime;//2. * (xp(0) - x(0)) / dTime - up(0);
        u[1] = (x(1) - xp(1)) / dTime;//2. * (xp(1) - x(1)) / dTime - up(1);
        
        nodesFluid_[i] -> setMeshVelocity(u);
    };


};


//------------------------------------------------------------------------------
//---------SETS FLUID AND SOLID MODELS AND GETS ITS BASIC INFORMATIONS----------
//------------------------------------------------------------------------------
template<>
void FSInteraction<2>::transferSolidVelocity(){
    

    for (int i = 0; i < numInterfaces; i++){

        int interf = groupInterfaces[i];       

        for (int ibound = 0; ibound < numElemFluidBoundary; ibound++){          
            if (boundaryFluid_[ibound] -> getBoundaryGroup() == interf){
                
                typename Boundary::BoundConnect connec;
                connec = boundaryFluid_[ibound] -> getBoundaryConnectivity();
              
                typename Nodes::VecLocD xsi,x,Acc;
                double u[2];

                for (int k=0; k<3; k++){
                    
                    int elem = nodesFluid_[connec(k)] ->
                        getNodalElemCorrespondence();
                    xsi = nodesFluid_[connec(k)] -> 
                        getNodalXsiCorrespondence();
                    
                    
                    if (rank == 0) 
                        getinterpolatedvelocity_(&u[0],&u[1],&elem,&xsi(0));
                        
                    MPI_Bcast(&u[0],1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
                    MPI_Bcast(&u[1],1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
               
                    nodesFluid_[connec(k)] -> setVelocity(u);
  
                    // x = nodesFluid_[connec(k)] -> getCoordinates();
                    // Acc(0) = nodesFluid_[connec(k)] -> getAcceleration(0);
                    // Acc(1) = nodesFluid_[connec(k)] -> getAcceleration(1);

                    // if(rank == 0) std::cout << "Accel " << Acc(0) << " " << Acc(1) << " " << std::endl; 
                };
            };
        };
    };

   

};

//------------------------------------------------------------------------------
//---------SETS FLUID AND SOLID MODELS AND GETS ITS BASIC INFORMATIONS----------
//------------------------------------------------------------------------------
template<>
void FSInteraction<2>::transferFluidLoad(){
    
    clearcouplingloads_();
    
    for (int iInterf = 0; iInterf < numInterfaces; iInterf++){
        for (int isolid = 0; isolid < numNodesSolid; isolid++){
            
            typename Nodes::VecLocD xsi, x;
            ublas::bounded_vector<double,2> load;

            int ielem = nodesSolid_[iInterf][isolid] ->
                getNodalElemCorrespondence();
            xsi = nodesSolid_[iInterf][isolid] -> getNodalXsiCorrespondence();
            
            load = elementsFluid_[ielem] -> getBoundaryLoad(xsi);

            int inode = isolid+1;
            setcouplingload_(&load(0),&load(1),&inode);
                       
        };//isolid
    };//iInterf
    
};


//------------------------------------------------------------------------------
//----------------SOLVES THE FLUID-STRUCTURE INTERACTION PROBLEM----------------
//------------------------------------------------------------------------------
template<>
void FSInteraction<2>::solveFSIProblem(int numTimeSteps){

    for (int iTimeStep = 0; iTimeStep < numTimeSteps; iTimeStep++){
        
        boost::posix_time::ptime t1 =                                 
            boost::posix_time::microsec_clock::local_time();

        if (rank == 0) {
            std::cout << std::endl;
            std::cout << "****************************************"
                      << "****************************************"
                      << std::endl;
            std::cout << 
                "                               TIME STEP = "
                      << iTimeStep << std::endl;
        };


        // // Updates fluid velocity profile
        for (int ibound = 0; ibound < numElemFluidBoundary; ibound++){
        
            Boundary::BoundConnect connectB;
            connectB = boundaryFluid_[ibound] -> getBoundaryConnectivity();
            int no1 = connectB(0);
            int no2 = connectB(1);
            int no3 = connectB(2);
        
            if (boundaryFluid_[ibound] -> getConstrain(0) == 1){
                
                double value = boundaryFluid_[ibound] -> getConstrainValue(0) * 
                    (1. - cos(0.4 * pi * dTime * iTimeStep));
                nodesFluid_[no1] -> setConstrains(0,boundaryFluid_[ibound] -> 
                                                  getConstrain(0),value);
                nodesFluid_[no2] -> setConstrains(0,boundaryFluid_[ibound] -> 
                                                  getConstrain(0),value);
                nodesFluid_[no3] -> setConstrains(0,boundaryFluid_[ibound] ->
                                                  getConstrain(0),value);
            };
            
        // };
            // if(boundaryFluid_[ibound] -> getBoundaryGroup() == 3){
            //     typename Nodes::VecLocD x1,x2,x3;
                
            //     x1 = nodesFluid_[no1] -> getCoordinates();
            //     x2 = nodesFluid_[no2] -> getCoordinates();
            //     x3 = nodesFluid_[no3] -> getCoordinates();

            //     double value1 = (1. - cos(0.4 * pi * dTime * iTimeStep)) * 
            //         (8. * x1(1) - 7.);
            //     double value2 = (1. - cos(0.4 * pi * dTime * iTimeStep)) * 
            //         (8. * x2(1) - 7.);
            //     double value3 = (1. - cos(0.4 * pi * dTime * iTimeStep)) * 
            //         (8. * x3(1) - 7.);

            //     nodesFluid_[no1] -> setConstrains(0,boundaryFluid_[ibound] -> 
            //                                       getConstrain(0),value1);
            //     nodesFluid_[no2] -> setConstrains(0,boundaryFluid_[ibound] -> 
            //                                       getConstrain(0),value2);
            //     nodesFluid_[no3] -> setConstrains(0,boundaryFluid_[ibound] ->
            //                                       getConstrain(0),value3);
            //     nodesFluid_[no1] -> setConstrains(1,boundaryFluid_[ibound] -> 
            //                                       getConstrain(1),0);
            //     nodesFluid_[no2] -> setConstrains(1,boundaryFluid_[ibound] -> 
            //                                       getConstrain(1),0);
            //     nodesFluid_[no3] -> setConstrains(1,boundaryFluid_[ibound] ->
            //                                       getConstrain(1),0);

            // };

               };  

        // Solves solid

        if (rank == 0) transferFluidLoad();
        
        if (rank == 0) updateqsrs_();
        
        if (rank == 0) solveframestructure_(&iTimeStep);

        if (rank == 0) updatesolid_(&iTimeStep);



        for (int i = 0; i < numNodesFluid; i++){
            typename Nodes::VecLocD x, um;
            
            x = nodesFluid_[i] -> getCoordinates();
            um(0) = nodesFluid_[i] -> getMeshVelocity(0);
            um(1) = nodesFluid_[i] -> getMeshVelocity(1);

            nodesFluid_[i] -> setPreviousCoordinates(0,x(0));
            nodesFluid_[i] -> setPreviousCoordinates(1,x(1));

            nodesFluid_[i] -> setPreviousMeshVelocity(0,um(0));
            nodesFluid_[i] -> setPreviousMeshVelocity(1,um(1));           
        };

        
        updateFluidMesh();
        
     for (int i = 0; i < numNodesFluid; i++){
            double accel[2], u[2], uprev[2];
            
            //Compute acceleration
            u[0] = nodesFluid_[i] -> getVelocity(0);
            u[1] = nodesFluid_[i] -> getVelocity(1);
            
            uprev[0] = nodesFluid_[i] -> getPreviousVelocity(0);
            uprev[1] = nodesFluid_[i] -> getPreviousVelocity(1);
            
            accel[0] = (u[0] - uprev[0]) / dTime;
            accel[1] = (u[1] - uprev[1]) / dTime;
            
            nodesFluid_[i] -> setAcceleration(accel);
            
            //Updates velocity
            nodesFluid_[i] -> setPreviousVelocity(u);
            
        };

        transferSolidVelocity();
        
        fluidModel.solveFSIFluid(3, 1.e-5, 2);
        
        boost::posix_time::ptime t2 =                                   \
            boost::posix_time::microsec_clock::local_time();

        if (rank == 0) {
            boost::posix_time::time_duration diff = t2 - t1;
            std::cout << "****************************************"
                      << "****************************************"
                      << std::endl;
            std::cout << "********************** PROCESSING TIME = " << 
                std::fixed << diff.total_milliseconds()/1000.
                      << " seconds **********************" 
                       << std::endl;
            std::cout << "****************************************"
                      << "****************************************"
                      << std::endl;
        };

        if (rank == 0) {
            //            if(iTimeStep % 5 == 0){
            fluidModel.printVelocity(iTimeStep);
            printstructure_();
            //            };
        };




    };//Time Steps
};

//------------------------------------------------------------------------------
//----------------SOLVES THE FLUID-STRUCTURE INTERACTION PROBLEM----------------
//------------------------------------------------------------------------------
template<>
void FSInteraction<2>::solveFSIProblemGaussSeidel(int numTimeSteps){

    std::string om = "omega.txt";
    std::ofstream saidaOmega(om.c_str());


    double sizeSolid = 3 * numNodesSolid;

    ublas::vector<double> X_k(sizeSolid), Y_k(sizeSolid), deltaXi(sizeSolid),
        deltaXii(sizeSolid);
    

    double omega = 1.;
    double mu = 0.;

    X_k.clear();Y_k.clear();deltaXi.clear();deltaXii.clear();

    for (int iTimeStep = 0; iTimeStep < numTimeSteps; iTimeStep++){  

        boost::posix_time::ptime t1 =                                 
            boost::posix_time::microsec_clock::local_time();

        if (rank == 0) {
            std::cout << std::endl;
            std::cout << "****************************************"
                      << "****************************************"
                      << std::endl;
            std::cout << 
                "                               TIME STEP = "
                      << iTimeStep << std::endl;
        };
        
        //SOMENTE PARA EXEMPLO DA CAVIDADE - INICIO
        for (int ibound = 0; ibound < numElemFluidBoundary; ibound++){
            
            Boundary::BoundConnect connectB;
            connectB = boundaryFluid_[ibound] -> getBoundaryConnectivity();
            int no1 = connectB(0);
            int no2 = connectB(1);
            int no3 = connectB(2);
            
            if (boundaryFluid_[ibound] -> getConstrain(0) == 1){
                
                double value = boundaryFluid_[ibound] -> getConstrainValue(0) * 
                    (1. - cos(0.4 * pi * dTime * iTimeStep));
                nodesFluid_[no1] -> setConstrains(0,boundaryFluid_[ibound] -> 
                                                  getConstrain(0),value);
                nodesFluid_[no2] -> setConstrains(0,boundaryFluid_[ibound] -> 
                                                  getConstrain(0),value);
                nodesFluid_[no3] -> setConstrains(0,boundaryFluid_[ibound] ->
                                                  getConstrain(0),value);
            };
        };
        //SOMENTE PARA EXEMPLO DA CAVIDADE - FIM

        for (int i = 0; i < numNodesFluid; i++){
            double accel[2], u[2], uprev[2];
            //Compute acceleration
            u[0] = nodesFluid_[i] -> getVelocity(0);
            u[1] = nodesFluid_[i] -> getVelocity(1);
            
            uprev[0] = nodesFluid_[i] -> getPreviousVelocity(0);
            uprev[1] = nodesFluid_[i] -> getPreviousVelocity(1);
            
            accel[0] = (u[0] - uprev[0]) / dTime;
            accel[1] = (u[1] - uprev[1]) / dTime;
            
            nodesFluid_[i] -> setAcceleration(accel);
            
            //Updates velocity
            nodesFluid_[i] -> setPreviousVelocity(u);
        };

        if (rank == 0) updateqsrs_();

        for (int i = 0; i < numNodesFluid; i++){
            typename Nodes::VecLocD x;
            
            x = nodesFluid_[i] -> getCoordinates();
            nodesFluid_[i] -> setPreviousCoordinates(0,x(0));
            nodesFluid_[i] -> setPreviousCoordinates(1,x(1));
        };

        for (int i = 0; i < numNodesSolid; i++){
            int dof = 3*i+1;
            if (rank == 0) getposition_(&dof,&Y_k(3*i  ));
            dof++;
            if (rank == 0) getposition_(&dof,&Y_k(3*i+1));
            MPI_Bcast(&Y_k(3*i  ),1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
            MPI_Bcast(&Y_k(3*i+1),1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
        };
      
        saidaOmega << std::endl << "Passo de tempo " << iTimeStep << std::endl;

        //COMPUTING PREDICTOR
        for (int i = 0; i < numNodesSolid; i++){
            if (rank == 0) {
                double v_ = 0.;
                double v_prev = 0.;
                int dof = 3*i+1;
                getposition_(&dof,&Y_k(3*i  ));
                getvelocity_(&dof,&v_);
                getpreviousvelocity_(&dof,&v_prev);

                Y_k(3*i  ) += dTime * (1.5 * v_ - 0.5 * v_prev);

                setposition_(&dof,&Y_k(3*i  ));
                
                dof++;

                getposition_(&dof,&Y_k(3*i+1));
                getvelocity_(&dof,&v_);
                getpreviousvelocity_(&dof,&v_prev);

                Y_k(3*i+1) += dTime * (1.5 * v_ - 0.5 * v_prev);

                setposition_(&dof,&Y_k(3*i+1));
            };
            MPI_Bcast(&Y_k(3*i  ),1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
            MPI_Bcast(&Y_k(3*i+1),1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
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

            X_k = Y_k;
            
            updateFluidMesh();
             
            transferSolidVelocity();
            
            fluidModel.solveFSIFluid(2, 1.e-5, 2);
            
            if (rank == 0) transferFluidLoad();
            
            if (rank == 0) solveframestructure_(&iTimeStep);
            
            
            for (int i = 0; i < numNodesSolid; i++){
                int dof = 3*i+1;
                if (rank == 0) getposition_(&dof,&Y_k(3*i  ));
                dof++;
                if (rank == 0) getposition_(&dof,&Y_k(3*i+1));
                MPI_Bcast(&Y_k(3*i  ),1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
                MPI_Bcast(&Y_k(3*i+1),1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
            };
            
            deltaXii = X_k - Y_k;
            
            //Aitken Relaxation
            if (iterations > 0){
                mu = mu + (mu - 1.0) * inner_prod(deltaXi-deltaXii,deltaXii) / 
                    norm_2(deltaXi-deltaXii);
              
            };        
            
            deltaXi = deltaXii;
            residual = norm_2(deltaXii);

            omega = 1. - mu;

            Y_k = (1. - omega) * X_k + omega * Y_k;
            

            for (int i = 0; i < numNodesSolid; i++){
                if (rank == 0) {
                    int dof = 3*i+1;
                    setposition_(&dof,&Y_k(3*i  ));
                    dof++;
                    setposition_(&dof,&Y_k(3*i+1));
                };
                MPI_Bcast(&Y_k(3*i  ),1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
                MPI_Bcast(&Y_k(3*i+1),1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
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
                
        boost::posix_time::ptime t2 =                                   \
            boost::posix_time::microsec_clock::local_time();

        if (rank == 0) {
            boost::posix_time::time_duration diff = t2 - t1;
            std::cout << "****************************************"
                      << "****************************************"
                      << std::endl;
            std::cout << "********************** PROCESSING TIME = " << 
                std::fixed << diff.total_milliseconds()/1000.
                      << " seconds **********************" 
                       << std::endl;
            std::cout << "****************************************"
                      << "****************************************"
                      << std::endl;
        };

        // Printing Results
        if (rank == 0) {
            if(iTimeStep % 1 == 0){
                fluidModel.printVelocity(iTimeStep);
                printstructure_();
            };
        };


    };//Time Steps
};


#endif
