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

#include "Arlequin.hpp"

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
    typedef Fluid<DIM>                     FluidModel;
    /// Defines locally the class Element
    typedef typename FluidModel::Elements   Elements;
    /// Defines locally the class Node
    typedef typename FluidModel::Node       Nodes;
    /// Defines locally the class Boundary
    typedef typename FluidModel::Boundaries Boundary;

    /// Defines locally the class Arlequin
    typedef Arlequin<DIM>                   ArlequinModel;
   
private:
    FluidModel         fluidModel;
    ArlequinModel      arlequinModel;
    
    std::vector<Nodes *>     nodesFluid_;
    //    std::vector<Nodes *>     nodesSolid_;
    std::vector<std::vector<Nodes *> > nodesSolid_;
    std::vector<Elements *>  elementsFluid_;
    std::vector<Boundary *>  boundaryFluid_;

    std::vector<Nodes *>     nodesArlequinCoarse_;
    std::vector<Nodes *>     nodesArlequinFine_;
    std::vector<Elements *>  elementsArlequinCoarse_;
    std::vector<Elements *>  elementsArlequinFine_;
    std::vector<Boundary *>  boundaryArlequinCoarse_;
    std::vector<Boundary *>  boundaryArlequinFine_;

    //Fluid variables
    int numElemFluid;
    int numElemFluidBoundary;
    int numNodesFluid;
    int numNodesSolid;
    int numElemSolid;
    int numInterfaces;
    double dTime;
    static const int one=1;

    //Arlequin variables
    int numElemArlequinFine;
    int numElemArlequinCoarse;
    int numNodesArlequinFine;
    int numNodesArlequinCoarse;
    int numElemArlequinBoundaryCoarse;
    int numElemArlequinBoundaryFine;
   


    int rank, size;

    double pi = M_PI;

    std::vector<int>         groupInterfaces;

    std::pair<idx_t*,idx_t*> domDecompFluid; //Fluid Model Domain Decomposition

    //Arlequin Model Domain Decomposition
    std::pair<idx_t*,idx_t*> domDecompArlequinFine; 
    std::pair<idx_t*,idx_t*> domDecompArlequinCoarse; 

public:

    /// Sets the fluid and solid models and perform the preprocessing tasks
    /// @param Fluid fluid model @param char* solid input file
    void setFluidAndSolidModels(FluidModel fluid, char *in_solid);

    /// Sets the Arlequin and solid models and perform the preprocessing tasks
    /// @param Arlequin Arlequin fluid model @param char* solid input file
    void setArlequinAndSolidModels(ArlequinModel arlq, char *in_solid);
    
    /// Perform the fluid preprocessing tasks
    void preProcessFluid();

    /// Perform the fluid preprocessing tasks
    void preProcessArlequin();

    /// Compute and store the element boxes for improving the element 
    /// correspondence searching process
    void setElementBoxes();

    /// Searchs the solid node correspondence into the fluid mesh
    /// @param int fluid boundary interface index @param int interface index
    void searchSolidNodeCorrespondence(int interface, int i);

    /// Searchs the solid node correspondence into the Arlequin model
    /// @param int fluid boundary interface index @param int interface index
    void searchSolidNodeCorrespondenceArlequin(int interface, int i);

    /// Searchs the solid node correspondence into the fluid mesh
    /// @param int fluid boundary interface index
    void searchFluidNodeCorrespondence(int interface);

    /// Searchs the solid node correspondence into the Arlequin model
    /// @param int Arlequin boundary interface index
    void searchArlequinNodeCorrespondence(int interface);

    /// Updates the fluid mesh solving the Laplace problem
    void updateFluidMesh();

    /// Updates the Arlequin fine mesh solving the Laplace problem
    void updateArlequinMesh();

    /// Transfer solid velocity to the fluid nodes
    void transferSolidVelocity();

    /// Transfer solid velocity to the Arlequin model nodes
    void transferSolidVelocityArlequin();

    /// Transfer fluid loads to the solid nodes
    void transferFluidLoad();

    /// Transfer Arlequin model fluid loads to the solid nodes
    void transferArlequinLoad();

    /// Solves the partitioned weakly coupled Dirichlet-Neumann fluid-structure
    /// interaction problem 
    /// @param int number of time steps
    void solveFSIProblem(int numTimeSteps);

    /// Solves the partitioned strong coupled fixed-point block Gauss-Seidel
    /// with Aitken relaxation fluid-structure interaction problem 
    /// @param int number of time steps
    void solveFSIProblemGaussSeidel(int numTimeSteps);

    /// Solves the partitioned strong coupled fixed-point block Gauss-Seidel
    /// with Aitken relaxation fluid-structure interaction problem with the
    /// Arlequin method 
    /// @param int number of time steps
    void solveFSIProblemGaussSeidelArlequin(int numTimeSteps);


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

        int *connec;
        QuadShapeFunction<2>                       shapeQuad;
        double phi_[6] = {};
        double ainv[2][2] = {};
        double xsiCC[3];
        // std::pair<double*,double*> XK;
        int elemC;

        double* x = nodesSolid_[iSol][isolid] -> getCoordinates();
        
        elemC = 150000;
        double xsiC[2];
        xsiC[0] = 1.e50;
        xsiC[1] = 1.e50; 
        double xsi[2];
        double x_[2];
        double deltaX[2];
        double deltaXsi[2];
        
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
                
                xsi[0] = 1. / 3.;
                xsi[1] = 1. / 3.;
                
                shapeQuad.evaluate(xsi,phi_);
                
                x_[0] = 0.;
                x_[1] = 0.;
                
                for (int i = 0; i < 6; i++){
                    double* xint = nodesFluid_[connec[i]] -> getCoordinates();
                    x_[0] += xint[0] * phi_[i];
                    x_[1] += xint[1] * phi_[i];                    
                };
                
                double error = 1.e6;
                
                xsi[0] = 1. / 3.;
                xsi[1] = 1. / 3.;
                
                int iterations = 0;
                
                while ((error > 1.e-8) && (iterations < 4)) {
                    
                    iterations++;
                    
                    deltaX[0] = x[0] - x_[0];
                    deltaX[1] = x[1] - x_[1];
                    
                    deltaXsi[0] = 0.;
                    deltaXsi[1] = 0.;
                    
                    elementsFluid_[jel] -> getJacobianMatrix(xsi,ainv);
                    
                    // noalias(deltaXsi) = prod(trans(ainv),deltaX);

                    for (int i = 0; i < 2; i++)
                        for (int j = 0; j < 2; j++)
                            deltaXsi[i] += ainv[j][i] * deltaX[j];


                    xsi[0] += deltaXsi[0];
                    xsi[1] += deltaXsi[1];
                    
                    x_[0] = 0.;
                    x_[1] = 0.; 
                    
                    shapeQuad.evaluate(xsi,phi_);
                    
                    for (int i=0; i<6; i++){
                        double* xint = nodesFluid_[connec[i]] -> getCoordinates();
                        x_[0] += xint[0] * phi_[i];
                        x_[1] += xint[1] * phi_[i];   
                    };                   
                    error = sqrt(deltaXsi[0]*deltaXsi[0] + deltaXsi[1]*deltaXsi[1]);
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
                    break;
                };           
            };
        };

        nodesSolid_[iSol][isolid] -> setNodalCorrespondence(elemC,xsiC);
        
        //  std::cout << "isolid " << isolid << " " << interface << " " << elemC << " " << x(0) << " " << x(1) << " " << xsiC(0) << " " << xsiC(1) << std::endl;
    };
};

//------------------------------------------------------------------------------
//-------------------COMPUTE NODAL CORRESPONDECE WITH ELEMENTS------------------
//------------------------------------------------------------------------------
template<>
void FSInteraction<2>::searchSolidNodeCorrespondenceArlequin(int interface, 
                                                             int iSol){
    
    for (int isolid = 0; isolid < numNodesSolid; isolid++){

        int *connec;
        QuadShapeFunction<2>                       shapeQuad;
        double phi_[6] = {};
        double ainv[2][2] = {};
        double xsiCC[3];
        // std::pair<typename Elements::DimVector,typename Elements::DimVector> XK;
        int elemC;

        double* x = nodesSolid_[iSol][isolid] -> getCoordinates();
        
        elemC = 150000;
        double xsiC[2];
        xsiC[0] = 1.e50;
        xsiC[1] = 1.e50; 
        double xsi[2];
        double x_[2];
        double deltaX[2];
        double deltaXsi[2];
        
        for (int ibound = 0; ibound < numElemArlequinBoundaryFine; ibound++){
            
            if (boundaryArlequinFine_[ibound] -> 
                getBoundaryGroup() == interface){
                
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
                
                xsi[0] = 1. / 3.;
                xsi[1] = 1. / 3.;
                
                shapeQuad.evaluate(xsi,phi_);
                
                x_[0] = 0.;
                x_[1] = 0.;
                
                for (int i = 0; i < 6; i++){
                    double *xint = nodesArlequinFine_[connec[i]] -> getCoordinates();
                    x_[0] += xint[0] * phi_[i];
                    x_[1] += xint[1] * phi_[i];                    
                };
                
                double error = 1.e6;
                
                xsi[0] = 1. / 3.;
                xsi[1] = 1. / 3.;
                
                int iterations = 0;
                
                while ((error > 1.e-8) && (iterations < 4)) {
                    
                    iterations++;
                    
                    deltaX[0] = x[0] - x_[0];
                    deltaX[1] = x[1] - x_[1];
                    
                    deltaXsi[0] = 0.;
                    deltaXsi[1] = 0.;
                    
                    elementsArlequinFine_[jel] -> getJacobianMatrix(xsi,ainv);
                    
                    // noalias(deltaXsi) = prod(trans(ainv),deltaX);

                    for (int i = 0; i < 2; i++)
                        for (int j = 0; j < 2; j++)
                            deltaXsi[i] += ainv[j][i] * deltaX[j];
                    
                    xsi[0] += deltaXsi[0];
                    xsi[1] += deltaXsi[1];
                    
                    x_[0] = 0.;
                    x_[1] = 0.;
                    
                    shapeQuad.evaluate(xsi,phi_);
                    
                    for (int i=0; i<6; i++){
                        double* xint = nodesArlequinFine_[connec[i]] -> getCoordinates();
                        x_[0] += xint[0] * phi_[i];
                        x_[1] += xint[1] * phi_[i];   
                    };                   
                    error = sqrt(deltaXsi[0]*deltaXsi[0] + deltaXsi[1]*deltaXsi[1]);
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
                    break;
                };           
            };
        };

        nodesSolid_[iSol][isolid] -> setNodalCorrespondence(elemC,xsiC);
        
        // if (rank == 0) std::cout << "isolid " << isolid << " " << interface << " " << elemC << " " << x(0) << " " << x(1) << " " << xsiC(0) << " " << xsiC(1) << std::endl;
    };
};

//------------------------------------------------------------------------------
//-------------------COMPUTE NODAL CORRESPONDECE WITH ELEMENTS------------------
//------------------------------------------------------------------------------
template<>
void FSInteraction<2>::searchFluidNodeCorrespondence(int interface){
    
    for (int ibound = 0; ibound < numElemFluidBoundary; ibound++){

        if (boundaryFluid_[ibound] -> getBoundaryGroup() == interface){    
            typename Boundary::BoundConnect connec;
            int elemC;
            double xsiC;

            connec = boundaryFluid_[ibound] -> getBoundaryConnectivity();
            
            for (int inode = 0; inode < 3; inode++){
                double* x = nodesFluid_[connec(inode)] -> getCoordinates();
                
                searchcorrespondencefluid_(&x[0], &x[1], &elemC, &xsiC);
                double xsi[2] = {};
                xsi[0] = xsiC;

                //std::cout << "asdasd " << elemC << " " << xsiC << std::endl; 

                nodesFluid_[connec(inode)] -> setNodalCorrespondence(elemC,xsi);
                
                // std::cout << "isolid " << connec(inode) << " " << elemC << " " << xsi(0) << std::endl;
            };
        };        
    };
};

//------------------------------------------------------------------------------
//-------------------COMPUTE NODAL CORRESPONDECE WITH ELEMENTS------------------
//------------------------------------------------------------------------------
template<>
void FSInteraction<2>::searchArlequinNodeCorrespondence(int interface){
    
    for (int ibound = 0; ibound < numElemArlequinBoundaryFine; ibound++){

        if (boundaryArlequinFine_[ibound] -> getBoundaryGroup() == interface){

            typename Boundary::BoundConnect connec;
            int elemC;
            double xsiC;

            connec = boundaryArlequinFine_[ibound] -> getBoundaryConnectivity();
            
            for (int inode = 0; inode < 3; inode++){
                double* x = nodesArlequinFine_[connec(inode)] -> getCoordinates();
                
                searchcorrespondencefluid_(&x[0], &x[1], &elemC, &xsiC);
                double xsi[2] = {};
                xsi[0] = xsiC;

                // if (rank == 0) std::cout << "asdasd " << elemC << " " << xsiC << " " << interface << std::endl; 

                nodesArlequinFine_[connec(inode)] -> setNodalCorrespondence(elemC,xsi);
                
                // std::cout << "isolid " << connec(inode) << " " << elemC << " " << xsi(0) << std::endl;
            };
        };        
    };
};

//------------------------------------------------------------------------------
//-------------------------COMPUTE ELEMENT REGIONS/BOXES------------------------
//------------------------------------------------------------------------------
template<>
void FSInteraction<2>::setElementBoxes() {
    
    int *connec;
    double xk[2], Xk[2];
    double dCk[3], dck[3];
    std::vector<ublas::bounded_vector<double,2> > di1, di2;

    di1.reserve(3);

    //Compute element boxes for coarse model
    //Only function for straight elements
    for (int jel = 0; jel < numElemFluid; jel++){
        connec = elementsFluid_[jel] -> getConnectivity();
        double* x1 = nodesFluid_[connec[0]] -> getCoordinates();
        double* x2 = nodesFluid_[connec[1]] -> getCoordinates();
        double* x3 = nodesFluid_[connec[2]] -> getCoordinates();      

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
template<>
void FSInteraction<2>::preProcessFluid(){
    
    
    numInterfaces = fluidModel.getNumberofFSIInterfaces();
    
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

            typename Boundary::BoundConnect connectB;
            connectB = boundaryFluid_[i] -> getBoundaryConnectivity();

            for (int j=0; j<numElemFluid; j++){
                int *connect;
                connect = elementsFluid_[j] -> getConnectivity();
                
                int flag = 0;
                int side[3];
                for (int k=0; k<6; k++){
                    if ((connectB(0) == connect[k]) || 
                        (connectB(1) == connect[k]) ||
                        (connectB(2) == connect[k])){
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
            double x[2];
            int inode = i+1;
            getsolidposition_(&inode,&x[0],&x[1]);
            
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
void FSInteraction<2>::preProcessArlequin(){
    
    
    numInterfaces = arlequinModel.fineModel.getNumberofFSIInterfaces();
    
    std::cout << "Number of interfaces " << numInterfaces << std::endl;

    groupInterfaces.reserve(numInterfaces);
    nodesSolid_.reserve(numInterfaces);

    int fl=0;
    //Sets fluid elements and sides on interface boundaries
    for (int i=0; i<numElemArlequinBoundaryFine; i++){
        if ((boundaryArlequinFine_[i] -> getConstrain(0) == 3) ||
            (boundaryArlequinFine_[i] -> getConstrain(1) == 3)) {

            typename Boundary::BoundConnect connectB;
            connectB = boundaryArlequinFine_[i] -> getBoundaryConnectivity();

            for (int j=0; j<numElemArlequinFine; j++){
                
                int *connect = elementsArlequinFine_[j] -> getConnectivity();
                
                int flag = 0;
                int side[3];
                for (int k=0; k<6; k++){
                    if ((connectB(0) == connect[k]) || 
                        (connectB(1) == connect[k]) ||
                        (connectB(2) == connect[k])){
                        side[flag] = k;
                        flag++;
                    };
                };
                
                if (flag == 3){
                    boundaryArlequinFine_[i] -> setElement(j);
                    elementsArlequinFine_[j] -> setFSIInterface();
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

                    //Sets element index and side
                    if ((side[0]==4) || (side[1]==4) || (side[2]==4)){
                        boundaryArlequinFine_[i] -> setElementSide(0);
                        elementsArlequinFine_[boundaryArlequinFine_[i] -> 
                                              getElement()] -> 
                            setElemSideInBoundary(0);
                    };
                    if ((side[0]==5) || (side[1]==5) || (side[2]==5)){
                        boundaryArlequinFine_[i] -> setElementSide(1);
                        elementsArlequinFine_[boundaryArlequinFine_[i] ->
                                              getElement()] -> 
                            setElemSideInBoundary(1);
                    };
                    if ((side[0]==3) || (side[1]==3) || (side[2]==3)){
                        boundaryArlequinFine_[i] -> setElementSide(2);
                        elementsArlequinFine_[boundaryArlequinFine_[i] ->
                                              getElement()] -> 
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
            double x[2];
            int inode = i+1;
            getsolidposition_(&inode,&x[0],&x[1]);
            
            Nodes *node = new Nodes(x,index++);
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
       
    domDecompArlequinCoarse =arlequinModel.coarseModel.getDomainDecomposition();
    domDecompArlequinFine = arlequinModel.fineModel.getDomainDecomposition();

};


//------------------------------------------------------------------------------
//---------SETS FLUID AND SOLID MODELS AND GETS ITS BASIC INFORMATIONS----------
//------------------------------------------------------------------------------
template<>
void FSInteraction<2>::setFluidAndSolidModels(FluidModel fluid, char *in_solid){
    
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

    if(rank == 0) std::cout << "NumElemSolid " << numElemSolid << std::endl;

    //Pre Processing data
    preProcessFluid();

};

//------------------------------------------------------------------------------
//---------SETS FLUID AND SOLID MODELS AND GETS ITS BASIC INFORMATIONS----------
//------------------------------------------------------------------------------
template<>
void FSInteraction<2>::setArlequinAndSolidModels(ArlequinModel arlq, 
                                                 char *in_solid){
    
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &size);

    arlequinModel = arlq;

    numElemArlequinCoarse = arlequinModel.elementsCoarse_.size();
    numElemArlequinFine = arlequinModel.elementsFine_.size();
    numNodesArlequinCoarse = arlequinModel.nodesCoarse_.size();
    numNodesArlequinFine = arlequinModel.nodesFine_.size();
    numElemArlequinBoundaryCoarse = arlequinModel.boundaryCoarse_.size();
    numElemArlequinBoundaryFine = arlequinModel.boundaryFine_.size();
    
    nodesArlequinCoarse_ = arlequinModel.nodesCoarse_;
    nodesArlequinFine_ = arlequinModel.nodesFine_;
    elementsArlequinCoarse_ = arlequinModel.elementsCoarse_;
    elementsArlequinFine_ = arlequinModel.elementsFine_;
    boundaryArlequinCoarse_ = arlequinModel.boundaryCoarse_;
    boundaryArlequinFine_ = arlequinModel.boundaryFine_;

    // Reads Solid input file
    preprocessing_(in_solid);

    getnumberofnodessolid_(&numNodesSolid);
    getnumberofelementssolid_(&numElemSolid);

    dTime = arlequinModel.fineModel.getTimeStep();

    if(rank == 0) std::cout << "NumElemSolid " << numElemSolid << std::endl;

    //Pre Processing data
    preProcessArlequin();

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
              
                for (int k=0; k<3; k++){

                    double x[2];
                    
                    int elem = nodesFluid_[connec(k)] -> getNodalElemCorrespondence();
                    double* xsi = nodesFluid_[connec(k)] -> getNodalXsiCorrespondence();
                    
                    
                    if (rank == 0) getupdatedcoordinates_(&x[0],&x[1],&elem,&xsi[0]);
              
                    MPI_Bcast(&x[0],1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
                    MPI_Bcast(&x[1],1,MPI_DOUBLE,0,PETSC_COMM_WORLD);

                    nodesFluid_[connec(k)] -> setUpdatedCoordinates(x);
                    nodesFluid_[connec(k)] -> setCoordinates(x);

                    //if(connec(k) == 29)std::cout << "Updated Coord " << x(0) << " " << x(1) << std::endl;
                };
                
            };//if interface
        };//ibound
    };//i


    //Solves Laplace Smoothing mesh moving scheme
    fluidModel.solveSteadyLaplaceProblem(5,1.e-6);
   
    for (int i = 0; i < numNodesFluid; i++){
        double u[2], up[2];
            
        double *x = nodesFluid_[i] -> getCoordinates();
        double *xp = nodesFluid_[i] -> getPreviousCoordinates();
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
template<>
void FSInteraction<2>::updateArlequinMesh(){

    // std::cout << "AQUI6.1 " << rank << std::endl;
    // MPI_Barrier(PETSC_COMM_WORLD);

    for (int i = 0; i < numInterfaces; i++){

        int interf = groupInterfaces[i];       

        for (int ibound = 0; ibound < numElemArlequinBoundaryFine; ibound++){
            if (boundaryArlequinFine_[ibound] -> getBoundaryGroup() == interf){

                typename Boundary::BoundConnect connec;
                connec = boundaryArlequinFine_[ibound] -> getBoundaryConnectivity();
              
                for (int k=0; k<3; k++){

                    double x[2];
                    
                    int elem = nodesArlequinFine_[connec(k)] -> getNodalElemCorrespondence();
                    double* xsi = nodesArlequinFine_[connec(k)] -> getNodalXsiCorrespondence();

                    // if (rank == 0) std::cout << "AQUI6.1.1 " << rank << " " << elem << " " << xsi << " " << x(0) << " " << x(1) << std::endl;
                    
                    if (rank == 0) getupdatedcoordinates_(&x[0],&x[1],&elem,&xsi[0]);

                    MPI_Bcast(&x[0],1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
                    MPI_Bcast(&x[1],1,MPI_DOUBLE,0,PETSC_COMM_WORLD);

                    // if (rank == 0) std::cout << "AQUI6.1.2 " << rank << " " << elem << " " << xsi << " " << x(0) << " " << x(1) << std::endl;
                    MPI_Barrier(PETSC_COMM_WORLD);

                    nodesArlequinFine_[connec(k)] -> setUpdatedCoordinates(x);
                    nodesArlequinFine_[connec(k)] -> setCoordinates(x);

                    // std::cout << "Updated Coord " << x(0) << " " << x(1) << std::endl;
                };
                
            };//if interface
        };//ibound
    };//i

    // std::cout << "AQUI6.2 " << rank << std::endl;
    // MPI_Barrier(PETSC_COMM_WORLD);


    //Solves Laplace Smoothing mesh moving scheme
    arlequinModel.fineModel.solveSteadyLaplaceProblem(1,1.e-4);
   
    // std::cout << "AQUI6.3 " << rank << std::endl;
    // MPI_Barrier(PETSC_COMM_WORLD);

    for (int i = 0; i < numNodesArlequinFine; i++){
        double up[2];
        double u[2];
            
        double* x = nodesArlequinFine_[i] -> getCoordinates();
        double* xp = nodesArlequinFine_[i] -> getPreviousCoordinates();
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
template<>
void FSInteraction<2>::transferSolidVelocity(){
    

    for (int i = 0; i < numInterfaces; i++){

        int interf = groupInterfaces[i];       

        for (int ibound = 0; ibound < numElemFluidBoundary; ibound++){          
            if (boundaryFluid_[ibound] -> getBoundaryGroup() == interf){
                
                typename Boundary::BoundConnect connec;
                connec = boundaryFluid_[ibound] -> getBoundaryConnectivity();
              
                double u[2];

                for (int k=0; k<3; k++){
                    
                    int elem = nodesFluid_[connec(k)] -> getNodalElemCorrespondence();
                    double* xsi = nodesFluid_[connec(k)] -> getNodalXsiCorrespondence();
                    
                    
                    if (rank == 0) 
                        getinterpolatedvelocity_(&u[0],&u[1],&elem,&xsi[0]);
                        
                    MPI_Bcast(&u[0],1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
                    MPI_Bcast(&u[1],1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
               
                    nodesFluid_[connec(k)] -> setVelocity(u);
  
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
template<>
void FSInteraction<2>::transferSolidVelocityArlequin(){
    

    for (int i = 0; i < numInterfaces; i++){

        int interf = groupInterfaces[i];       

        for (int ibound = 0; ibound < numElemArlequinBoundaryFine; ibound++){
            if (boundaryArlequinFine_[ibound] -> getBoundaryGroup() == interf){
                
                typename Boundary::BoundConnect connec;
                connec = boundaryArlequinFine_[ibound] -> 
                    getBoundaryConnectivity();
              
                double u[2];

                for (int k=0; k<3; k++){
                    
                    int elem = nodesArlequinFine_[connec(k)] -> getNodalElemCorrespondence();
                    double* xsi = nodesArlequinFine_[connec(k)] ->  getNodalXsiCorrespondence();
                    
                    if (rank == 0) 
                        getinterpolatedvelocity_(&u[0],&u[1],&elem,&xsi[0]);
                        
                    MPI_Bcast(&u[0],1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
                    MPI_Bcast(&u[1],1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
               
                    nodesArlequinFine_[connec(k)] -> setVelocity(u);
  
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
template<>
void FSInteraction<2>::transferFluidLoad(){
    
    clearcouplingloads_();
    
    for (int iInterf = 0; iInterf < numInterfaces; iInterf++){
        for (int isolid = 0; isolid < numNodesSolid; isolid++){
            
            int ielem = nodesSolid_[iInterf][isolid] -> getNodalElemCorrespondence();
            double* xsi = nodesSolid_[iInterf][isolid] -> getNodalXsiCorrespondence();

            double load[2] = {};
            elementsFluid_[ielem] -> getBoundaryLoad(xsi,load);

            int inode = isolid+1;
            setcouplingload_(&load[0],&load[1],&inode);
                       
        };//isolid
    };//iInterf
    
};

//------------------------------------------------------------------------------
//---------SETS FLUID AND SOLID MODELS AND GETS ITS BASIC INFORMATIONS----------
//------------------------------------------------------------------------------
template<>
void FSInteraction<2>::transferArlequinLoad(){
    
    clearcouplingloads_();
    
    for (int iInterf = 0; iInterf < numInterfaces; iInterf++){
        for (int isolid = 0; isolid < numNodesSolid; isolid++){
        
            int ielem = nodesSolid_[iInterf][isolid] -> getNodalElemCorrespondence();
            double* xsi = nodesSolid_[iInterf][isolid] -> getNodalXsiCorrespondence();
            
            double load[2] = {};
            if (sqrt(xsi[0]*xsi[0]+xsi[1]*xsi[1]) < 3) elementsArlequinFine_[ielem] -> getBoundaryLoad(xsi,load);
        
            int inode = isolid+1;
            setcouplingload_(&load[0],&load[1],&inode);
                       
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
            double um[2];
            double *x = nodesFluid_[i] -> getCoordinates();
            
            um[0] = nodesFluid_[i] -> getMeshVelocity(0);
            um[1] = nodesFluid_[i] -> getMeshVelocity(1);

            nodesFluid_[i] -> setPreviousCoordinates(0,x[0]);
            nodesFluid_[i] -> setPreviousCoordinates(1,x[1]);

            nodesFluid_[i] -> setPreviousMeshVelocity(0,um[0]);
            nodesFluid_[i] -> setPreviousMeshVelocity(1,um[1]);           
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
            fluidModel.printResults(iTimeStep);
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

    double &alpha_f = fluidModel.fluidParameters.getAlphaF();
    double &alpha_m = fluidModel.fluidParameters.getAlphaM();
    double &gamma = fluidModel.fluidParameters.getGamma();

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

        // if (iTimeStep == 20){
        //     double integ = 1.0;
        //     fluidModel.fluidParameters.setSpectralRadius(integ);  
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
            double accel[2], u[2], uprev[2];
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
            double* x = nodesFluid_[i] -> getCoordinates();
            nodesFluid_[i] -> setPreviousCoordinates(0,x[0]);
            nodesFluid_[i] -> setPreviousCoordinates(1,x[1]);
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
            
            fluidModel.solveFSIFluid(3, 1.e-5, 2);
            
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
                fluidModel.printResults(iTimeStep);
                printstructure_();
            };
        };


    };//Time Steps
};



//------------------------------------------------------------------------------
//----------------SOLVES THE FLUID-STRUCTURE INTERACTION PROBLEM----------------
//------------------------------------------------------------------------------
template<>
void FSInteraction<2>::solveFSIProblemGaussSeidelArlequin(int numTimeSteps){


    std::string om = "omega.txt";
    std::ofstream saidaOmega(om.c_str());


    double sizeSolid = 3 * numNodesSolid;

    ublas::vector<double> X_k(sizeSolid), Y_k(sizeSolid), deltaXi(sizeSolid),
        deltaXii(sizeSolid);
    

    double omega = 1.;
    double mu = 0.;

    if (rank == 0) {
        arlequinModel.printResults(0);
        printstructure_();
    };

    X_k.clear();Y_k.clear();deltaXi.clear();deltaXii.clear();

    double &alpha_f = arlequinModel.fineModel.fluidParameters.getAlphaF();
    double &alpha_m = arlequinModel.fineModel.fluidParameters.getAlphaM();
    double &gamma = arlequinModel.fineModel.fluidParameters.getGamma();

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
        //     arlequinModel.fineModel.fluidParameters.setSpectralRadius(spec);
        //     arlequinModel.coarseModel.fluidParameters.setSpectralRadius(spec);
        //     std::cout << "AQUI " << rank << std::endl;
        // }

        for (int i = 0; i < numNodesArlequinCoarse; i++){
            double accel[2], u[2], uprev[2];
            
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
            double accel[2], u[2], uprev[2], lag[2];
            
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
            double* x = nodesArlequinFine_[i] -> getCoordinates();
            nodesArlequinFine_[i] -> setPreviousCoordinates(0,x[0]);
            nodesArlequinFine_[i] -> setPreviousCoordinates(1,x[1]);
        };


        // std::cout << "AQUI4 " << rank << std::endl;
        // MPI_Barrier(PETSC_COMM_WORLD);

        for (int i = 0; i < numNodesSolid; i++){
            int dof = 3*i+1;
            if (rank == 0) getposition_(&dof,&Y_k(3*i  ));
            dof++;
            if (rank == 0) getposition_(&dof,&Y_k(3*i+1));
            MPI_Bcast(&Y_k(3*i  ),1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
            MPI_Bcast(&Y_k(3*i+1),1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
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

            X_k = Y_k;
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
            arlequinModel.solveFSIArlequin(4, 1.e-3, 2, iTimeStep);

                // std::cout << "AQUI4 " << rank << std::endl;
                // MPI_Barrier(PETSC_COMM_WORLD);
            // }else{
            //     arlequinModel.solveFSIArlequin(5, 1.e-3, 2, iTimeStep);
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
            if (rank == 0) solveframestructure_(&iTimeStep);
            // };
            // std::cout << "AQUI6 " << rank << std::endl;
            // std::cout << "AQUI11 " << rank << std::endl;
            // MPI_Barrier(PETSC_COMM_WORLD);
            
            for (int i = 0; i < numNodesSolid; i++){
                int dof = 3*i+1;
                if (rank == 0) getposition_(&dof,&Y_k(3*i  ));
                dof++;
                if (rank == 0) getposition_(&dof,&Y_k(3*i+1));
                MPI_Bcast(&Y_k(3*i  ),1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
                MPI_Bcast(&Y_k(3*i+1),1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
            };
// std::cout << "AQUI7 " << rank << std::endl;
            // std::cout << "AQUI12 " << rank << std::endl;
            // MPI_Barrier(PETSC_COMM_WORLD);
            
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
                arlequinModel.printResults(iTimeStep);
                printstructure_();
            };
        };


    };//Time Steps
};

#endif
