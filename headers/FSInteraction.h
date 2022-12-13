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
#include "hdf5.h"
#include "Arlequin.h"

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

template<int DIM, int DEG>
class FSInteraction{
public:
    /// Defines locally the class Fluid
    typedef Fluid<DIM,DEG>                  FluidModel;
    /// Defines locally the class Element
    typedef typename FluidModel::Elements   Elements;
    /// Defines locally the class Node
    typedef typename FluidModel::Node       Nodes;
    /// Defines locally the class Boundary
    typedef typename FluidModel::Boundaries Boundary;

    /// Defines locally the class Arlequin
    typedef Arlequin<DIM,DEG>               ArlequinModel;
   
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

    int nElNodes = 3*(DIM*DEG-DEG)-2*DIM+4;
    int nLocDOF = -8*DIM -21*DEG + 15*DIM*DEG + 16;
    int nBdNodes = 3*(1-DEG)+DIM*(2*DEG-1);

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
    /// with Aitken relaxation fluid-structure interaction problem with the
    /// Arlequin method 
    /// @param int number of time steps
    void solveFSIProblemArlequin(int numTimeSteps);


};


#endif
