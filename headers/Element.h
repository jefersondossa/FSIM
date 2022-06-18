
// 
//                   Jeferson W D Fernandes and Rodolfo A K Sanches
//                             University of Sao Paulo
//                           (C) 2017 All Rights Reserved
//
// <LicenseText>
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//-----------------------------------ELEMENT------------------------------------
//------------------------------------------------------------------------------

#ifndef ELEMENT_H
#define ELEMENT_H

#include "Node.h"
#include "BoundaryIntegrationQuadrature.h"
#include "FluidParameters.h"

#include "IntegrationQuadrature.h"
#include "IntegrationQuadrature11.h"
#include "DomainIntegration.h"
// #include "PartitionedQuadrature.hpp"

/// Defines the fluid element object and all the element information

template<int DIM, int DEG>
class Element{
 
public:
    /// Defines the class Node locally
    typedef Node<DIM,DEG>                                       Nodes;

    ///Defines the partitioned integration quadrature rule class locally
    typedef IntegQuadratureSpecial<DIM,DEG>                     SpecialQuad;

    /// Defines the normal integration quadrature rule class locally
    typedef IntegQuadrature<DIM,DEG>                            NormalQuad;

    /// Defines the boundary integration quadrature rule class locally
    typedef BoundaryIntegQuadrature<DIM,DEG>                    BoundaryQuad;

    typedef FluidParameters<DIM,DEG>                            FParameters;

    typedef QuadShapeFunction<DIM,DEG>                          ShapeFunction;

    /// Defines the domain integration class locally
    typedef DomainIntegration<DIM,DEG>                          DIntegration;

private:
    const int     nElNodes = 3*(DIM*DEG-DEG)-2*DIM+4;
    const int     nLocDOF = -8*DIM -21*DEG + 15*DIM*DEG + 16;
    std::vector<Nodes *>   *nodes_;    //Velocity nodes
    FParameters   parameters;
    VecInt        connect_; //Velocity mesh connectivity 
    int           index_;             //Element index
    VecDouble     xK, XK;
    int           sideBoundary_;
    double        meshMovingParameter;
    std::vector<int> neighborElements;

    VecDouble intPointWeightFunction;
    VecDouble intPointWeightFunctionPrev;
    VecDouble intPointWeightFunctionSpecial;
    VecDouble intPointWeightFunctionSpecialPrev;
    VecDouble intPointDistGlueZone;
    VecBool   intPointGlueZone;

    VecInt intPointCorrespElem;
    MatrixDouble intPointCoordinates;
    MatrixDouble intPointCorrespXsi;

    bool          glueZone;
    bool          model; //true for local and false for global
    bool          FSIInterface;
    DIntegration  *DI;   

public:
    /// fluid element constructor
    /// @param int element index @param Connectivity element connectivity
    /// @param vector<Nodes> 
    Element(int index, VecInt &connect, std::vector<Nodes *> &nodes, FParameters &param, DIntegration* &dInt){
        connect_.resize(nElNodes);

        index_ = index;
        for (int i = nElNodes; i--; ) connect_[i] = connect[i];
        nodes_ = &nodes;
        parameters = param;

        // std::cout <<"AASASD 1"<<std::endl;
        DI = dInt;
        // std::cout <<"AASASD 2"<<std::endl;
        glueZone = false;        FSIInterface = false;
        sideBoundary_ = -1;
        neighborElements.clear();

        NormalQuad nQuad = NormalQuad();
        // std::cout <<"AASASD 3"<<std::endl;
        intPointWeightFunction.resize(nQuad.getNumberOfIntegrationPoints());
        intPointWeightFunctionPrev.resize(nQuad.getNumberOfIntegrationPoints());

        intPointWeightFunction.fill(1.);
        intPointWeightFunctionPrev.fill(1.);
        
        SpecialQuad sQuad = SpecialQuad();

        intPointWeightFunctionSpecial.resize(sQuad.getNumberOfIntegrationPoints());
        intPointWeightFunctionSpecialPrev.resize(sQuad.getNumberOfIntegrationPoints());
        intPointDistGlueZone.resize(sQuad.getNumberOfIntegrationPoints());
        intPointGlueZone.resize(sQuad.getNumberOfIntegrationPoints());
        intPointCorrespElem.resize(sQuad.getNumberOfIntegrationPoints());

        intPointCoordinates.resize(sQuad.getNumberOfIntegrationPoints(),2);
        intPointCorrespXsi.resize(sQuad.getNumberOfIntegrationPoints(),2);
        intPointCoordinates.setZero();
        intPointCorrespXsi.setZero();
        intPointDistGlueZone.setZero();
        intPointWeightFunctionSpecial.fill(1.);
        intPointWeightFunctionSpecialPrev.fill(1.);
        intPointGlueZone.fill(false);

        // std::cout <<"AASASD 5"<<std::endl;
        getIntegPointCoordinates();

        // double xsi[2] = {};        
        // double ainv_[2][2] = {};       
        // getJacobianMatrix(xsi, ainv_);
        // std::cout <<"AASASD 6"<<std::endl;

    };

    //........................Element basic information.........................
    /// Clear all element variables
    void clearVariables();

    /// Sets the element connectivity
    /// @param int* element connectivity
    void setConnectivity(VecInt &connect){connect_ = connect;};

    /// Gets the element connectivity
    /// @return element connectivity
    VecInt &getConnectivity(){return connect_;};

    /// Compute and store the spatial jacobian matrix
    /// @param bounded_vector integration point adimensional coordinates
    void getJacobianMatrix(VecDouble &xsi, MatrixDouble &ainv_, double &djac_);

    /// Compute and store the shape function spatial derivatives
    /// @param bounded_vector integration point adimensional coordinates
    void getSpatialDerivatives(VecDouble &xsi, MatrixDouble &ainv_, MatrixDouble &dphi_dx);

    /// Interpolate pressure and its derivatives
    /// @param int integration point index
    /// @param double** shape function spatial derivatives
    /// @param double pressure @param double* pressure derivatives
    void interpolatePressure(int &index, MatrixDouble &dphi_dx, double &p_, VecDouble &dp_dx);

    /// Interpolate mesh velocity
    /// @param int integration point index
    /// @param double* mesh velocity 
    /// @param double* previous time step mesh velocity
    void interpolateMeshVelocity(int &index, VecDouble &umesh_, VecDouble &umeshPrev_);

    /// Interpolate acceleration
    /// @param int integration point index @param double* acceleration
    /// @param double* previous time step acceleration
    void interpolateAcceleration(int &index, VecDouble &a_, VecDouble &aPrev_);

    /// Interpolate velocity
    /// @param int integration point index @param velocity @param previous time step velocity
    void interpolateVelocity(int &index, VecDouble &u_, VecDouble &uPrev_);

    /// Interpolate Coordinates
    /// @param int integration point index @param coordinates @param previous time step coordinates
    void interpolateCoordinates(int &index, VecDouble &x_, VecDouble &xPrev_);

    /// Interpolate Lagrange Multiplier
    /// @param int integration point index @param interpolated lagrange multipliers @param shape functions
    void interpolateLagMultiplier(int &index, VecDouble &lagM_);

    /// Interpolate Lagrange Multiplier
    /// @param int integration point index @param interpolated lagrange multipliers @param shape functions
    void interpolateLagMultiplierDerivatives(MatrixDouble &dphi_dx, MatrixDouble &dL_dx);

    /// Interpolate velocity derivatives
    /// @param double** shape function spatial derivatives @param double** velocity derivatives
    /// @param double** previous time step velocity derivatives
    void interpolateVelDerivatives(MatrixDouble &dphi_dx, MatrixDouble &du_dx, MatrixDouble &duprev_dx);

    /// Compute and store the SUPG, PSPG and LSIC stabilization parameters
    void getParameterSUPG(int &index, double &tSUPG_, double &tPSPG_, double &tLSIC_, MatrixDouble &dphi_dx);
    void getParameterArlequin(int &index, double &tARLQ_, double &tSUPG_, double &tPSPG_, double &tLSIC_, MatrixDouble &dphi_dx);
    void getParameterArlequin2();

    /// Gets the element jacobian determinant
    /// @return element jacobinan determinant
    double getJacobian(){
        VecDouble xsi(DIM);
        MatrixDouble ainv_(DIM,DIM);
        MatrixDouble dphi_dx(nElNodes,DIM);
        
        QuadShapeFunction<DIM,DEG>    shapeQuad;

        xsi[0] = 0.5;
        xsi[1] = 0.5;
        
        // std::cout << "AAA 1 "<< std::endl;

        double djac_ = 0.;
        //Computes the jacobian matrix
        getJacobianMatrix(xsi, ainv_, djac_);
        // std::cout << "AAA 2 "<< std::endl;        
        //Computes spatial derivatives
        // getSpatialDerivatives(xsi, ainv_, dphi_dx);
        // std::cout << "AAA 3 "<< std::endl;

        return djac_;
    };

    /// Compute and store the drag and lift forces at the element boundary
    void computeDragAndLiftForces(double &pressureDragForce, double &pressureLiftForce, double &frictionDragForce,
                                  double &frictionLiftForce, double &dragForce, double &liftForce,
                                  double &pitchingMoment, double & perimeter);

    /// Compute and store the boundary forces
    void getBoundaryLoad(VecDouble &xsi, VecDouble &load);

    /// Gets the spatial jacobian matrix
    /// @param bounded_vector integration point coordinates
    /// @return spatial jacobian matrix
    // void getJacobianMatrixValues(double xsi[], double ainv_[][2]){
    //     getJacobianMatrix(xsi, ainv_);
    //     return;};

    /// Pushs back a term of the inverse incidence, i.e., an element which
    /// contains the node
    /// @param int element
    void pushNeighborElement(int el) {neighborElements.push_back(el);}

    /// Gets the number of elements which contains the node
    /// @return int number of elements which contains the node
    int getNumberOfNeighborElements(){return neighborElements.size();}

    /// Gets an specific member of the inverse incidence
    /// @param int index @return int element of the inverse incidence
    int getNeighborElement(int i){return neighborElements[i];}
    void clearInverseIncidence(){
        neighborElements.clear();
        neighborElements.shrink_to_fit();
    }

    void sortEraseNeighborElements(){
        // std::cout << "AA1 " << index_ << " " << neighborElements.size() << std::endl;
        std::sort(neighborElements.begin(), neighborElements.end());
        neighborElements.erase(std::unique(neighborElements.begin(),neighborElements.end()), neighborElements.end());
        // std::cout << "AA2 " << index_ <<" "<< neighborElements.size() << std::endl;
        // neighborElements.shrink_to_fit();
    }

    /// Sets the element side in boundary
    /// @param int side in boundary
    void setElemSideInBoundary(int side){sideBoundary_ = side;};

    /// Gets the element side in boundary
    /// @return side in boundary
    int getElemSideInBoundary(){return sideBoundary_;};

    /// Sets the mesh moving weighting parameter for solving the Laplace problem
    /// @param double parameter value
    void setMeshMovingParameter(double &value) {meshMovingParameter = value;};

    /// Gets the mesh moving weighting parameter
    /// @return mesh moving weighting parameter
    double getMeshMovingParameter(){return meshMovingParameter;};

    /// Gets the boundary connectivity for boundary integration
    /// @param int* boundary connectivity
    void getBoundaryNodes(int *nodesb_);

    //.................Element intersection and correspondence..................
    /// Gets the element intersection parameters 
    /// @param minimum coordinates @param maximum coordinates 
    /// @param minimum inner product @param maximum inner product
    /// @param side lenght
    void setIntersectionParameters(VecDouble &x, VecDouble &X);

    /// Gets the coordinates intersection parameters
    /// @return minimum and maximum coordinates
    std::pair<VecDouble,VecDouble> getXIntersectionParameter() {return std::make_pair(xK,XK);};

    //.............................Model functions..............................
    /// Sets if the element is in the gluing zone
    /// @param glueZone: if true is in the glue zone
    void setGlueZone(){glueZone = true;}
    bool getGlueZone(){return glueZone;}
    bool getGlueZoneInt(){
        int aux = 0;
        if (glueZone) aux = 1;
        return aux;}

    /// Sets which model the fluid element belongs
    /// @param bool model: true = fine; false = coarse.
    void setModel(bool m){model = m;};

    /// Sets if the element belongs to the fluid structure interface
    void setFSIInterface(){FSIInterface = true;};
    bool getFSIInterface(){return FSIInterface;};

    //......................Integration Points Information......................
    /// Gets the number of integration points of the special quadrature rule
    /// @retunr number of integration point of the special quadrature rule
    int getNumberOfIntegrationPoints(){SpecialQuad sQuad = SpecialQuad(); return sQuad.getNumberOfIntegrationPoints();};

    /// Sets the integration point correspondence to the overlapped mesh
    /// @param int element correspondent @param VecLoc Adimensional coordinates
    void setIntegrationPointCorrespondence(int ipoint, int elem, VecDouble &x){
        intPointCorrespElem[ipoint] = elem;
        intPointCorrespXsi(ipoint,0) = x[0];
        intPointCorrespXsi(ipoint,1) = x[1]; };

    /// Gets the integration point correspondence - element
    /// @return overlapped element correspondence
    int getIntegPointCorrespondenceElement(int index){return intPointCorrespElem[index];};

    /// Compute and store the integration points global coordinates
    void getIntegPointCoordinates();

    /// Gets the integration point global coordinates
    /// @param int integration point index @return integration point coordinates
    VecDouble getIntegPointCoordinatesValue(int index){
        VecDouble aux(2);
        aux[0] = intPointCoordinates(index,0); 
        aux[1] = intPointCoordinates(index,1); 
        return aux;
    };

    /// Sets the integration point energy weight function
    /// @param int integration point index 
    /// @param double energy weight function value
    void setIntegPointWeightFunction();

    /// Gets the integration point energy weight function
    /// @param int integration point index @return energy weight function value
    double getIntegPointWeightFunction(int index)
    {return intPointWeightFunction[index];};

    /// Sets if the integration point is the gluing zone
    /// @param int integration point index 
    void setIntegPointInGlueZone(int index){intPointGlueZone[index] = true;};

    /// Gets true if the integration point is in the gluing zone
    /// @param int integration point index @return If is in the gluing zone
    bool getIntegPointInGlueZone(int index){return intPointGlueZone[index];};

    /// Sets the integration point signaled distance function
    /// @param int integration point index
    /// @param signaled distance function valye
    void setIntegPointDistFunction(int index, double val)
    {intPointDistGlueZone[index] = val;};

    /// Gets the integration point signaled distance function value
    /// @return integration point signaled distance function value
    double getIntegPointDistFunction(int index)
    {return intPointDistGlueZone[index];};

    //.......................Element vectors and matrices.......................
    /// Compute and store the element matrix for the incompressible flow problem
    /// @param int integration point index
    void getElemMatrix(int &index, MatrixDouble &dphi_dx, double &tSUPG_, double &tPSPG_, double &tLSIC_, double &weight_, double &djac_, MatrixDouble &jacobianNRMatrix);

    /// Compute and store the element matrix for the Laplace/Poisson problem
    void getElemLaplMatrix(double &weight_, double &djac_, MatrixDouble &dphi_dx, double** jacobianNRMatrix);

    /// Sets the boundary conditions for the incompressible flow problem
    void setBoundaryConditions(MatrixDouble &jacobianNRMatrix, VecDouble &rhsVector);

    /// Sets the boundary conditions for the Laplace/Poisson problem
    void setBoundaryConditionsLaplace(double** jacobianNRMatrix, double* rhsVector);

    ///Compute and store the residual vector for the incompressible flow problem
    /// @param int integration point index
    void getResidualVector(int &index, MatrixDouble &dphi_dx, double &tSUPG_, double &tPSPG_, double &tLSIC_, double &weight_, double &djac_, VecDouble &rhsVector);

    /// Compute and store the residual vector for the Laplace/Poisson problem
    void getResidualVectorLaplace(double *rhsVector);

    /// Apply the boundary conditions and returns the matrix and residual vector
    /// for the Lagrange multiplier operator matrix, used when computing the 
    /// operator term from different meshes
    /// @param LocalMatrix Lagrange multiplier operator matrix
    /// @return residual vector and Lagrange multiplier operator matrix with
    /// boundary conditions applied
    void getRhsVectorAndBoundaryConditions(double** Ajac, double* rhsVector){
        setBoundaryConditionsLagrangeMultipliers(Ajac, rhsVector);       
        return;
    };

    /// Sets the boundary conditions to the Lagrange multiplier operator
    void setBoundaryConditionsLagrangeMultipliers(double** jacobianNRMatrix, double* rhsVector);

    /// Compute and store the Lagrange multiplier operator when integrating 
    /// the same mesh portion
    void getLagrangeMultipliersSameMesh(double **lagrMultMatrix, double *lagrMultVector, double *rhsVector);
    void getLagrangeMultipliersSUPG_PSPG_SameMesh(double **jacobianNRMatrix, double *rhsVector);
    void getLagrangeMultipliersArlequinSameMesh(double **arlequinStab, double **laplMatrix, double *arlequinStabVector);

    /// Compute and store the Lagrange multiplier operator when integrationg
    /// the different mesh portion
    /// @param int element of the coarse mesh (used to verify which integration
    /// point belongs to the coarse mesh element)
    void getLagrangeMultipliersDifferentMesh(int &ielem, double &tPSPG2_,double* press, double* velx, double* vely,
                                             double* velxPrev, double* velyPrev, double **lagrMultMatrix, double *rhsVectorLM, double *rhsVector);
    void getLagrangeMultipliersSUPG_PSPG_DifferentMesh(int &ielem, double &tPSPG2_,double* press, double* velx, double* vely,
                                                       double **jacobianNRMatrix, double *rhsVector);
    void getLagrangeMultipliersArlequinDifferentMesh(int &ielem, double &tPSPG2_,double* press, double* velx, double* vely,
                                                     double **arlequinStab, double **laplMatrix, double *arlequinStabVector);

    //...............................Problem type...............................
    /// Compute the Transient Navier-Stokes problem matrices and vectors
    void getTransientNavierStokes(MatrixDouble &jacobianNRMatrix, VecDouble &rhsVector);

    /// Compute the Steady Laplace problem matrices and vectors 
    /// (usually for the mesh moving step)
    void getSteadyLaplace(double** jacobianNRMatrix, double* rhsVector);
    void getSteadyLaplace2(double** jacobianNRMatrix, double* rhsVector);

};


#endif

