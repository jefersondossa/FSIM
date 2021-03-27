
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

#include "Node.hpp"
#include "BoundaryIntegrationQuadrature.hpp"
#include "FluidParameters.hpp"

#include "IntegrationQuadrature.hpp"
#include "IntegrationQuadrature11.hpp"
#include "DomainIntegration.hpp"
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
    int           connect_[3*(DIM*DEG-DEG)-2*DIM+4]; //Velocity mesh connectivity 
    int           index_;             //Element index
    double        xK[DIM], XK[DIM];
    int           sideBoundary_;
    double        meshMovingParameter;

    double* intPointWeightFunction;
    double* intPointWeightFunctionPrev;
    double* intPointWeightFunctionSpecial;
    double* intPointWeightFunctionSpecialPrev;
    double* intPointDistGlueZone;
    bool*   intPointGlueZone;

    int* intPointCorrespElem;
    double** intPointCoordinates;
    double** intPointCorrespXsi;

    bool          glueZone;
    bool          model; //true for local and false for global
    bool          FSIInterface;
    DIntegration  *DI;   

public:
    /// fluid element constructor
    /// @param int element index @param Connectivity element connectivity
    /// @param vector<Nodes> 
    Element(int index, int *connect, std::vector<Nodes *> &nodes, FParameters &param, DIntegration* &dInt){
        index_ = index;
        for (int i = nElNodes; i--; ) connect_[i] = connect[i];
        nodes_ = &nodes;
        parameters = param;

        // std::cout <<"AASASD 1"<<std::endl;
        DI = dInt;
        // std::cout <<"AASASD 2"<<std::endl;
        glueZone = false;        FSIInterface = false;
        sideBoundary_ = -1;

        NormalQuad nQuad = NormalQuad();
        // std::cout <<"AASASD 3"<<std::endl;
        intPointWeightFunction = new double[nQuad.getNumberOfIntegrationPoints()];
        intPointWeightFunctionPrev = new double[nQuad.getNumberOfIntegrationPoints()];

        for (int i=0; i < (nQuad.end() - nQuad.begin()); i++){
            intPointWeightFunction[i] = 1.;
            intPointWeightFunctionPrev[i] = 1.;
        };
// std::cout <<"AASASD 4"<<std::endl;
        SpecialQuad sQuad = SpecialQuad();

        intPointWeightFunctionSpecial = new double[sQuad.getNumberOfIntegrationPoints()]{};
        intPointWeightFunctionSpecialPrev = new double[sQuad.getNumberOfIntegrationPoints()]{};
        intPointDistGlueZone = new double[sQuad.getNumberOfIntegrationPoints()]{};
        intPointGlueZone = new bool[sQuad.getNumberOfIntegrationPoints()]{};
        intPointCorrespElem = new int[sQuad.getNumberOfIntegrationPoints()]{};

        intPointCoordinates = new double*[sQuad.getNumberOfIntegrationPoints()]{};
        intPointCorrespXsi = new double*[sQuad.getNumberOfIntegrationPoints()]{};
        for(int i = 0; i< sQuad.getNumberOfIntegrationPoints(); i++) {
            intPointCoordinates[i] = new double[2]();
            intPointCorrespXsi[i] = new double[2]();
        }
        for (int i = 0; i < (sQuad.end() - sQuad.begin()); i++){
            intPointWeightFunctionSpecial[i] = 1.;
            intPointWeightFunctionSpecialPrev[i] = 1.;
            intPointGlueZone[i] = false; 
        };
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
    void setConnectivity(int *connect){for (int i=nElNodes; i--; ) connect_[i] = connect[i];};

    /// Gets the element connectivity
    /// @return element connectivity
    int* getConnectivity(){return connect_;};

    /// Compute and store the spatial jacobian matrix
    /// @param bounded_vector integration point adimensional coordinates
    void getJacobianMatrix(double *xsi, double **ainv_, double &djac_);

    /// Compute and store the shape function spatial derivatives
    /// @param bounded_vector integration point adimensional coordinates
    void getSpatialDerivatives(double *xsi, double **ainv_, double **dphi_dx);

    /// Interpolate pressure and its derivatives
    /// @param int integration point index
    /// @param double** shape function spatial derivatives
    /// @param double pressure @param double* pressure derivatives
    void interpolatePressure(int &index, double **dphi_dx, double &p_, double *dp_dx);

    /// Interpolate mesh velocity
    /// @param int integration point index
    /// @param double* mesh velocity 
    /// @param double* previous time step mesh velocity
    void interpolateMeshVelocity(int &index, double *umesh_, double *umeshPrev_);

    /// Interpolate acceleration
    /// @param int integration point index @param double* acceleration
    /// @param double* previous time step acceleration
    void interpolateAcceleration(int &index, double *a_, double *aPrev_);

    void computeInverseAndDeterminant(double &djac_, double **ainv_, double dx_dxsi[DIM][DIM]);

    /// Interpolate velocity
    /// @param int integration point index @param velocity @param previous time step velocity
    void interpolateVelocity(int &index, double *u_, double *uPrev_);

    /// Interpolate Coordinates
    /// @param int integration point index @param coordinates @param previous time step coordinates
    void interpolateCoordinates(int &index, double *x_, double *xPrev_);

    /// Interpolate Lagrange Multiplier
    /// @param int integration point index @param interpolated lagrange multipliers @param shape functions
    void interpolateLagMultiplier(int &index, double *lagM_);

    /// Interpolate Lagrange Multiplier
    /// @param int integration point index @param interpolated lagrange multipliers @param shape functions
    void interpolateLagMultiplierDerivatives(double **dphi_dx, double dL_dx[DIM][DIM]);

    /// Interpolate velocity derivatives
    /// @param double** shape function spatial derivatives @param double** velocity derivatives
    /// @param double** previous time step velocity derivatives
    void interpolateVelDerivatives(double **dphi_dx, double du_dx[DIM][DIM], double duprev_dx[DIM][DIM]);

    /// Compute and store the SUPG, PSPG and LSIC stabilization parameters
    void getParameterSUPG(int &index, double &tSUPG_, double &tPSPG_, double &tLSIC_, double **dphi_dx);
    void getParameterArlequin(int &index, double &tARLQ_, double &tSUPG_, double &tPSPG_, double &tLSIC_, double **dphi_dx);
    void getParameterArlequin2();

    /// Gets the element jacobian determinant
    /// @return element jacobinan determinant
    double getJacobian(){
        double xsi[DIM] = {};

        double **ainv_;
        ainv_ = new double*[DIM];
        for (int i = DIM; i--; ) ainv_[i] = new double[DIM];
        
        double **dphi_dx;
        dphi_dx = new double*[DIM];
        for (int i = DIM; i--; ) dphi_dx[i] = new double[nElNodes];

        QuadShapeFunction<DIM,DEG>    shapeQuad;

        xsi[0] = 0.5;
        xsi[1] = 0.5;
                        
        double djac_ = 0.;
        //Computes the jacobian matrix
        getJacobianMatrix(xsi, ainv_, djac_);
        
        //Computes spatial derivatives
        getSpatialDerivatives(xsi, ainv_, dphi_dx);

        for (int i = DIM; i--; ) delete [] dphi_dx[i];
        delete [] dphi_dx;
        for (int i = DIM; i--; ) delete [] ainv_[i];
        delete [] ainv_;

        return djac_;};

    /// Compute and store the drag and lift forces at the element boundary
    void computeDragAndLiftForces(double &pressureDragForce, double &pressureLiftForce, double &frictionDragForce,
                                  double &frictionLiftForce, double &dragForce, double &liftForce,
                                  double &pitchingMoment, double & perimeter);

    /// Compute and store the boundary forces
    void getBoundaryLoad(double* xsi, double* load);

    /// Gets the spatial jacobian matrix
    /// @param bounded_vector integration point coordinates
    /// @return spatial jacobian matrix
    // void getJacobianMatrixValues(double xsi[], double ainv_[][2]){
    //     getJacobianMatrix(xsi, ainv_);
    //     return;};

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
    void setIntersectionParameters(double x[], double X[]);

    /// Gets the coordinates intersection parameters
    /// @return minimum and maximum coordinates
    std::pair<double*,double*> getXIntersectionParameter() {return std::make_pair(xK,XK);};

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

    //......................Integration Points Information......................
    /// Gets the number of integration points of the special quadrature rule
    /// @retunr number of integration point of the special quadrature rule
    int getNumberOfIntegrationPoints(){SpecialQuad sQuad = SpecialQuad(); return (sQuad.end() - sQuad.begin());};

    /// Sets the integration point correspondence to the overlapped mesh
    /// @param int element correspondent @param VecLoc Adimensional coordinates
    void setIntegrationPointCorrespondence(int ipoint, int elem, double* x){
        intPointCorrespElem[ipoint] = elem;
        intPointCorrespXsi[ipoint][0] = x[0];
        intPointCorrespXsi[ipoint][1] = x[1]; };

    /// Gets the integration point correspondence - element
    /// @return overlapped element correspondence
    int getIntegPointCorrespondenceElement(int index){return intPointCorrespElem[index];};

    /// Compute and store the integration points global coordinates
    void getIntegPointCoordinates();

    /// Gets the integration point global coordinates
    /// @param int integration point index @return integration point coordinates
    double* getIntegPointCoordinatesValue(int index){return intPointCoordinates[index];};

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
    void getElemMatrix(int &index, double **dphi_dx, double &tSUPG_, double &tPSPG_, double &tLSIC_, double &weight_, double &djac_, double **jacobianNRMatrix);

    /// Compute and store the element matrix for the Laplace/Poisson problem
    void getElemLaplMatrix(double &weight_, double &djac_, double **dphi_dx, double** jacobianNRMatrix);

    /// Sets the boundary conditions for the incompressible flow problem
    void setBoundaryConditions(double **jacobianNRMatrix, double *rhsVector);

    /// Sets the boundary conditions for the Laplace/Poisson problem
    void setBoundaryConditionsLaplace(double** jacobianNRMatrix, double* rhsVector);

    ///Compute and store the residual vector for the incompressible flow problem
    /// @param int integration point index
    void getResidualVector(int &index, double **dphi_dx, double &tSUPG_, double &tPSPG_, double &tLSIC_, double &weight_, double &djac_, double *rhsVector);

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
    void getTransientNavierStokes(double **jacobianNRMatrix, double *rhsVector);

    /// Compute the Steady Laplace problem matrices and vectors 
    /// (usually for the mesh moving step)
    void getSteadyLaplace(double** jacobianNRMatrix, double* rhsVector);
    void getSteadyLaplace2(double** jacobianNRMatrix, double* rhsVector);

};


//------------------------------------------------------------------------------
//--------------------------------IMPLEMENTATION--------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//-----------------INVERSE AND DETERMINANT OF A DIMxDIM MATRIX------------------
//------------------------------------------------------------------------------
template<>
void Element<2,1>::computeInverseAndDeterminant(double &djac_, double **ainv_, double dx_dxsi[2][2]) {

    djac_ = dx_dxsi[0][0] * dx_dxsi[1][1] - dx_dxsi[0][1] * dx_dxsi[1][0];

    ainv_[0][0] =  dx_dxsi[1][1] / djac_;
    ainv_[0][1] = -dx_dxsi[1][0] / djac_;
    ainv_[1][0] = -dx_dxsi[0][1] / djac_;
    ainv_[1][1] =  dx_dxsi[0][0] / djac_;

    return;
}

template<>
void Element<2,2>::computeInverseAndDeterminant(double &djac_, double **ainv_, double dx_dxsi[2][2]) {

    djac_ = dx_dxsi[0][0] * dx_dxsi[1][1] - dx_dxsi[0][1] * dx_dxsi[1][0];

    ainv_[0][0] =  dx_dxsi[1][1] / djac_;
    ainv_[0][1] = -dx_dxsi[1][0] / djac_;
    ainv_[1][0] = -dx_dxsi[0][1] / djac_;
    ainv_[1][1] =  dx_dxsi[0][0] / djac_;

    return;
}

template<>
void Element<3,1>::computeInverseAndDeterminant(double &djac_, double **ainv_, double dx_dxsi[3][3]) {

    djac_ = dx_dxsi[0][0] * dx_dxsi[1][1] * dx_dxsi[2][2] - dx_dxsi[0][0] * dx_dxsi[1][2] * dx_dxsi[2][1]
          + dx_dxsi[0][1] * dx_dxsi[1][2] * dx_dxsi[2][0] - dx_dxsi[0][1] * dx_dxsi[1][0] * dx_dxsi[2][2]
          + dx_dxsi[0][2] * dx_dxsi[1][0] * dx_dxsi[2][1] - dx_dxsi[0][2] * dx_dxsi[1][1] * dx_dxsi[2][0];

    ainv_[0][0] = (dx_dxsi[1][1] * dx_dxsi[2][2] - dx_dxsi[1][2] * dx_dxsi[2][1]) / djac_;
    ainv_[1][0] = (dx_dxsi[0][2] * dx_dxsi[2][1] - dx_dxsi[0][1] * dx_dxsi[2][2]) / djac_;
    ainv_[2][0] = (dx_dxsi[0][1] * dx_dxsi[1][2] - dx_dxsi[0][2] * dx_dxsi[1][1]) / djac_;
    ainv_[0][1] = (dx_dxsi[1][2] * dx_dxsi[2][0] - dx_dxsi[1][0] * dx_dxsi[2][2]) / djac_;
    ainv_[1][1] = (dx_dxsi[0][0] * dx_dxsi[2][2] - dx_dxsi[0][2] * dx_dxsi[2][0]) / djac_;
    ainv_[2][1] = (dx_dxsi[0][2] * dx_dxsi[1][0] - dx_dxsi[0][0] * dx_dxsi[1][2]) / djac_;
    ainv_[0][2] = (dx_dxsi[1][0] * dx_dxsi[2][1] - dx_dxsi[1][1] * dx_dxsi[2][0]) / djac_;
    ainv_[1][2] = (dx_dxsi[0][1] * dx_dxsi[2][0] - dx_dxsi[0][0] * dx_dxsi[2][1]) / djac_;
    ainv_[2][2] = (dx_dxsi[0][0] * dx_dxsi[1][1] - dx_dxsi[0][1] * dx_dxsi[1][0]) / djac_;        

    return;
}

template<>
void Element<3,2>::computeInverseAndDeterminant(double &djac_, double **ainv_, double dx_dxsi[3][3]) {

    djac_ = dx_dxsi[0][0] * dx_dxsi[1][1] * dx_dxsi[2][2] - dx_dxsi[0][0] * dx_dxsi[1][2] * dx_dxsi[2][1]
          + dx_dxsi[0][1] * dx_dxsi[1][2] * dx_dxsi[2][0] - dx_dxsi[0][1] * dx_dxsi[1][0] * dx_dxsi[2][2]
          + dx_dxsi[0][2] * dx_dxsi[1][0] * dx_dxsi[2][1] - dx_dxsi[0][2] * dx_dxsi[1][1] * dx_dxsi[2][0];

    ainv_[0][0] = (dx_dxsi[1][1] * dx_dxsi[2][2] - dx_dxsi[1][2] * dx_dxsi[2][1]) / djac_;
    ainv_[1][0] = (dx_dxsi[0][2] * dx_dxsi[2][1] - dx_dxsi[0][1] * dx_dxsi[2][2]) / djac_;
    ainv_[2][0] = (dx_dxsi[0][1] * dx_dxsi[1][2] - dx_dxsi[0][2] * dx_dxsi[1][1]) / djac_;
    ainv_[0][1] = (dx_dxsi[1][2] * dx_dxsi[2][0] - dx_dxsi[1][0] * dx_dxsi[2][2]) / djac_;
    ainv_[1][1] = (dx_dxsi[0][0] * dx_dxsi[2][2] - dx_dxsi[0][2] * dx_dxsi[2][0]) / djac_;
    ainv_[2][1] = (dx_dxsi[0][2] * dx_dxsi[1][0] - dx_dxsi[0][0] * dx_dxsi[1][2]) / djac_;
    ainv_[0][2] = (dx_dxsi[1][0] * dx_dxsi[2][1] - dx_dxsi[1][1] * dx_dxsi[2][0]) / djac_;
    ainv_[1][2] = (dx_dxsi[0][1] * dx_dxsi[2][0] - dx_dxsi[0][0] * dx_dxsi[2][1]) / djac_;
    ainv_[2][2] = (dx_dxsi[0][0] * dx_dxsi[1][1] - dx_dxsi[0][1] * dx_dxsi[1][0]) / djac_;        

    return;
}

//------------------------------------------------------------------------------
//---------------CREATES AN AUXILIARY FINITE ELEMENT OF DIMENSION---------------
//----------------------DIM-1 FOR THE BOUNDARY INTEGRATION----------------------
//------------------------------------------------------------------------------
template<>
void Element<2,1>::getBoundaryNodes(int *nodesb_){

    if(sideBoundary_ == 0){
        nodesb_[0] = connect_[2]; 
        nodesb_[1] = connect_[1]; 
    }else{
        if(sideBoundary_ == 1){
            nodesb_[0] = connect_[0]; 
            nodesb_[1] = connect_[2]; 
        }else{
            nodesb_[0] = connect_[1];
            nodesb_[1] = connect_[0];
        };        
    };

    return;
}
template<>
void Element<2,2>::getBoundaryNodes(int *nodesb_){
    //!!!!!!!ATENÇÃO ESSAS CONECTIVIDADES FORAM ALTERADAS EM RELAÇÃO AO PROGRAMA Fluid.
    if(sideBoundary_ == 0){
        nodesb_[0] = connect_[1]; 
        nodesb_[1] = connect_[2]; 
        nodesb_[2] = connect_[4];         
    }else{
        if(sideBoundary_ == 1){
            nodesb_[0] = connect_[2]; 
            nodesb_[1] = connect_[0]; 
            nodesb_[2] = connect_[5]; 
        }else{
            nodesb_[0] = connect_[0];
            nodesb_[1] = connect_[1];
            nodesb_[2] = connect_[3];
        };        
    };

    return;
}

template<>
void Element<3,1>::getBoundaryNodes(int *nodesb_){

    if(sideBoundary_ == 0){
        nodesb_[0] = connect_[1]; 
        nodesb_[1] = connect_[2]; 
        nodesb_[2] = connect_[3];         
    }else{
        if(sideBoundary_ == 1){
            nodesb_[0] = connect_[0]; 
            nodesb_[1] = connect_[3]; 
            nodesb_[2] = connect_[2]; 
        }else{
            if(sideBoundary_ == 2){
                nodesb_[0] = connect_[0];
                nodesb_[1] = connect_[1];
                nodesb_[2] = connect_[3];
            }else{
                nodesb_[0] = connect_[0];
                nodesb_[1] = connect_[2];
                nodesb_[2] = connect_[1];
            }
        };        
    };

    return;
}

template<>
void Element<3,2>::getBoundaryNodes(int *nodesb_){

    if(sideBoundary_ == 0){
        nodesb_[0] = connect_[2]; 
        nodesb_[1] = connect_[3]; 
        nodesb_[2] = connect_[1];
        nodesb_[3] = connect_[9]; 
        nodesb_[4] = connect_[8]; 
        nodesb_[5] = connect_[5];
    }else{
        if(sideBoundary_ == 1){
            nodesb_[0] = connect_[0]; 
            nodesb_[1] = connect_[3]; 
            nodesb_[2] = connect_[2];
            nodesb_[3] = connect_[7]; 
            nodesb_[4] = connect_[9]; 
            nodesb_[5] = connect_[6];
        }else{
            if(sideBoundary_ == 2){
                nodesb_[0] = connect_[3]; 
                nodesb_[1] = connect_[0]; 
                nodesb_[2] = connect_[1];
                nodesb_[3] = connect_[7]; 
                nodesb_[4] = connect_[4]; 
                nodesb_[5] = connect_[8];
            }else{
                nodesb_[0] = connect_[0]; 
                nodesb_[1] = connect_[2]; 
                nodesb_[2] = connect_[1];
                nodesb_[3] = connect_[6]; 
                nodesb_[4] = connect_[5]; 
                nodesb_[5] = connect_[4];
            }
        };        
    };

    return;
}


//------------------------------------------------------------------------------
//----------------------SET ELEMENT INTERSECTION PARAMETERS---------------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Element<DIM,DEG>::setIntersectionParameters(double *x, double *X) {
    xK[0] = x[0]; xK[1] = x[1]; 
    XK[0] = X[0]; XK[1] = X[1];
    return;
};

//------------------------------------------------------------------------------
//----------------------SET ELEMENT INTERSECTION PARAMETERS---------------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Element<DIM,DEG>::setIntegPointWeightFunction() {
    
    double xsi[DIM] = {};
    ShapeFunction shapeQuad;
    double phi_[nElNodes] = {};
    
    NormalQuad nQuad = NormalQuad();
    for(int i = 0; i < nQuad.getNumberOfIntegrationPoints(); i++) {
        intPointWeightFunctionPrev[i] = intPointWeightFunction[i];
        intPointWeightFunction[i] = 0.;
    }

    for(int i = 0; i < nQuad.getNumberOfIntegrationPoints(); i++) {
        intPointWeightFunctionSpecialPrev[i] = intPointWeightFunctionSpecial[i];
        intPointWeightFunctionSpecial[i] = 0.;
    }

    int index=0;

    for(double* it = nQuad.begin(); it != nQuad.end(); it++){
        
       xsi[0] = nQuad.PointList(index,0);
       xsi[1] = nQuad.PointList(index,1);
            
       //Computes the velocity shape functions
       shapeQuad.evaluate(xsi,phi_);

       for (int j=0; j<nElNodes; j++){
           intPointWeightFunction[index] += phi_[j] * (*nodes_)[connect_[j]] -> getWeightFunction();
       };
       // intPointWeightFunction(index) = 1.;
       index++;
    }; 

    // if (index_ == 2117){
    //     for (int i = 0; i < intPointWeightFunction.size(); ++i)
    //     {
    //         std::cout << "Int point " << i << " " << intPointWeightFunction(i) << " " << intPointWeightFunctionPrev(i) << std::endl;      
    //     }
    // } 

    index = 0;
    SpecialQuad sQuad = SpecialQuad();

    for(double* it = sQuad.begin(); it != sQuad.end(); it++){
        
       xsi[0] = sQuad.PointList(index,0);
       xsi[1] = sQuad.PointList(index,1);
            
       //Computes the velocity shape functions
       shapeQuad.evaluate(xsi,phi_);

       for (int j=0; j<nElNodes; j++){
           intPointWeightFunctionSpecial[index] += phi_[j] * (*nodes_)[connect_[j]] -> getWeightFunction();
       };
       // intPointWeightFunction(index) = 1.;
       index++;
    }; 

    //     typename Nodes::VecLocD x;
    // int cont = 0;
    // if (model == false){
    //     for (int i = 0; i < 6; i++){
    //         x = (*nodes_)[connect_[i]] -> getCoordinates();
    //         if((x(0) < 1.501) || (x(0) > 6.499)){            
    //             cont++;
    //         };  
    //     };
    //     if (cont == 6){
    //         for (int i = 0; i < intPointWeightFunction.size(); i++){
    //             intPointWeightFunction(i) = 1.;
    //         };
    //         std::cout << "aqui " << cont << " " << index_ << std::endl;
    //     }
    // };

    // if (norm_2(intPointWeightFunction - intPointWeightFunctionPrev) > 1.e-8) std::cout << "AQUI " << index_ << std::endl
    return;
};

//------------------------------------------------------------------------------
//------------------COMPUTES THE INTEGRATION POINT COORDINATE-------------------
//------------------------------------------------------------------------------
template<int DIM,int DEG>
void Element<DIM,DEG>::getIntegPointCoordinates(){

    SpecialQuad sQuad = SpecialQuad();
    double xsi[DIM] = {};
    ShapeFunction   shapeQuad;
    double phi_[nElNodes] = {};

    for (int i = 0; i < sQuad.end() - sQuad.begin(); i++){
        double x[DIM] = {};

        for (int k = DIM; k--; ) xsi[k] = sQuad.PointList(i,k);

        shapeQuad.evaluate(xsi,phi_);

        for (int k = DIM; k--; ) intPointCoordinates[i][k] = 0.;

        for (int j = 0; j < nElNodes; j++)
            for (int k = DIM; k--; )
                intPointCoordinates[i][k] += (*nodes_)[connect_[j]] -> getCoordinateValue(k) * phi_[j];
        
    };

    return;
};

//------------------------------------------------------------------------------
//---------------------------CLEAR ELEMENT VARIABLES----------------------------
//------------------------------------------------------------------------------
template<>
void Element<2,2>::clearVariables(){
        
    glueZone = false;

    SpecialQuad sQuad = SpecialQuad();

    for (int i=0; i < (sQuad.end() - sQuad.begin()); i++){
        intPointWeightFunction[i] = 1.;
        intPointGlueZone[i] = false;
        intPointCorrespElem[i] = 0;
        for (int j =0; j<2; j++) {
            intPointCoordinates[i][j] = 0.;
            intPointCorrespXsi[i][j] = 0.;
        }
    };
    
    getIntegPointCoordinates();

    return;
}; 


//------------------------------------------------------------------------------
//-------------------------SPATIAL TRANSFORM - JACOBIAN-------------------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Element<DIM,DEG>::getJacobianMatrix(double *xsi, double **ainv_, double &djac_) {

    //Computes the spatial Jacobian matrix and its inverse
    double **dphi;
    dphi = new double*[nElNodes];
    for (int i = nElNodes; i--; ) dphi[i] = new double[DIM];
    double dx_dxsi[DIM][DIM] = {};
    double xna_[DIM] = {};

    ShapeFunction shapeQuad;
    shapeQuad.evaluateGradient(xsi,dphi);

    double &alpha_f = parameters.getAlphaF();

    for (int i = nElNodes; i--; ){
        for (int j = DIM; j--; ){
            // Approximate the integration space
            xna_[j] = alpha_f * (*nodes_)[connect_[i]] -> getCoordinateValue(j) + 
                      (1. - alpha_f) * (*nodes_)[connect_[i]] -> getPreviousCoordinateValue(j);

            for (int k = DIM; k--; ){
                dx_dxsi[j][k] += xna_[j] * dphi[i][k];
            };
        };
    };

    //Computing the jacobian determinant and Inverse
    computeInverseAndDeterminant(djac_, ainv_, dx_dxsi);

    for (int i = nElNodes; i--; ) delete [] dphi[i];
    delete [] dphi;

    return;
};

//------------------------------------------------------------------------------
//-----------------------------SPATIAL DERIVATIVES------------------------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Element<DIM,DEG>::getSpatialDerivatives(double *xsi, double **ainv_, double **dphi_dx) {
    
    // typename QuadShapeFunction<2,2>::ValueDDeriv ddphi;
    
    double **dphi;
    dphi = new double*[nElNodes];
    for (int i = nElNodes; i--; ) dphi[i] = new double[DIM];

    ShapeFunction shapeQuad;
    
    shapeQuad.evaluateGradient(xsi,dphi);
    // shapeQuad.evaluateHessian(xsi,ddphi);
    
    for (int i = nElNodes; i--; )
        for (int j = DIM; j--; )
            dphi_dx[i][j] = 0.;

    //Quadratic shape functions spatial first derivatives
    // noalias(dphi_dx) = prod(ainv_,dphi);
    for (int i = nElNodes; i--; )
        for (int k = DIM; k--; )
            for (int j = DIM; j--; )
                dphi_dx[i][k] += ainv_[k][j] * dphi[i][j];

    for (int i = nElNodes; i--; ) delete [] dphi[i];
    delete [] dphi;

    return;
};

//------------------------------------------------------------------------------
//------------------INTERPOLATES PRESSURE AND ITS DERIVATIVES-------------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Element<DIM,DEG>::interpolatePressure(int &index, double **dphi_dx, double &p_, double *dp_dx) {
    p_ = 0.;
    for (int i = DIM; i--; ) dp_dx[i] = 0.;

    for (int i = nElNodes; i--; ){
        double shapeFi = DI -> phi_[i][index];

        p_ += (*nodes_)[connect_[i]] -> getPressure() * shapeFi;
        
        for (int j = DIM; j--; ) dp_dx[j] += (*nodes_)[connect_[i]] -> getPressure() * dphi_dx[i][j];
    }
    return;
}


//------------------------------------------------------------------------------
//---------------INTERPOLATES MESH VELOCITY AND ITS DERIVATIVES-----------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Element<DIM,DEG>::interpolateMeshVelocity(int &index, double *umesh_, double *umeshPrev_) {

    for (int i = DIM; i--; ) {
        umesh_[i] = 0.;
        umeshPrev_[i] = 0.;
    }

    for (int i = nElNodes; i--; ){
        double shapeFi = DI -> phi_[i][index];
        for (int j = DIM; j--; ){
            umesh_[j] += (*nodes_)[connect_[i]] -> getMeshVelocity(j) * shapeFi;
            umeshPrev_[j] += (*nodes_)[connect_[i]] -> getPreviousMeshVelocity(j) * shapeFi;
        }
    }

    return;
}

//------------------------------------------------------------------------------
//--------------------------INTERPOLATES ACCELERATION---------------------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Element<DIM,DEG>::interpolateAcceleration(int &index, double *a_, double *aPrev_) {

    for (int i = DIM; i--; ) {
        a_[i] = 0.;
        aPrev_[i] = 0.;
    }

    for (int i = nElNodes; i--; ){
        double shapeFi = DI -> phi_[i][index];
        for (int j = DIM; j--; ){
            a_[j] += (*nodes_)[connect_[i]] -> getAcceleration(j) * shapeFi;
            aPrev_[j] += (*nodes_)[connect_[i]] -> getPreviousAcceleration(j) * shapeFi;
        }
    }

    return;
}

//------------------------------------------------------------------------------
//----------------------------INTERPOLATES VELOCITY-----------------------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Element<DIM,DEG>::interpolateCoordinates(int &index, double *x_, double *xPrev_) {

    for (int i = DIM; i--; ) {
        x_[i] = 0.;
        xPrev_[i] = 0.;
    }
    for (int i = nElNodes; i--; ){
        double shapeFi = DI -> phi_[i][index];
        for (int j = DIM; j--; ){
            x_[j] += (*nodes_)[connect_[i]] -> getCoordinateValue(j) * shapeFi;
            xPrev_[j] += (*nodes_)[connect_[i]] -> getPreviousCoordinateValue(j) * shapeFi;
        }
    }
    return;
}


//------------------------------------------------------------------------------
//----------------------------INTERPOLATES VELOCITY-----------------------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Element<DIM,DEG>::interpolateVelocity(int &index, double *u_, double *uPrev_) {

    for (int i = DIM; i--; ) {
        u_[i] = 0.;
        uPrev_[i] = 0.;
    }
    for (int i = nElNodes; i--; ){
        double shapeFi = DI -> phi_[i][index];
        for (int j = DIM; j--; ){
            u_[j] += (*nodes_)[connect_[i]] -> getVelocity(j) * shapeFi;
            uPrev_[j] += (*nodes_)[connect_[i]] -> getPreviousVelocity(j) * shapeFi;
        }
    }
    return;
}

//------------------------------------------------------------------------------
//----------------------------INTERPOLATES VELOCITY-----------------------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Element<DIM,DEG>::interpolateLagMultiplier(int &index, double *lagM_) {

    for (int i = DIM; i--; )
        lagM_[i] = 0.;

    for (int i = nElNodes; i--; )
        for (int j = DIM; j--; )
            lagM_[j] += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(j) * DI -> phi_[i][index];

    return;
}

//------------------------------------------------------------------------------
//----------------------------INTERPOLATES VELOCITY-----------------------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Element<DIM,DEG>::interpolateVelDerivatives(double **dphi_dx, double du_dx[DIM][DIM], double duprev_dx[DIM][DIM]) {

    for (int i = DIM; i--; ) {
        for (int j = DIM; j--; ) {
            du_dx[i][j] = 0.;
            duprev_dx[i][j] = 0.;
        }
    }

    for (int i = nElNodes; i--; ){
        for (int j = DIM; j--; ){
            for (int k = DIM; k--; ){
                du_dx[k][j] += (*nodes_)[connect_[i]] -> getVelocity(k) * dphi_dx[i][j];
                duprev_dx[k][j] += (*nodes_)[connect_[i]] -> getPreviousVelocity(k) * dphi_dx[i][j];
            }
        }
    }

    return;
}

//------------------------------------------------------------------------------
//----------------------------INTERPOLATES VELOCITY-----------------------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Element<DIM,DEG>::interpolateLagMultiplierDerivatives(double **dphi_dx, double dL_dx[DIM][DIM]) {

    for (int i = DIM; i--; )
        for (int j = DIM; j--; )
            dL_dx[i][j] = 0.;
            

    for (int i = nElNodes; i--; )
        for (int j = DIM; j--; )
            for (int k = DIM; k--; )
                dL_dx[k][j] += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(k) * dphi_dx[i][j];
            
    return;
}


//------------------------------------------------------------------------------
//-------------INTERPOLATES VELOCITY, PRESSURE AND ITS DERIVATIVES--------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Element<DIM,DEG>::getBoundaryLoad(double* xsi, double* load) {

    double &visc_ = parameters.getViscosity();
    double &alpha_f = parameters.getAlphaF();

    double phi_[nElNodes] = {};

    double **dphi_dx;
    dphi_dx = new double*[nElNodes];
    for (int i = nElNodes; i--; ) dphi_dx[i] = new double[DIM];

    double shearStress[DIM][DIM] = {};
    
    double **dphi;
    dphi = new double*[nElNodes];
    for (int i = nElNodes; i--; ) dphi[i] = new double[DIM];
    double **ainv_;
    ainv_ = new double*[DIM];
    for (int i = DIM; i--; ) ainv_[i] = new double[DIM];

    double t_vector[2] = {};
    double n_vector[2] = {};
    double ident[2][2] = {}; ident[0][0] = 1.; ident[1][1] = 1.;
    ShapeFunction shapeQuad;
    

    //Computes the shape functions        
    shapeQuad.evaluate(xsi,phi_);
    double djac_ = 0.;
    //Computes the jacobian matrix
    getJacobianMatrix(xsi, ainv_, djac_);
    
    //Computes spatial derivatives
    getSpatialDerivatives(xsi, ainv_, dphi_dx);

    shapeQuad.evaluateGradient(xsi, dphi); 

    //Velocity Derivatives
    double du_dx[DIM][DIM], duprev_dx[DIM][DIM], duna_dx[DIM][DIM];
    interpolateVelDerivatives(dphi_dx, du_dx, duprev_dx);
    for (int i=DIM; i--; )
        for (int j = DIM; j--; )
            duna_dx[i][j] = alpha_f * du_dx[i][j] + (1. - alpha_f) * duprev_dx[i][j];

    //Pressure
    double p_;
    double dp_dx[DIM] = {};
    int index = 0;
    interpolatePressure(index, dphi_dx, p_, dp_dx);
   
    shearStress[0][0] = 2. * visc_ * duna_dx[0][0];
    shearStress[0][1] = visc_ * (duna_dx[0][1] + duna_dx[1][0]);
    shearStress[1][0] = visc_ * (duna_dx[0][1] + duna_dx[1][0]);
    shearStress[1][1] = 2. * visc_ * duna_dx[1][1];

    int nodesb_[3];
    int aux[3];

    if(sideBoundary_ == 0){
        nodesb_[0] = connect_[2]; 
        nodesb_[1] = connect_[4]; 
        nodesb_[2] = connect_[1];
        aux[0] = 2;
        aux[1] = 4;
        aux[2] = 1;   
    }else{
        if(sideBoundary_ == 1){
            nodesb_[0] = connect_[0]; 
            nodesb_[1] = connect_[5]; 
            nodesb_[2] = connect_[2]; 
            aux[0] = 0;
            aux[1] = 5;
            aux[2] = 2;
        }else{
            nodesb_[0] = connect_[1];
            nodesb_[1] = connect_[3];
            nodesb_[2] = connect_[0];
            aux[0] = 1;
            aux[1] = 3;
            aux[2] = 0;
        };        
    };
    double dx_dxsiB[3][1] = {};
    double xna_[3] = {};

    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 2; j++){
            // Approximate the integration space
            xna_[j] = (*nodes_)[nodesb_[i]] -> getCoordinateValue(j);

            dx_dxsiB[j][0] += xna_[j] * dphi[aux[i]][j];    
        };
    };
    double Maux[1][1] = {};
    for (int k = 0; k < 3; k++)
        Maux[0][0] += dx_dxsiB[k][0] * dx_dxsiB[k][0];
    double djacb_ = 0.;
    djacb_ = std::sqrt(Maux[0][0]);
    n_vector[0] = -dx_dxsiB[1][0] / djacb_;
    n_vector[1] =  dx_dxsiB[0][0] / djacb_;

    // std::cout << "CONNEc " << nodesb_[0] << " " << nodesb_[1] << " " << nodesb_[2] << std::endl;
    // std::cout << "N vector " << xna_[0] << " " << xna_[1] << " " << n_vector[0] << " " << n_vector[1] << std::endl;

    // if (sideBoundary_ == 1){
    //     for (int i = 0; i < 6; i++){
    //         t_vector[0] -= dphi[1][i] * (*nodes_)[connect_[i]] -> getCoordinateValue(0);
    //         t_vector[1] -= dphi[1][i] * (*nodes_)[connect_[i]] -> getCoordinateValue(1);
    //     };        
    // };

    // if (sideBoundary_ == 2){
    //     for (int i = 0; i < 6; i++){
    //         t_vector[0] += dphi[0][i] * (*nodes_)[connect_[i]] -> getCoordinateValue(0);
    //         t_vector[1] += dphi[0][i] * (*nodes_)[connect_[i]] -> getCoordinateValue(1);
    //     };        
    // };

    // if (sideBoundary_ == 0){
    //     std::cout << "VERIFICAR VETOR NORMAL - getBoundaryLoad" << std::endl;
    // };

    // n_vector[0] =  t_vector[1] / std::sqrt(t_vector[0]*t_vector[0] + t_vector[1]*t_vector[1]);
    // n_vector[1] = -t_vector[0] / std::sqrt(t_vector[0]*t_vector[0] + t_vector[1]*t_vector[1]);

    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            load[i] += -p_ * ident[i][j] * n_vector[j] + shearStress[i][j] * n_vector[j];
    // load = -p_ * prod(ident,n_vector) + prod(shearStress,n_vector);

    //std::cout << "N Vector " << sideBoundary_ << " " <<  n_vector(0) << " " << n_vector(1) << " " << load(0) << " " << load(1) << " " << p_ << std::endl;
    

    for (int i = nElNodes; i--; ) delete [] dphi_dx[i];
    delete [] dphi_dx;
    for (int i = nElNodes; i--; ) delete [] dphi[i];
    delete [] dphi;
    for (int i = DIM; i--; ) delete [] ainv_[i];
    delete [] ainv_;

    return;
};

//------------------------------------------------------------------------------
//-------------INTERPOLATES VELOCITY, PRESSURE AND ITS DERIVATIVES--------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Element<DIM,DEG>::computeDragAndLiftForces(double &pressureDragForce, double &pressureLiftForce, double &frictionDragForce,
                                          double &frictionLiftForce, double &dragForce, double &liftForce,
                                          double &pitchingMoment, double & perimeter) {
    
    double localNodesBoundary_[3][2] = {};
    BoundShapeFunction<2,2> shapeBound;//Boundary shape function 

    int nodesb_[3]; 
    if(sideBoundary_ == 0){
        nodesb_[0] = connect_[1]; 
        nodesb_[1] = connect_[4]; 
        nodesb_[2] = connect_[2]; 
        for (int i=0; i<2; i++){
            localNodesBoundary_[0][i] = (*nodes_)[connect_[1]] -> getCoordinateValue(i);
            localNodesBoundary_[1][i] = (*nodes_)[connect_[4]] -> getCoordinateValue(i);
            localNodesBoundary_[2][i] = (*nodes_)[connect_[2]] -> getCoordinateValue(i);
        };
    }else{
        if(sideBoundary_ == 1){
            nodesb_[0] = connect_[2]; 
            nodesb_[1] = connect_[5]; 
            nodesb_[2] = connect_[0]; 
            for (int i=0; i<2; i++){
                localNodesBoundary_[0][i] = (*nodes_)[connect_[2]] -> getCoordinateValue(i);
                localNodesBoundary_[1][i] = (*nodes_)[connect_[5]] -> getCoordinateValue(i);
                localNodesBoundary_[2][i] = (*nodes_)[connect_[0]] -> getCoordinateValue(i);
            };
        }else{
            nodesb_[0] = connect_[0];
            nodesb_[1] = connect_[3];
            nodesb_[2] = connect_[1];
            for (int i=0; i<2; i++){
                localNodesBoundary_[0][i] = (*nodes_)[connect_[0]] -> getCoordinateValue(i);
                localNodesBoundary_[1][i] = (*nodes_)[connect_[3]] -> getCoordinateValue(i);
                localNodesBoundary_[2][i] = (*nodes_)[connect_[1]] -> getCoordinateValue(i);
            };
        };        
    };

    BoundaryQuad           bQuad;     //Boundary Integration Quadrature
    ShapeFunction          shapeQuad;
    double phi_[nElNodes] = {};
    
    double **dphi_dx;
    dphi_dx = new double*[nElNodes];
    for (int i = nElNodes; i--; ) dphi_dx[i] = new double[DIM];



    std::pair<double*,double*> gaussQuad;
    double n_vector[2] = {};
    double shearStress[2][2] = {};
    double ident[2][2] = {}; ident[0][0] = 1.; ident[1][1] = 1.;
    double load_friction[2] = {};
    double load_pressure[2] = {};

    double **ainv_;
    ainv_ = new double*[DIM];
    for (int i = DIM; i--; ) ainv_[i] = new double[DIM];
    
    double xsi[2] = {};

    double &visc_ = parameters.getViscosity();

    gaussQuad = bQuad.GaussQuadrature();
    double moment = 0.;
    double per = 0.;

    int index = 0;
    for(double* it = bQuad.begin(); it != bQuad.end(); it++){
        
        double xsiB = gaussQuad.first[index];
        double weightB = gaussQuad.second[index];

        if(sideBoundary_ == 2){
            xsi[0] = (-xsiB + 1.) / 2.;
            xsi[1] = 0.;
        };
        if(sideBoundary_ == 1){
            xsi[1] = (xsiB + 1.) / 2.;
            xsi[0] = 0.;
        };
        if(sideBoundary_ == 0){
            xsi[0] = (xsiB + 1.) / 2.;
            xsi[1] = 1. - xsi[0];
        };

        //Computes the velocity shape functions
        shapeQuad.evaluate(xsi,phi_);
        double djac_ = 0.;
        //Computes the jacobian matrix
        getJacobianMatrix(xsi, ainv_, djac_);

        //Computes spatial derivatives
        getSpatialDerivatives(xsi, ainv_, dphi_dx);

        //Velocity Derivatives
        double du_dx[DIM][DIM], duprev_dx[DIM][DIM], duna_dx[DIM][DIM];
        interpolateVelDerivatives(dphi_dx, du_dx, duprev_dx);

        //Pressure
        double p_;
        double dp_dx[DIM] = {};
        interpolatePressure(index, dphi_dx, p_, dp_dx);

        double x_[DIM] = {}; double xPrev_[DIM] = {};
        interpolateCoordinates(index, x_, xPrev_);


        double* phib_ = shapeBound.getShapeFunction(gaussQuad.first[index]);
        double* dphib_ = shapeBound.getShapeFunctionDerivative(gaussQuad.first[index]);

        double Tx=0.; double Ty = 0.;
        
        for (int i=0; i<3; i++){
            Tx += localNodesBoundary_[i][0] * dphib_[i];
            Ty += localNodesBoundary_[i][1] * dphib_[i];
        };

        double jacb_ = std::sqrt(Tx*Tx + Ty*Ty);
        
        n_vector[0] =  Ty / jacb_;
        n_vector[1] = -Tx / jacb_;

        shearStress[0][0] = 2. * visc_ * du_dx[0][0];
        shearStress[0][1] = visc_ * (du_dx[0][1] + du_dx[1][0]);
        shearStress[1][0] = visc_ * (du_dx[0][1] + du_dx[1][0]);
        shearStress[1][1] = 2. * visc_ * du_dx[1][1];

        for (int i = 0; i < DIM; i++){
            for (int j = 0; j < DIM; j++){
                load_pressure[i] += -p_ * ident[i][j] * n_vector[j] * jacb_ * weightB;
                load_friction[i] += shearStress[i][j] * n_vector[j] * jacb_ * weightB;
            }
        }

        moment += ((-p_ + shearStress[0][0] + shearStress[1][0]) * x_[1] 
                 -(-p_ + shearStress[0][1] + shearStress[1][1]) * (x_[0] - 0.248792267683901))
                 * jacb_ * weightB;
        per += jacb_ * weightB;
        //std::cout << "n_vector " << load_pressure(0) << " " << load_pressure(1) << std::endl;
        index++;
    };

    perimeter = per;
    pitchingMoment = moment;

    pressureDragForce = -load_pressure[0];
    pressureLiftForce = -load_pressure[1];
    
    frictionDragForce = -load_friction[0];
    frictionLiftForce = -load_friction[1];

    dragForce = pressureDragForce + frictionDragForce;
    liftForce = pressureLiftForce + frictionLiftForce;

    for (int i = nElNodes; i--; ) delete [] dphi_dx[i];
    delete [] dphi_dx;
    for (int i = DIM; i--; ) delete [] ainv_[i];
    delete [] ainv_;

    return;
};




//------------------------------------------------------------------------------
//------------------COMPUTES THE SUPG STABILIZATION PARAMETER-------------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Element<DIM,DEG>::getParameterSUPG(int &index, double &tSUPG_, double &tPSPG_, double &tLSIC_, double **dphi_dx) {

    double tSUGN1_ = 0.;
    double tSUGN2_ = 0.;
    double tSUGN3_ = 0.;
    double hRGN_ = 0.;
    double r[DIM] = {};
    double s[DIM] = {};
    double hUGN_ = 0.;

    tSUPG_ = 0.;  

    double u__[DIM] = {};
    double aux = 0.;
    double aux2 = 0.;

    double &alpha_f = parameters.getAlphaF();
    double &visc_ = parameters.getViscosity();
    double &dens_ = parameters.getDensity();
    double &dTime_ = parameters.getTimeStep();

    for (int i = nElNodes; i--; ){
        double a1 = 0.;
        for (int j = DIM; j--; ){
            double ua = alpha_f * (*nodes_)[connect_[i]] -> getVelocity(j) + (1. - alpha_f) * (*nodes_)[connect_[i]] -> getPreviousVelocity(j);
            double uma = alpha_f * (*nodes_)[connect_[i]] -> getMeshVelocity(j) + (1. - alpha_f) * (*nodes_)[connect_[i]] -> getPreviousMeshVelocity(j);

            ua -= uma;
            u__[j] += ua * DI -> phi_[i][index];
            a1 += ua*ua;
        }

        // a1 = std::sqrt(a1);
        for (int j = DIM; j--; ) r[j] += std::sqrt(a1) * dphi_dx[i][j];
    };


    for (int j = DIM; j--; ) aux += u__[j]*u__[j];
    double uNorm = std::sqrt(aux);

    // for (int i = nElNodes; i--; ){
    //     for (int j = DIM; j--; ){
    //         r[j] += uNorm * dphi_dx[j][i];
    //     }
    // }

    for (int j = DIM; j--; ) aux2 += r[j]*r[j];
    double rNorm = std::sqrt(aux2);

    if(uNorm > 1.e-10){
        for (int j = DIM; j--; ) s[j] = u__[j] / uNorm;
    }else{
        for (int j = DIM; j--; ) s[j] = 1. / std::sqrt(2.);
    };

    if (rNorm >= 1.e-10){
        for (int j = DIM; j--; ) r[j] /= rNorm;
    }else{
        for (int j = DIM; j--; ) r[j] = 1. / std::sqrt(2.);
    };
    
    for (int i = nElNodes; i--; ){
        for (int j = DIM; j--; ){
            hRGN_ += r[j] * dphi_dx[i][j];
            hUGN_ += s[j] * dphi_dx[i][j];
        }
    };
    hRGN_ = std::fabs(hRGN_);
    hUGN_ = std::fabs(hUGN_);

    if (hRGN_ >= 1.e-10){
        hRGN_ = 2. / hRGN_;
    }else{
        hRGN_ = 2. / 1.e-10;
    };

    if (hUGN_ >= 1.e-10){
        hUGN_ = 2. / hUGN_;
    }else{
        hUGN_ = 2. / 1.e-10;
    };    

    if (uNorm >= 1.e-10){
        tSUGN1_ = hUGN_ / (2. * uNorm);
    }else{
        tSUGN1_ = hUGN_ / 2.e-10;
    };
              
    tSUGN2_ = dTime_ / 2.;

    tSUGN3_ = hRGN_ * hRGN_ / (4. * visc_ / dens_);
   
    if (std::fabs(tSUGN1_) <= 1.e-10) tSUGN1_ = 1.e-10;
    if (std::fabs(tSUGN3_) <= 1.e-10) tSUGN3_ = 1.e-10;

    //if(model == false) std::cout << "SUPG " << tSUGN1_ << " " << tSUGN3_ << std::endl;
    //Computing tSUPG parameter
    tSUPG_ = 1. / std::sqrt(1. / (tSUGN1_ * tSUGN1_) + 
                       1. / (tSUGN2_ * tSUGN2_) + 
                       1. / (tSUGN3_ * tSUGN3_));
    //tSUPG_ = 0.;


    tPSPG_ = 1*tSUPG_;

    tLSIC_ = tSUPG_ * uNorm * uNorm;


    return;
};

//------------------------------------------------------------------------------
//------------------COMPUTES THE SUPG STABILIZATION PARAMETER-------------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Element<DIM,DEG>::getParameterArlequin(int &index, double &tARLQ_, double &tSUPG_, double &tPSPG_, double &tLSIC_, double **dphi_dx) {

    double        tSUGN1_;
    double        tSUGN2_;
    double        tSUGN3_;
    double        hRGN_;
    
    double   r[2] = {};
    double   s[2] = {};
    double   sl[2] = {};
    double   rl[2] = {};

    tSUPG_ = 0.;
    tSUGN1_ = 0.;
    tSUGN2_ = 0.;
    tSUGN3_ = 0.;
    hRGN_ = 0.;
    double hUGN_ = 0.;
    
    double u__ = 0.;
    double v__ = 0.;
    double lx__ = 0.;
    double ly__ = 0.;

    double &alpha_f = parameters.getAlphaF();
    double &visc_ = parameters.getViscosity();
    double &dens_ = parameters.getDensity();
    double &dTime_ = parameters.getTimeStep();
    double &k1 = parameters.getArlequinK1();

    for (int i = 0; i < nElNodes; i++){
        double ua = alpha_f * (*nodes_)[connect_[i]] -> getVelocity(0) + (1. - alpha_f) * (*nodes_)[connect_[i]] -> getPreviousVelocity(0);
        double va = alpha_f * (*nodes_)[connect_[i]] -> getVelocity(1) + (1. - alpha_f) * (*nodes_)[connect_[i]] -> getPreviousVelocity(1);
        // double ua = (*nodes_)[connect_[i]] -> getVelocity(0);
        // double va = (*nodes_)[connect_[i]] -> getVelocity(1);

        double uma = alpha_f * (*nodes_)[connect_[i]] -> getMeshVelocity(0) + (1. - alpha_f) * (*nodes_)[connect_[i]] -> getPreviousMeshVelocity(0);
        double vma = alpha_f * (*nodes_)[connect_[i]] -> getMeshVelocity(1) + (1. - alpha_f) * (*nodes_)[connect_[i]] -> getPreviousMeshVelocity(1);
            
        double lxa = (*nodes_)[connect_[i]] -> getLagrangeMultiplier(0);
        double lya = (*nodes_)[connect_[i]] -> getLagrangeMultiplier(1);
            
        ua -= uma;
        va -= vma;
        
        u__ += ua * DI -> phi_[i][index];
        v__ += va * DI -> phi_[i][index];

        lx__ += lxa * DI -> phi_[i][index];
        ly__ += lya * DI -> phi_[i][index];
    };

    double uNorm = std::sqrt(u__ * u__ + v__ * v__);
    double lNorm = std::sqrt(lx__ * lx__ + ly__ * ly__);

    if(uNorm > 1.e-10){
        s[0] = u__ / uNorm;
        s[1] = v__ / uNorm;
    }else{
        s[0] = 1. / std::sqrt(2.);
        s[1] = 1. / std::sqrt(2.);
    };
    if(lNorm > 1.e-10){
        sl[0] = lx__ / lNorm;
        sl[1] = ly__ / lNorm;
    }else{
        sl[0] = 1. / std::sqrt(2.);
        sl[1] = 1. / std::sqrt(2.);
    };


    for (int i = 0; i < nElNodes; i++){
        double ua = alpha_f * (*nodes_)[connect_[i]] -> getVelocity(0) + (1. - alpha_f) * (*nodes_)[connect_[i]] -> getPreviousVelocity(0);
        double va = alpha_f * (*nodes_)[connect_[i]] -> getVelocity(1) + (1. - alpha_f) * (*nodes_)[connect_[i]] -> getPreviousVelocity(1);
        // double ua = (*nodes_)[connect_[i]] -> getVelocity(0);
        // double va = (*nodes_)[connect_[i]] -> getVelocity(1);

        double uma = alpha_f * (*nodes_)[connect_[i]] -> getMeshVelocity(0) + (1. - alpha_f) * (*nodes_)[connect_[i]] -> getPreviousMeshVelocity(0);
        double vma = alpha_f * (*nodes_)[connect_[i]] -> getMeshVelocity(1) + (1. - alpha_f) * (*nodes_)[connect_[i]] -> getPreviousMeshVelocity(1);

        double lxa = (*nodes_)[connect_[i]] -> getLagrangeMultiplier(0);
        double lya = (*nodes_)[connect_[i]] -> getLagrangeMultiplier(1);

        ua -= uma;
        va -= vma;
        
        r[0] += std::sqrt(ua * ua + va * va) * dphi_dx[i][0];
        r[1] += std::sqrt(ua * ua + va * va) * dphi_dx[i][1];

        rl[0] += std::sqrt(lxa * lxa + lya * lya) * dphi_dx[i][0];
        rl[1] += std::sqrt(lxa * lxa + lya * lya) * dphi_dx[i][1];
    };

    double rNorm = std::sqrt(r[0]*r[0] + r[1]*r[1]);
    double rlNorm = std::sqrt(rl[0]*rl[0] + rl[1]*rl[1]);

    if (rNorm >= 1.e-10){
        r[0] /= rNorm;
        r[1] /= rNorm;
    }else{
        r[0] = 1. / std::sqrt(2.);
        r[1] = 1. / std::sqrt(2.);
    };

    if (rlNorm >= 1.e-10){
        rl[0] /= rlNorm;
        rl[1] /= rlNorm;
    }else{
        rl[0] = 1. / std::sqrt(2.);
        rl[1] = 1. / std::sqrt(2.);
    };
    
    double hUGNL_ = 0.;
    double hRGNL_ = 0.;
    for (int i = 0; i < nElNodes; i++){
        hRGN_ += std::fabs(r[0] * dphi_dx[i][0] + r[1] * dphi_dx[i][1]);
        hUGN_ += std::fabs(s[0] * dphi_dx[i][0] + s[1] * dphi_dx[i][1]);        

        hRGNL_ += std::fabs(rl[0] * dphi_dx[i][0] + rl[1] * dphi_dx[i][1]);
        hUGNL_ += std::fabs(sl[0] * dphi_dx[i][0] + sl[1] * dphi_dx[i][1]);        
    };

    if (hRGN_ >= 1.e-10){
        hRGN_ = 2. / hRGN_;
    }else{
        hRGN_ = 2. / 1.e-10;
    };
    if (hRGNL_ >= 1.e-10){
        hRGNL_ = 2. / hRGNL_;
    }else{
        hRGNL_ = 2. / 1.e-10;
    };


    if (hUGN_ >= 1.e-10){
        hUGN_ = 2. / hUGN_;
    }else{
        hUGN_ = 2. / 1.e-10;
    };

    if (hUGNL_ >= 1.e-10){
        hUGNL_ = 2. / hUGNL_;
    }else{
        hUGNL_ = 2. / 1.e-10;
    };    

    if (uNorm >= 1.e-10){
        tSUGN1_ = hUGN_ / (2. * uNorm);
    }else{
        tSUGN1_ = hUGN_ / 2.e-10;
    };

    double tSUGN1L_ = 0.;
    if (lNorm >= 1.e-10){
        tSUGN1L_ = hUGNL_ / (2. * lNorm);
    }else{
        tSUGN1L_ = hUGNL_ / 2.e-10;
    };    
              
    tSUGN2_ = dTime_ / 2.;

    tSUGN3_ = hRGN_ * hRGN_ / (4. * visc_ / dens_);

    double tSUGN3L_ = 0.;
    tSUGN3L_ = hRGNL_ * hRGNL_ / (4. * visc_ / dens_);

   
    if (std::fabs(tSUGN1_) <= 1.e-10) tSUGN1_ = 1.e-10;
    if (std::fabs(tSUGN3_) <= 1.e-10) tSUGN3_ = 1.e-10;

    if (std::fabs(tSUGN1L_) <= 1.e-10) tSUGN1L_ = 1.e-10;
    if (std::fabs(tSUGN3L_) <= 1.e-10) tSUGN3L_ = 1.e-10;

    //if(model == false) std::cout << "SUPG " << tSUGN1_ << " " << tSUGN3_ << std::endl;
    //Computing tSUPG parameter
    tSUPG_ = 1. / std::sqrt(1. / (tSUGN1_ * tSUGN1_) + 
                       1. / (tSUGN2_ * tSUGN2_) + 
                       1. / (tSUGN3_ * tSUGN3_));
    //tSUPG_ = 0.;

    tARLQ_ = -1. / std::sqrt(1. / (tSUGN1L_ * tSUGN1L_) + 
                       1. / (tSUGN2_ * tSUGN2_) + 
                       1. / (tSUGN3L_ * tSUGN3L_))*0;


    //if (tSUPG_ > 10) tSUPG_ = 0.;

    tPSPG_ = 1*tSUPG_;
    // if (std::fabs(lagMx_) > 0){
    //  //   std::cout << "aqe" << std::endl;
    //     tARLQ_ = djac_ * std::sqrt(u_ * u_ + v_ * v_) / std::sqrt(lagMx_ * lagMx_ + lagMy_ * lagMy_);//-tSUPG_*1;
    // }else{
    tARLQ_ = -1. * k1 * tSUPG_ * 1.e-2;
    //}
    //tARLQ_ = 0.;

    //tARLQ_ = 0.;
    //tSUPG_ = 0.;
    tLSIC_ = tSUPG_ * uNorm * uNorm;
 

     //tARLQ_ = 0.;
    // tSUPG_ = 0.;



    // LocalMatrix lambda, inercia, conveccao, visccc;
    // lambda.clear(); inercia.clear(); conveccao.clear(); visccc.clear();
    // double le = 0.;
    // double li = 0.;
    // double lc = 0.;
    // double lv = 0.;

    // for (int i = 0; i < 6; i++){
    //     for (int j = 0; j < 6; j++){
            
    //         // lambda(2*i  ,2*j  ) += phi_[i] * phi_[j] * weight_ * djac_;
    //         // lambda(2*i+1,2*j+1) += phi_[i] * phi_[j] * weight_ * djac_;

    //         // conveccao(2*i  ,2*j  ) += dphi_dx[0][i] * (du_dx*dphi_dx[0][j] + dv_dx*dphi_dx[1][j]) * weight_ * djac_;
    //         // conveccao(2*i  ,2*j+1) += dphi_dx[1][i] * (du_dx*dphi_dx[0][j] + dv_dx*dphi_dx[1][j]) * weight_ * djac_;
    //         // conveccao(2*i+1,2*j  ) += dphi_dx[0][i] * (du_dy*dphi_dx[0][j] + dv_dy*dphi_dx[1][j]) * weight_ * djac_;
    //         // conveccao(2*i+1,2*j+1) += dphi_dx[1][i] * (du_dy*dphi_dx[0][j] + dv_dy*dphi_dx[1][j]) * weight_ * djac_;

    //         // inercia(2*i  ,2*j  ) += dphi_dx[0][i] * dphi_dx[0][j] * weight_ * djac_;
    //         // inercia(2*i  ,2*j+1) += dphi_dx[1][i] * dphi_dx[1][j] * weight_ * djac_;
    //         // inercia(2*i+1,2*j  ) += dphi_dx[0][i] * dphi_dx[0][j] * weight_ * djac_;
    //         // inercia(2*i+1,2*j+1) += dphi_dx[1][i] * dphi_dx[1][j] * weight_ * djac_;

    //         // visccc(2*i  ,2*j  ) += ((ddphi_dx(0,0)(i)+ddphi_dx(0,1)(i)) * (2. * ddphi_dx(0,0)(j) + ddphi_dx(1,1)(j))) * weight_ * djac_;
    //         // visccc(2*i  ,2*j+1) += ((ddphi_dx(0,0)(i)+ddphi_dx(0,1)(i)) * (ddphi_dx(0,1)(j))) * weight_ * djac_;
    //         // visccc(2*i+1,2*j  ) += ((ddphi_dx(1,1)(i)+ddphi_dx(0,1)(i)) * (ddphi_dx(0,1)(j))) * weight_ * djac_;
    //         // visccc(2*i+1,2*j+1) += ((ddphi_dx(1,1)(i)+ddphi_dx(0,1)(i)) * (2. * ddphi_dx(1,1)(j) + ddphi_dx(0,0)(j))) * weight_ * djac_;
    //     }
    // }

    // LocalVector U_, DU_, UG_, DUG_;
    // U_.clear(); DU_.clear(); UG_.clear();

    // for (int i=0; i<6; i++){
    //     U_(2*i  ) = (*nodes_)[connect_[i]] -> getVelocity(0);
    //     U_(2*i+1) = (*nodes_)[connect_[i]] -> getVelocity(1);
    //     UG_(2*i  ) = (*nodes_)[connect_[i]] -> getVelocityGlobal(0);
    //     UG_(2*i+1) = (*nodes_)[connect_[i]] -> getVelocityGlobal(1);
    //     DU_(2*i  ) = (*nodes_)[connect_[i]] -> getAcceleration(0);
    //     DU_(2*i+1) = (*nodes_)[connect_[i]] -> getAcceleration(1);
    //     DUG_(2*i  ) = (*nodes_)[connect_[i]] -> getAccelerationGlobal(0);
    //     DUG_(2*i+1) = (*nodes_)[connect_[i]] -> getAccelerationGlobal(1);
    // };

    // le = norm_2(prod(lambda,U_-UG_));
    // lc = norm_2(prod(conveccao,U_));
    // li = norm_2(prod(inercia,DU_));
    // lv = norm_2(prod(visccc,U_));

    // //std::cout << "asd " << index_ << " " << le << " " << lc << " " << li << " " << tARLQ_ << std::endl;
    //  // std::cout << "antes " << index_ << " " << tARLQ_ << std::endl;

    // double tp1 = std::fabs(le) / std::fabs(lc);
    // double tp2 = std::fabs(le) / std::fabs(li);
    // double tp3 = tp1 / visc_;
    // tp3 = le / (visc_*lv);

    // if (std::fabs(tp1) <= 1.e-10) tp1 = 1.e-10;
    // if (std::fabs(tp2) <= 1.e-10) tp2 = 1.e-10;
    // if (std::fabs(tp3) <= 1.e-10) tp3 = 1.e-10;

    // tARLQ_ = -0.0 * k1 / std::sqrt(1. / (tp1*tp1) +
    //                         1. / (tp2*tp2) +
    //                         1. / (tp3*tp3));

    // std::cout << "depois " << index_ << " " << le << " " << lv << " " << tARLQ_ << std::endl;



    return;
};

// //------------------------------------------------------------------------------
// //------------------COMPUTES THE SUPG STABILIZATION PARAMETER-------------------
// //------------------------------------------------------------------------------
// template<int DIM, int DEG>
// void Element<DIM,DEG>::getParameterArlequin2() {

//     double xsi[2] = {};
//     double phi_[6] = {};
    
//     double **dphi_dx;
//     dphi_dx = new double*[DIM];
//     for (int i = DIM; i--; ) dphi_dx[i] = new double[nElNodes];


//     ShapeFunction           shapeQuad;
//     double **ainv_;
//     ainv_ = new double*[DIM];
//     for (int i = DIM; i--; ) ainv_[i] = new double[DIM];
//     NormalQuad nQuad = NormalQuad();

//     double &visc_ = parameters.getViscosity();
//     double &dens_ = parameters.getDensity();
//     double &dTime_ = parameters.getTimeStep();
//     double &k1 = parameters.getArlequinK1();

//     double lambda[18][18] = {};
//     double inercia[18][18] = {};
//     double conveccao[18][18] = {};
//     double visccc[18][18] = {};

//     double le = 0.;
//     double li = 0.;
//     double lc = 0.;
//     double lv = 0.;

//     tARLQ_ = 0.;
    
//     int index = 0;

//     for(double* it = nQuad.begin(); it != nQuad.end(); it++){

//         //Defines the integration points adimentional coordinates
//         xsi[0] = nQuad.PointList(index,0);
//         xsi[1] = nQuad.PointList(index,1);

//         //Computes the velocity shape functions
//         shapeQuad.evaluate(xsi,phi_);

//         //Returns the quadrature integration weight
//         double weight_ = nQuad.WeightList(index);

//         //Computes the jacobian matrix
//         getJacobianMatrix(xsi, ainv_);

//         getSpatialDerivatives(xsi, ainv_, dphi_dx);
        
//         getVelAndDerivatives(phi_, dphi_dx);

//         for (int i = 0; i < 6; i++){
//             for (int j = 0; j < 6; j++){
                
//                 lambda[2*i  ][2*j  ] += phi_[i] * phi_[j] * weight_ * djac_;
//                 lambda[2*i+1][2*j+1] += phi_[i] * phi_[j] * weight_ * djac_;

//                 conveccao[2*i  ][2*j  ] += dphi_dx[0][i] * (du_dx*dphi_dx[0][j] + dv_dx*dphi_dx[1][j]) * weight_ * djac_;
//                 conveccao[2*i  ][2*j+1] += dphi_dx[1][i] * (du_dx*dphi_dx[0][j] + dv_dx*dphi_dx[1][j]) * weight_ * djac_;
//                 conveccao[2*i+1][2*j  ] += dphi_dx[0][i] * (du_dy*dphi_dx[0][j] + dv_dy*dphi_dx[1][j]) * weight_ * djac_;
//                 conveccao[2*i+1][2*j+1] += dphi_dx[1][i] * (du_dy*dphi_dx[0][j] + dv_dy*dphi_dx[1][j]) * weight_ * djac_;

//                 inercia[2*i  ][2*j  ] += dphi_dx[0][i] * dphi_dx[0][j] * weight_ * djac_;
//                 inercia[2*i  ][2*j+1] += dphi_dx[1][i] * dphi_dx[1][j] * weight_ * djac_;
//                 inercia[2*i+1][2*j  ] += dphi_dx[0][i] * dphi_dx[0][j] * weight_ * djac_;
//                 inercia[2*i+1][2*j+1] += dphi_dx[1][i] * dphi_dx[1][j] * weight_ * djac_;

//                 // visccc(2*i  ,2*j  ) += ((ddphi_dx(0,0)(i)+ddphi_dx(0,1)(i)) * (2. * ddphi_dx(0,0)(j) + ddphi_dx(1,1)(j))) * weight_ * djac_;
//                 // visccc(2*i  ,2*j+1) += ((ddphi_dx(0,0)(i)+ddphi_dx(0,1)(i)) * (ddphi_dx(0,1)(j))) * weight_ * djac_;
//                 // visccc(2*i+1,2*j  ) += ((ddphi_dx(1,1)(i)+ddphi_dx(0,1)(i)) * (ddphi_dx(0,1)(j))) * weight_ * djac_;
//                 // visccc(2*i+1,2*j+1) += ((ddphi_dx(1,1)(i)+ddphi_dx(0,1)(i)) * (2. * ddphi_dx(1,1)(j) + ddphi_dx(0,0)(j))) * weight_ * djac_;
//             }
//         }
//         index++;
//     }

//     double U_[18] = {};
//     double DU_[18] = {};

//     for (int i=0; i<6; i++){
//         U_[2*i  ] = (*nodes_)[connect_[i]] -> getVelocity(0);
//         U_[2*i+1] = (*nodes_)[connect_[i]] -> getVelocity(1);
//         DU_[2*i  ] = (*nodes_)[connect_[i]] -> getAcceleration(0);
//         DU_[2*i+1] = (*nodes_)[connect_[i]] -> getAcceleration(1);
//     };

//     double a1[18] = {};
//     double a2[18] = {};
//     double a3[18] = {};
//     double a4[18] = {};

//     for (int i = 0; i < 18; i++){
//         for (int j = 0; j < 18; j++){
//             a1[i] += lambda[i][j] * U_[j];
//             a2[i] += conveccao[i][j] * U_[j];
//             a3[i] += inercia[i][j] * DU_[j];
//             a4[i] += visccc[i][j] * U_[j];
//         }
//     }

//     le = 0.; //norm_2(prod(lambda,U_));
//     lc = 0.; //norm_2(prod(conveccao,U_));
//     li = 0.; //norm_2(prod(inercia,DU_));
//     lv = 0.; //norm_2(prod(visccc,U_));
    
//     for (int i = 0; i < 18; i++){
//         le += a1[i] * a1[i];
//         lc += a2[i] * a2[i];
//         li += a3[i] * a3[i];
//         lv += a4[i] * a4[i];
//     }


//     // //std::cout << "asd " << index_ << " " << le << " " << lc << " " << li << " " << tARLQ_ << std::endl;
//     //  // std::cout << "antes " << index_ << " " << tARLQ_ << std::endl;

//     double tp1 = std::fabs(le) / std::fabs(lc);
//     double tp2 = std::fabs(le) / std::fabs(li);
//     double tp3 = tp1 / visc_;
//     tp3 = le / (visc_*lv);

//     if (std::fabs(tp1) <= 1.e-10) tp1 = 1.e-10;
//     if (std::fabs(tp2) <= 1.e-10) tp2 = 1.e-10;
//     if (std::fabs(tp3) <= 1.e-10) tp3 = 1.e-10;

//     tARLQ_ = -1000.0 * k1 / std::sqrt(1. / (tp1*tp1) +
//                             1. / (tp2*tp2) +
//                             1. / (tp3*tp3));

//     // std::cout << "depois " << index_ << " " << le << " " << lv << " " << tARLQ_ << std::endl;

//     for (int i = DIM; i--; ) delete [] dphi_dx[i];
//     delete [] dphi_dx;
//     for (int i = DIM; i--; ) delete [] ainv_[i];
//     delete [] ainv_;

//     return;
// };

//------------------------------------------------------------------------------
//----------------------ELEMENT DIFFUSION/VISCOSITY MATRIX----------------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Element<DIM,DEG>::getElemMatrix(int &index, double **dphi_dx, double &tSUPG_, double &tPSPG_, double &tLSIC_, double &weight_, double &djac_, double **jacobianNRMatrix){

    double &dTime_ = parameters.getTimeStep();
    double &visc_ = parameters.getViscosity();
    double &dens_ = parameters.getDensity();
    double &alpha_f = parameters.getAlphaF();
    double &alpha_m = parameters.getAlphaM();
    double &gamma = parameters.getGamma();
    int &iTimeStep = parameters.getTimeInstant();

    //Velocity
    double u_[DIM], uPrev_[DIM], una_[DIM];
    interpolateVelocity(index, u_, uPrev_);
    for(int i=DIM; i--; ) una_[i] = alpha_f * u_[i] + (1. - alpha_f) * uPrev_[i];

    //Velocity Derivatives
    double du_dx[DIM][DIM], duprev_dx[DIM][DIM], duna_dx[DIM][DIM];
    interpolateVelDerivatives(dphi_dx, du_dx, duprev_dx);
    for (int i=DIM; i--; )
        for (int j = DIM; j--; )
            duna_dx[i][j] = alpha_f * du_dx[i][j] + (1. - alpha_f) * duprev_dx[i][j];

    //Mesh Velocity
    double umesh_[DIM], umeshPrev_[DIM], umeshna_[DIM];
    interpolateMeshVelocity(index, umesh_, umeshPrev_);
    for(int i=DIM; i--; ) umeshna_[i] = alpha_f * umesh_[i] + (1. - alpha_f) * umeshPrev_[i];

    double wna_ = alpha_f * intPointWeightFunction[index] + (1. - alpha_f) * intPointWeightFunctionPrev[index];

    // Trust me, it improves performance!
    double AGDT = alpha_f * gamma * dTime_; 
    double VAGDT = visc_ * AGDT; 
    double DAGDT = dens_ * AGDT; 
    double WJ = weight_ * djac_ * wna_;
    double DAM = dens_ * alpha_m;

    for (int i = nElNodes; i-- ; ){
        
        double shapeFi = DI -> phi_[i][index];
        double wSUPGi = 0.;
        for (int m=DIM; m--; ) wSUPGi += (una_[m] - umeshna_[m]) * dphi_dx[i][m];

        for (int j = nElNodes; j-- ; ){
            
            double shapeFj = DI -> phi_[j][index];
            double shapeFij = shapeFi * shapeFj;
            double wSUPGj = 0.;
            for (int m=DIM; m--; ) wSUPGj += (una_[m] - umeshna_[m]) * dphi_dx[j][m];
            
            //Mass matrix (use for both directions)
            double M = (shapeFij + wSUPGi * shapeFj * tSUPG_) * DAM;

            //Convection matrix
            double C = (wSUPGj * shapeFi + wSUPGi * wSUPGj * tSUPG_) * DAGDT;

            double aux1 = tSUPG_ * wSUPGi * shapeFj;
            double aux2 = tSUPG_ * shapeFj;

            for (int k = DIM; k--;  ){

                jacobianNRMatrix[DIM*i+k][DIM*j+k] += (M + C) * WJ;
                    
                double conv = 0.;
                for (int m = DIM; m--; ) conv += una_[m]*duna_dx[k][m];

                for (int l = DIM; l--; ){

                    //Diffusion matrix
                    double K = dphi_dx[i][l] * dphi_dx[j][k] * VAGDT;
                    if (k==l) for (int m = DIM; m--; ) K += dphi_dx[i][m] * dphi_dx[j][m] * VAGDT;

                    //Convection derivatives
                    double Cuu = (shapeFij * duna_dx[k][l] + 
                                  aux1 * duna_dx[k][l] +
                                  aux2 * conv * dphi_dx[i][l]) * DAGDT;

                    //LSIC
                    double KLS = dphi_dx[i][k] * dphi_dx[j][l] * tLSIC_ * DAGDT;

                    jacobianNRMatrix[DIM*i+k][DIM*j+l] += (K + KLS + Cuu) * WJ;
                }

                //SINAL DA PARCELA QUE MULTIPLICA O TSUPG ESTA COM SINAL TROCADO NA FORMULAÇAO DO TEZDUYAR
                //Gradient operator
                double Q_SUPG = - dphi_dx[i][k] * shapeFj + wSUPGi * dphi_dx[j][k] * tSUPG_;
                //Divergent operator
                double Q = dphi_dx[i][k] * shapeFj * AGDT;

                jacobianNRMatrix[DIM*i+k][DIM*(nElNodes)+j] += Q_SUPG * WJ;
                jacobianNRMatrix[DIM*(nElNodes)+j][DIM*i+k] += Q * WJ;


                //PSPG stabilization
                double H = dphi_dx[i][k] * shapeFj * tPSPG_ * alpha_m;
                double G = dphi_dx[i][k] * wSUPGj * tPSPG_ * AGDT;
                double Guu = 0.;
                // for (int m = DIM; m--; ) Guu += dphi_dx[m][i] * duna_dx[m][k] * shapeFj * tPSPG_ * AGDT;

                jacobianNRMatrix[DIM*(nElNodes)+j][DIM*i+k] += (H + G + Guu) * WJ;
            }

            double Q = 0.;
            for (int m = DIM; m--; ) Q += dphi_dx[i][m] * dphi_dx[j][m] * tPSPG_ / dens_;
            jacobianNRMatrix[DIM*(nElNodes)+j][DIM*(nElNodes)+i] += Q * WJ;
        };
    };

    return;
};

//------------------------------------------------------------------------------
//--------------------APPLY THE DIRICHLET BOUNDARY CONDITIONS-------------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Element<DIM,DEG>::setBoundaryConditions(double **jacobianNRMatrix, double *rhsVector){

    for (int i = nElNodes; i--; ){
        for (int k = DIM; k--; ){
            if (((*nodes_)[connect_[i]] -> getConstrains(k) == 1) ||
                ((*nodes_)[connect_[i]] -> getConstrains(k) == 3))  {
                for (int j = nLocDOF; j--; ){
                    jacobianNRMatrix[DIM*i+k][j] = 0.;
                    jacobianNRMatrix[j][DIM*i+k] = 0.;
                };
                jacobianNRMatrix[DIM*i+k][DIM*i+k] = 1.;
                rhsVector[DIM*i+k] = 0.;
            };
        }
    }

    return;
};

//------------------------------------------------------------------------------
//--------------------APPLY THE DIRICHLET BOUNDARY CONDITIONS-------------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Element<DIM,DEG>::setBoundaryConditionsLaplace(double **jacobianNRMatrix, double *rhsVector){

    for (int i = nElNodes; i--; ){
        for (int k = DIM; k--; ){
            if ((*nodes_)[connect_[i]] -> getConstrainsLaplace(k) == 1) {
                for (int j = nLocDOF; j--; ){
                    jacobianNRMatrix[DIM*i+k][j] = 0.;
                    jacobianNRMatrix[j][DIM*i+k] = 0.;
                };
                jacobianNRMatrix[DIM*i+k][DIM*i+k] = 1.;
                rhsVector[DIM*i+k] = 0.;
            };
        }
    }
    
    return;
};


//------------------------------------------------------------------------------
//---------------APPLY BOUNDARY CONDITIONS TO LAGRANGE MULTIPLIERS--------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Element<DIM,DEG>::setBoundaryConditionsLagrangeMultipliers(double** jacobianNRMatrix, double* rhsVector){

    for (int i = 0; i < nLocDOF; i++) rhsVector[i] = 0.;
    double U_[nLocDOF] = {};
    double &alpha_f = parameters.getAlphaF();

    for (int i = 0; i < nElNodes; i++)
        for (int k = 0; k < DIM; k++)
            U_[DIM*i+k] = alpha_f * (*nodes_)[connect_[i]] -> getVelocity(k) + (1. - alpha_f) * (*nodes_)[connect_[i]] -> getPreviousVelocity(k);
        
    for (int i = 0; i < nLocDOF; i++)
        for (int j = 0; j < nLocDOF; j++)
            rhsVector[i] -= jacobianNRMatrix[i][j] * U_[j];


    for (int i = nElNodes; i--; ){
        for (int k = DIM; k--; ){
            if ((*nodes_)[connect_[i]] -> getConstrains(k) == 1)  {
                for (int j = nLocDOF; j--; ){
                    jacobianNRMatrix[j][DIM*i+k] = 0.;
                };
                jacobianNRMatrix[DIM*i+k][DIM*i+k] = 1.;
                rhsVector[DIM*i+k] = 0.;
            };
        }
    }

    return;
};

//------------------------------------------------------------------------------
//-----------------------------RESIDUAL - RHS VECTOR----------------------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Element<DIM,DEG>::getResidualVector(int &index, double **dphi_dx, double &tSUPG_, double &tPSPG_, double &tLSIC_, double &weight_, double &djac_, double *rhsVector){

    double &dTime_ = parameters.getTimeStep();
    double &visc_ = parameters.getViscosity();
    double &dens_ = parameters.getDensity();
    double &alpha_f = parameters.getAlphaF();
    double &alpha_m = parameters.getAlphaM();
    double &gamma = parameters.getGamma();
    double* fieldForce = parameters.getFieldForce();

    //Velocity
    double u_[DIM], uPrev_[DIM], una_[DIM];
    interpolateVelocity(index, u_, uPrev_);
    for(int i=DIM; i--; ) una_[i] = alpha_f * u_[i] + (1. - alpha_f) * uPrev_[i];

    //Velocity Derivatives
    double du_dx[DIM][DIM], duprev_dx[DIM][DIM], duna_dx[DIM][DIM];
    interpolateVelDerivatives(dphi_dx, du_dx, duprev_dx);
    for (int i=DIM; i--; )
        for (int j=DIM; j--; )
            duna_dx[i][j] = alpha_f * du_dx[i][j] + (1. - alpha_f) * duprev_dx[i][j];

    //Acceleration
    double a_[DIM], aPrev_[DIM], am_[DIM];
    interpolateAcceleration(index, a_, aPrev_);
    for(int i=DIM; i--; ) am_[i] = alpha_m * a_[i] + (1. - alpha_m) * aPrev_[i];

    //Mesh Velocity
    double umesh_[DIM], umeshPrev_[DIM], umeshna_[DIM];
    interpolateMeshVelocity(index, umesh_, umeshPrev_);
    for(int i=DIM; i--; ) umeshna_[i] = alpha_f * umesh_[i] + (1. - alpha_f) * umeshPrev_[i];

    //Pressure
    double p_;
    double dp_dx[DIM] = {};
    interpolatePressure(index, dphi_dx, p_, dp_dx);

    double wna_ = alpha_f * intPointWeightFunction[index] + (1. - alpha_f) * intPointWeightFunctionPrev[index];
    
    double WJ = weight_ * djac_ * wna_;

    double divrU = 0.;
    for (int l=DIM; l--; ) divrU += duna_dx[l][l];

    for (int i = nElNodes; i--; ){
        double shapeFi = DI -> phi_[i][index];

        for (int k = DIM; k--; ){

            //Mass + SUPG mass
            double m = shapeFi * am_[k] * dens_;
            for (int l=DIM; l--; ) m += (una_[l] - umeshna_[l]) * dphi_dx[i][l] * am_[k] * tSUPG_ * dens_;
                    
            //Viscosity
            double K = 0.;
            for (int l=DIM; l--; ) K += dphi_dx[i][l] * duna_dx[k][l] * visc_;
            for (int l=DIM; l--; ) K += dphi_dx[i][l] * duna_dx[l][k] * visc_;

            //LSIC
            double KLS = dphi_dx[i][k] * divrU * tLSIC_ * dens_;

            //Convection + SUPG
            double C = 0.;
            for (int l=DIM; l--; ) C += duna_dx[k][l] * (una_[l] - umeshna_[l]) * shapeFi * dens_;
            double conv = 0.;
            for (int l=DIM; l--; ) conv += (una_[l] - umeshna_[l]) * dphi_dx[i][l];
            for (int l=DIM; l--; ) C += conv * (una_[l] - umeshna_[l]) * duna_dx[k][l] * tSUPG_ * dens_;

            //Pressure + SUPG
            double P = - (dphi_dx[i][k] * p_);
            for (int l=DIM; l--; ) P += dphi_dx[i][l] * (una_[l] - umeshna_[l]) * dp_dx[k] * tSUPG_;

            //External force
            double F = fieldForce[k] * shapeFi * dens_;
            
            rhsVector[DIM*i+k] += (-m -K - P - C - KLS + F) * WJ;
        }

        double Q = divrU * shapeFi;
        for (int l=DIM; l--; ) Q += dphi_dx[i][l] * dp_dx[l] * tPSPG_ / dens_
                                  + dphi_dx[i][l] * am_[l] * tPSPG_ 
                                  + dphi_dx[i][l] * fieldForce[l] * tPSPG_;
        for (int k=DIM; k--; )
            for (int l=DIM; l--; )
                Q += dphi_dx[i][k] * (una_[l] - umeshna_[l]) * duna_dx[k][l] * tPSPG_;

        rhsVector[DIM*(nElNodes)+i] += -Q * WJ;
                            
    };

    return;
};

//------------------------------------------------------------------------------
//-----------------------------RESIDUAL - RHS VECTOR----------------------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Element<DIM,DEG>::getResidualVectorLaplace(double *rhsVector){
    
    double U_[nLocDOF] = {};

    for (int i = 0; i < nElNodes; i++){
        double *x_up = (*nodes_)[connect_[i]] -> getUpdatedCoordinates();        
        for (int k = 0; k < DIM; k++)
            if ((*nodes_)[connect_[i]] -> getConstrainsLaplace(k) == 1)
                U_[DIM*i+k] = x_up[k] - (*nodes_)[connect_[i]] -> getCoordinateValue(k);
            
    };

    return;
};

//------------------------------------------------------------------------------
//---------------------------ELEMENT LAPLACIAN MATRIX---------------------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Element<DIM,DEG>::getElemLaplMatrix(double &weight_, double &djac_, double **dphi_dx, double **jacobianNRMatrix){

    double WJM = weight_ * djac_ * meshMovingParameter;
    for (int i = 0; i < nElNodes; i++){
        for (int j = 0; j < nElNodes; j++){  
            for (int k = DIM; k--; ){
                for (int l = DIM; l--; ){
                    //Stiffness matrix
                    // double K = dphi_dx[l][i] * dphi_dx[k][j];
                    double K = 0.;
                    if (k==l) for (int m = DIM; m--; ) K += dphi_dx[i][m] * dphi_dx[j][m];
                    
                    jacobianNRMatrix[DIM*i+k][DIM*j+l] += K * WJM;
                }
            }
        }
    }

    return;
};

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//----------------MOUNTS EACH TYPE OF INCOMPRESSIBLE FLOW PROBEM----------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//-----------------------TRANSIENT NAVIER-STOKES PROBEM-------------------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Element<DIM,DEG>::getTransientNavierStokes(double **jacobianNRMatrix, double *rhsVector){

    double xsi[DIM] = {};
    
    double **dphi_dx;
    dphi_dx = new double*[nElNodes];
    for (int i = nElNodes; i--; ) dphi_dx[i] = new double[DIM];

    double **ainv_;
    ainv_ = new double*[DIM];
    for (int i = DIM; i--; ) ainv_[i] = new double[DIM];

    ShapeFunction           shapeQuad;
    int index = 0;
    NormalQuad nQuad = NormalQuad();

    double tSUPG_;
    double tPSPG_;
    double tLSIC_;  

    for(double* it = nQuad.begin(); it != nQuad.end(); it++){

        //Defines the integration points adimentional coordinates
        for (int k = 0; k < DIM; k++) xsi[k] = nQuad.PointList(index,k);

        //Returns the quadrature integration weight
        double weight_ = nQuad.WeightList(index);

        double djac_ = 0.;
        //Computes the jacobian matrix
        getJacobianMatrix(xsi, ainv_, djac_);

        //Computes spatial derivatives
        getSpatialDerivatives(xsi, ainv_, dphi_dx);

        //Compute Stabilization Parameters
        getParameterSUPG(index, tSUPG_, tPSPG_, tLSIC_, dphi_dx);

        //Computes the element diffusion/viscosity matrix
        getElemMatrix(index, dphi_dx, tSUPG_, tPSPG_, tLSIC_, weight_, djac_, jacobianNRMatrix);

        //Computes the RHS vector
        getResidualVector(index, dphi_dx, tSUPG_, tPSPG_, tLSIC_, weight_, djac_, rhsVector); 

        index++;        
    };  
    
    //Apply boundary conditions
    setBoundaryConditions(jacobianNRMatrix, rhsVector);

    for (int i = nElNodes; i--; ) delete [] dphi_dx[i];
    delete [] dphi_dx;
    for (int i = DIM; i--; ) delete [] ainv_[i];
    delete [] ainv_;

    return;
};

//------------------------------------------------------------------------------
//----------------------------STEADY LAPLACE PROBEM-----------------------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Element<DIM,DEG>::getSteadyLaplace(double** jacobianNRMatrix, double* rhsVector){

    double xsi[DIM] = {};
    ShapeFunction           shapeQuad;
    
    double **dphi_dx;
    dphi_dx = new double*[nElNodes];
    for (int i = nElNodes; i--; ) dphi_dx[i] = new double[DIM];

    int index = 0;
    NormalQuad nQuad = NormalQuad();

    double **ainv_;
    ainv_ = new double*[DIM];
    for (int i = DIM; i--; ) ainv_[i] = new double[DIM];

    // ublas::bounded_matrix<double, 3, 3 > hooke;

    // hooke.clear();
    
    // // For EPT
    // double elastic_ = 10000.;
    // double poisson_ = 0.45;
    // double k = elastic_ / (1. - poisson_ * poisson_);
    // hooke(0,0) = k;
    // hooke(0,1) = k * poisson_;
    // hooke(1,0) = k * poisson_;
    // hooke(1,1) = k;
    // hooke(2,2) = k * (1. - poisson_) * 0.5;

    
    for(double* it = nQuad.begin(); it != nQuad.end(); it++){
        
        //Defines the integration points adimentional coordinates
        for (int k = 0; k < DIM; k++) xsi[k] = nQuad.PointList(index,k);

        //Returns the quadrature integration weight
        double weight_ = nQuad.WeightList(index);

        double djac_ = 0.;
        //Computes the jacobian matrix
        getJacobianMatrix(xsi, ainv_, djac_);

        //Computes spatial derivatives
        getSpatialDerivatives(xsi, ainv_, dphi_dx);

        getElemLaplMatrix(weight_, djac_, dphi_dx, jacobianNRMatrix);

        index++;        
    };  
    
    //Computes the RHS vector
    getResidualVectorLaplace(rhsVector);

    //Apply boundary conditions
    setBoundaryConditionsLaplace(jacobianNRMatrix, rhsVector);

    for (int i = nElNodes; i--; ) delete [] dphi_dx[i];
    delete [] dphi_dx;
    for (int i = DIM; i--; ) delete [] ainv_[i];
    delete [] ainv_;

    return;
};

//------------------------------------------------------------------------------
//----------------------------STEADY LAPLACE PROBEM-----------------------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Element<DIM,DEG>::getSteadyLaplace2(double** jacobianNRMatrix, double* rhsVector){

    double xsi[DIM] = {};
    double **dphi;
    dphi = new double*[nElNodes];
    for (int i = nElNodes; i--; ) dphi[i] = new double[DIM];
    
    ShapeFunction           shapeQuad;
    // ShapeFunctionDerivative dphi_dx;
    int index = 0;
    NormalQuad nQuad = NormalQuad();

    double &dTime_ = parameters.getTimeStep();
        
    for(double* it = nQuad.begin(); it != nQuad.end(); it++){
        
        //Defines the integration points adimentional coordinates
        for (int k = 0; k < DIM; k++) xsi[k] = nQuad.PointList(index,k);

        //Computes the velocity shape functions
        shapeQuad.evaluateGradient(xsi,dphi);

        //Returns the quadrature integration weight
        double weight_ = nQuad.WeightList(index);

        //COMPUTE A0
        double dx_dxsi[DIM][DIM];
        for (int i = 0; i < DIM; i++)
            for (int j = 0; j < DIM; j++)
                dx_dxsi[i][j] = 0.0;

        for (int i = 0; i < nElNodes; ++i){
            double* initialCoord = (*nodes_)[connect_[i]] -> getInitialCoordinates();
            for (int k = 0; k < DIM; k++)
                for (int l = 0; l < DIM; l++)
                    dx_dxsi[k][l] += initialCoord[k] * dphi[i][l];
        }

        double j0 = dx_dxsi[0][0] * dx_dxsi[1][1] - dx_dxsi[0][1] * dx_dxsi[1][0];
       
        //dxsi_dx
        double dxsi_dx[DIM][DIM] = {};

        dxsi_dx[0][0] = dx_dxsi[1][1] / j0;
        dxsi_dx[0][1] = -dx_dxsi[0][1] / j0;
        dxsi_dx[1][0] = -dx_dxsi[1][0] / j0;
        dxsi_dx[1][1] = dx_dxsi[0][0] / j0;

        //dphi_dx
        double dphi_dx[nElNodes][DIM]; 
        for (int i = 0; i < nElNodes; i++)
        {
            dphi_dx[i][0] = dxsi_dx[0][0] * dphi[i][0] + dxsi_dx[1][0] * dphi[i][1];
            dphi_dx[i][1] = dxsi_dx[0][1] * dphi[i][0] + dxsi_dx[1][1] * dphi[i][1];
        }


        //COMPUTE A1
        double dy_dxsi[DIM][DIM]; //row = cartesian, column = parametric
        for (int i = 0; i < DIM; i++)
            for (int j = 0; j < DIM; j++)
                dy_dxsi[i][j] = 0.0;
        
        for (int i = 0; i < nElNodes; i++){
            double* currentCoord = (*nodes_)[connect_[i]] -> getCoordinates();

            for (int k = 0; k < DIM; k++)
                for (int l = 0; l < DIM; l++)
                    dy_dxsi[k][l] += currentCoord[k] * dphi[i][l];
        }

        //dy_dx
        double dy_dx[DIM][DIM] = {};
        for (int i = 0; i < DIM; i++)
            for (int j = 0; j < DIM; j++)
                dy_dx[i][j] = dy_dxsi[i][0] * dxsi_dx[0][j] + dy_dxsi[i][1] * dxsi_dx[1][j];

        //jacobian
        double jac = dy_dx[0][0] * dy_dx[1][1] - dy_dx[0][1] * dy_dx[1][0];

        //Green-Lagrange strain tensor
        double E[DIM][DIM] = {};
        for (int i = 0; i < DIM; i++)
            for (int j = 0; j < DIM; j++)
                E[i][j] = 0.5 * (dy_dx[0][i] * dy_dx[0][j] + dy_dx[1][i] * dy_dx[1][j]);
        E[0][0] -= 0.5; E[1][1] -= 0.5;

        //Second Piola-Kirchhoff stress tensor
        double S[DIM][DIM] = {};
        double young = 1. / j0;
        double poisson = 0.3;
        S[0][0] = young / (1.0-(poisson*poisson)) * (E[0][0] + poisson * E[1][1]);
        S[0][1] = 2.0 * (young / (2.0 * (1.0+poisson))) * E[0][1];
        S[1][0] = S[0][1];
        S[1][1] = young / (1.0-(poisson*poisson)) * (E[1][1] + poisson * E[0][0]);
        //     S[0][0] = (young / ( (1.0+poisson) * (1.0-2.0*poisson))) * ((1.0-poisson) * E[0][0] + poisson * E[1][1]);
        //     S[1][1] = (young / ( (1.0+poisson) * (1.0-2.0*poisson))) * ((1.0-poisson) * E[1][1] + poisson * E[0][0]);


        //element rhs vector
        for (int a = 0; a < nElNodes; a++){
            for (int k = 0; k < DIM; k++){
                double dE_dyak[DIM][DIM];
                for (int i = 0; i < DIM; i++)
                    for (int j = 0; j < DIM; j++)
                        dE_dyak[i][j] = 0.5 * (dphi_dx[a][i] * dy_dx[k][j] + dy_dx[k][i] * dphi_dx[a][j]);
                
                //internal force
                double f = 0.0;
                for (int i = 0; i < DIM; i++)
                    for (int j = 0; j < DIM; j++)
                        f += S[i][j] * dE_dyak[i][j];

                double vel = 0.0;
                for (int i = 0; i < nElNodes; i++)
                    vel += DI -> phi_[i][index] * (*nodes_)[connect_[i]]->getMeshVelocity(k);

                double c =  DI -> phi_[a][index] * vel*0;

                //domain force
                // double b;
                // (k==1) ? b = phi(a) * density * gravity : 0.0;

                rhsVector[2 * a + k] -= (f+c) * weight_ * j0;

                //element tangent matrix
                for (int b = 0; b < nElNodes; b++){
                    for (int l = 0; l < DIM; l++){
                        double dE_dybl[DIM][DIM];
                        for (int i = 0; i < DIM; i++)
                            for (int j = 0; j < DIM; j++)
                                dE_dybl[i][j] = 0.5 * (dphi_dx[b][i] * dy_dx[l][j] + dy_dx[l][i] * dphi_dx[b][j]);

                        double d2E_dyakbl[DIM][DIM];
                        if (k==l)
                            for (int i = 0; i < DIM; i++)
                                for (int j = 0; j < DIM; j++)
                                    d2E_dyakbl[i][j] = 0.5 * (dphi_dx[a][i] * dphi_dx[b][j] + dphi_dx[b][i] * dphi_dx[a][j]);
                        else
                            for (int i = 0; i < DIM; i++)
                                for (int j = 0; j < DIM; j++)
                                    d2E_dyakbl[i][j] = 0.0; 

                        double dS_dybl[DIM][DIM];
                        dS_dybl[0][0] = young / (1.0-(poisson*poisson)) * (dE_dybl[0][0] + poisson * dE_dybl[1][1]);
                        dS_dybl[0][1] = 2.0 * (young / (2.0 * (1.0+poisson))) * dE_dybl[0][1];
                        dS_dybl[1][0] = 2.0 * (young / (2.0 * (1.0+poisson))) * dE_dybl[1][0];
                        dS_dybl[1][1] = young / (1.0-(poisson*poisson)) * (dE_dybl[1][1] + poisson * dE_dybl[0][0]);

                        //elastic and geometric componentes of tangent matrix
                        double e = 0.0;
                        for (int i = 0; i < DIM; i++)
                            for (int j = 0; j < DIM; j++)
                                e += dS_dybl[i][j] * dE_dyak[i][j] + S[i][j] * d2E_dyakbl[i][j];

                        //mass matrix
                        double m;
                        (k==l)? m = (1.0 / (0.25 * dTime_)) *0* DI -> phi_[a][index] * DI -> phi_[b][index] : m = 0.0;

                        jacobianNRMatrix[2 * a + k][2 * b + l] += (e+m) * j0 * weight_;
                    }
                }
            }
        }

        // for (int i = 0; i < 6; i++){
        //     for (int j = 0; j < 6; j++){        
        //         jacobianNRMatrix(2*i  ,2*j  ) += (dphi_dx[0][i] * dphi_dx[0][j] +
        //                                           dphi_dx[1][i] * dphi_dx[1][j]) 
        //                                     * weight_ * djac_ * meshMovingParameter;
        //         jacobianNRMatrix(2*i+1,2*j+1) += (dphi_dx[0][i] * dphi_dx[0][j] +
        //                                           dphi_dx[1][i] * dphi_dx[1][j]) 
        //                                     * weight_ * djac_ * meshMovingParameter;
        //     };
        // };

        index++;        
    };  

    //Apply boundary conditions
    setBoundaryConditionsLaplace(jacobianNRMatrix, rhsVector);

    for (int i = nElNodes; i--; ) delete [] dphi[i];
    delete [] dphi;

    return;
};

//------------------------------------------------------------------------------
//----------------------------STEADY LAPLACE PROBEM-----------------------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Element<DIM,DEG>::getLagrangeMultipliersSameMesh(double **lagrMultMatrix, double *lagrMultVector, double *rhsVector){

    double xsi[DIM] = {};
    
    double **dphi_dx;
    dphi_dx = new double*[nElNodes];
    for (int i = nElNodes; i--; ) dphi_dx[i] = new double[DIM];

    double **ainv_;
    ainv_ = new double*[DIM];
    for (int i = DIM; i--; ) ainv_[i] = new double[DIM];

    QuadShapeFunction<DIM,DEG> shapeQuad;
    int index = 0;
    NormalQuad nQuad = NormalQuad();

    double tSUPG_; double tPSPG_; double tLSIC_;

    double &k1 = parameters.getArlequinK1();
    double &k2 = parameters.getArlequinK2();
    double &alpha_f = parameters.getAlphaF();
    double &gamma = parameters.getGamma();
    double &dTime_ = parameters.getTimeStep();
    
    for(double* it = nQuad.begin(); it != nQuad.end(); it++){
        
        //Defines the integration points adimentional coordinates
        for (int k = 0; k < DIM; k++) xsi[k] = nQuad.PointList(index,k);

        //Returns the quadrature integration weight
        double weight_ = nQuad.WeightList(index);
        double djac_ = 0.;
        //Computes the jacobian matrix
        getJacobianMatrix(xsi, ainv_, djac_);

        getSpatialDerivatives(xsi, ainv_, dphi_dx);

        //Velocity
        double u_[DIM], uPrev_[DIM], una_[DIM];
        interpolateVelocity(index, u_, uPrev_);
        for(int i=DIM; i--; ) una_[i] = alpha_f * u_[i] + (1. - alpha_f) * uPrev_[i];

        //Velocity Derivatives
        double du_dx[DIM][DIM], duprev_dx[DIM][DIM], duna_dx[DIM][DIM];
        interpolateVelDerivatives(dphi_dx, du_dx, duprev_dx);
        for (int i=DIM; i--; )
            for (int j = DIM; j--; )
                duna_dx[i][j] = alpha_f * du_dx[i][j] + (1. - alpha_f) * duprev_dx[i][j];

        //Lagrange Multiplier
        double lagM_[DIM];
        interpolateLagMultiplier(index, lagM_);

        //Lagrange Multiplier Derivatives
        double dL_dx[DIM][DIM];
        interpolateLagMultiplierDerivatives(dphi_dx, dL_dx);
        
        double WJ = weight_ * djac_;
        for (int i = 0; i < nElNodes; i++){
            for (int j = 0; j < nElNodes; j++){
                double l2 = DI -> phi_[i][index] * DI -> phi_[j][index] * WJ * k1;
                for (int k = 0; k < DIM; k++){
                    // L2 COUPLING OPERATOR
                    lagrMultMatrix[DIM*i+k][DIM*j+k] -= l2;
                    for (int l = 0; l < DIM; l++){
                        //H1 COUPLING OPERATOR
                        double K = dphi_dx[i][l] * dphi_dx[j][k];
                        if (k==l) for (int m = DIM; m--; ) K += dphi_dx[i][m] * dphi_dx[j][m];

                        lagrMultMatrix[DIM*i+k][DIM*j+l] -= K * WJ * k2;
                    };
                };
            };
            // Lagrange multipliers residual
            for (int k = 0; k < DIM; k++){
                double L2 = lagM_[k] * DI -> phi_[i][index] * k1;

                double H1 = 0.;
                for (int l=DIM; l--; ) H1 += dphi_dx[i][l] * dL_dx[k][l] * k2;
                for (int l=DIM; l--; ) H1 += dphi_dx[i][l] * dL_dx[l][k] * k2;

                lagrMultVector[DIM*i+k] += (L2 + H1) * WJ;

                double L2u = una_[k] * DI -> phi_[i][index] * k1;

                double H1u = 0.;
                for (int l=DIM; l--; ) H1u += dphi_dx[i][l] * duna_dx[k][l] * k2;
                for (int l=DIM; l--; ) H1u += dphi_dx[i][l] * duna_dx[l][k] * k2;

                rhsVector[DIM*i+k] += (L2u + H1u) * WJ;
            };
        };
        index++; 
    }; 

    for (int i = nElNodes; i--; ) delete [] dphi_dx[i];
    delete [] dphi_dx;
    for (int i = DIM; i--; ) delete [] ainv_[i];
    delete [] ainv_;

    return;
};


//------------------------------------------------------------------------------
//----------------------------STEADY LAPLACE PROBEM-----------------------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Element<DIM,DEG>::getLagrangeMultipliersSUPG_PSPG_SameMesh(double **jacobianNRMatrix, double *rhsVector){

    // double xsi[DIM] = {};
    // double phi_[nElNodes] = {};
    
    // double **dphi_dx;
    // dphi_dx = new double*[DIM];
    // for (int i = DIM; i--; ) dphi_dx[i] = new double[nElNodes];

    // double **ainv_;
    // ainv_ = new double*[DIM];
    // for (int i = DIM; i--; ) ainv_[i] = new double[DIM];  

    // QuadShapeFunction<DIM,DEG> shapeQuad;
    // int index = 0;
    // NormalQuad nQuad = NormalQuad();

    // double tSUPG_; double tPSPG_; double tLSIC_;

    // double &dens_ = parameters.getDensity();
    // double &alpha_f = parameters.getAlphaF();
    // double &k1 = parameters.getArlequinK1();
    
    // for(double* it = nQuad.begin(); it != nQuad.end(); it++){
        
    //     //Defines the integration points adimentional coordinates
    //     for (int k = 0; k < DIM; k++) xsi[k] = nQuad.PointList(index,k);

    //     //Computes the velocity shape functions
    //     shapeQuad.evaluate(xsi,phi_);

    //     //Returns the quadrature integration weight
    //     double weight_ = nQuad.WeightList(index);

    //     //Computes the jacobian matrix
    //     getJacobianMatrix(xsi, ainv_);

    //     getSpatialDerivatives(xsi, ainv_, dphi_dx);
        
    //     getVelAndDerivatives(phi_, dphi_dx);

    //     getParameterSUPG(tSUPG_, tPSPG_, tLSIC_, phi_, dphi_dx);

    //     double una_ = alpha_f * u_ + (1. - alpha_f) * uPrev_;
    //     double vna_ = alpha_f * v_ + (1. - alpha_f) * vPrev_;

    //     double umeshna_ = alpha_f * umesh_ + (1. - alpha_f) * umeshPrev_;
    //     double vmeshna_ = alpha_f * vmesh_ + (1. - alpha_f) * vmeshPrev_;

    //     double wna_ = alpha_f * intPointWeightFunction[index] + (1. - alpha_f) * intPointWeightFunctionPrev[index];
        
    //     for (int i = 0; i < 6; i++){
    //         for (int j = 0; j < 6; j++){        

    //             //SUPG STABILIZATION TERM
    //             double LM = 0.;
    //             // LM = - ((una_ - umeshna_) * dphi_dx[0][i] + (vna_ - vmeshna_) * dphi_dx[1][i])
    //             //     * phi_[j] * tSUPG_;
            
    //             jacobianNRMatrix[2*i  ][2*j  ] += LM * weight_ * djac_ * k1;
    //             jacobianNRMatrix[2*i+1][2*j+1] += LM * weight_ * djac_ * k1;

    //             //PSPG STABILIZATION TERM
    //             double Lx = 0.;
    //             double Ly = 0.;
    //             double Tx = 0.;
    //             double Ty = 0.;

    //             // Lx = dphi_dx[0][i] * phi_[j] * tPSPG_ / dens_ * k1;
    //             // Ly = dphi_dx[1][i] * phi_[j] * tPSPG_ / dens_ * k1;

    //             jacobianNRMatrix[2*i  ][12+j] += Lx * weight_ * djac_;
    //             jacobianNRMatrix[2*i+1][12+j] += Ly * weight_ * djac_;
    //         };

    //         //SUPG AND PSPG STABILIZATION TERMS
    //         double LMx = 0.;
    //         double LMy = 0.;
    //         double LMp = 0.;
            
    //         // LMx = (((una_ - umeshna_) * dphi_dx[0][i] + (vna_ - vmeshna_) * dphi_dx[1][i])
    //         //        * lagMx_) * tSUPG_;
    //         // LMy = (((una_ - umeshna_) * dphi_dx[0][i] + (vna_ - vmeshna_) * dphi_dx[1][i])
    //         //        * lagMy_) * tSUPG_;
    //         // LMp = (dphi_dx[0][i] * lagMx_ + dphi_dx[1][i] * lagMy_) * tPSPG_/ dens_;

    //         rhsVector[2*i  ] += (-LMx) * weight_ * djac_ * k1;
    //         rhsVector[2*i+1] += (-LMy) * weight_ * djac_ * k1;
    //         rhsVector[12+i] += (-LMp) * weight_ * djac_ * k1;

    //     };         
        
    //     index++;        
    // }; 

    // for (int i = DIM; i--; ) delete [] dphi_dx[i];
    // delete [] dphi_dx;
    // for (int i = DIM; i--; ) delete [] ainv_[i];
    // delete [] ainv_;

    return;
};

//------------------------------------------------------------------------------
//----------------------------STEADY LAPLACE PROBEM-----------------------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Element<DIM,DEG>::getLagrangeMultipliersArlequinSameMesh(double **arlequinStab, double **laplMatrix, double *arlequinStabVector){

    double xsi[DIM] = {};
    
    double **dphi_dx;
    dphi_dx = new double*[nElNodes];
    for (int i = nElNodes; i--; ) dphi_dx[i] = new double[DIM];

    double **ainv_;
    ainv_ = new double*[DIM];
    for (int i = DIM; i--; ) ainv_[i] = new double[DIM];  

    QuadShapeFunction<DIM,DEG> shapeQuad;
    int index = 0;
    NormalQuad nQuad = NormalQuad();

    double tSUPG_, tPSPG_, tLSIC_, tARLQ_;

    double &dens_ = parameters.getDensity();
    double &alpha_f = parameters.getAlphaF();
    double &k1 = parameters.getArlequinK1();

    for(double* it = nQuad.begin(); it != nQuad.end(); it++){
        //Defines the integration points adimentional coordinates
        for (int k = 0; k < DIM; k++) xsi[k] = nQuad.PointList(index,k);

        //Returns the quadrature integration weight
        double weight_ = nQuad.WeightList(index);

        double djac_ = 0.; 
        //Computes the jacobian matrix
        getJacobianMatrix(xsi, ainv_, djac_);

        getSpatialDerivatives(xsi, ainv_, dphi_dx);
        
        getParameterArlequin(index, tARLQ_, tSUPG_, tPSPG_, tLSIC_, dphi_dx);

        double wna_ = alpha_f * intPointWeightFunction[index] + (1. - alpha_f) * intPointWeightFunctionPrev[index];
        
        //Lagrange Multiplier Derivatives
        double dL_dx[DIM][DIM];
        interpolateLagMultiplierDerivatives(dphi_dx, dL_dx);


        for (int i = 0; i < nElNodes; i++){
            for (int j = 0; j < nElNodes; j++){        

                //ARLEQUIN STABILIZATION TERMS
                double AM = 0.;
                double Lpx = 0.; double Lpy = 0.;
                double LC = 0.; double LL = 0.;

                // AM = phi_[i] * phi_[j] * tARLQ_ * intPointWeightFunction(index) * alpha_m;

                // LL = -2 * phi_[i] * phi_[j] * tARLQ_ / dens_;
                for (int m = DIM; m--; ) LL += dphi_dx[i][m] * dphi_dx[j][m] * tARLQ_ / dens_;

                // Lpx = -(dphi_dx[0][i] * dp_dxx + dphi_dx[1][i] * dp_dxy) * intPointWeightFunction(index)
                //      * tARLQ_ / dens_;
                // Lpy = -(dphi_dx[0][i] * dp_dxy + dphi_dx[1][i] * dp_dyy) * intPointWeightFunction(index)
                //      * tARLQ_ / dens_;

                // LC = phi_[i] * ((una_ - umesh_) * dphi_dx[0][i] + (vna_ - vmesh_) * dphi_dx[1][i]) * phi_[j] * tARLQ_ * intPointWeightFunction(index);

                for (int k = DIM; k--; )
                    arlequinStab[DIM*i+k][DIM*j+k] += LL * weight_ * djac_;
                

                // LC = -(dphi_dx[0][i]*(du_dx*dphi_dx[0][j] + dv_dx*dphi_dx[1][j]) +
                //        dphi_dx[1][i]*(du_dy*dphi_dx[0][j] + dv_dy*dphi_dx[1][j]) + 
                //        dphi_dx[0][i]*(u_*ddphi_dx(0,0)(j) + v_*ddphi_dx(0,1)(j)) + 
                //        dphi_dx[1][i]*(u_*ddphi_dx(1,0)(j) + v_*ddphi_dx(1,1)(j))) * tARLQ_ * intPointWeightFunction(index);

                // Lpx = 0.; Lpy = 0.;

                // Lpx = (dphi_dx[0][i] * ddphi_dx(0,0)(j) + 
                //        dphi_dx[1][i] * ddphi_dx(0,1)(j)) * tARLQ_ * intPointWeightFunction(index);
                // Lpy = (dphi_dx[0][i] * ddphi_dx(1,0)(j) + 
                //        dphi_dx[1][i] * ddphi_dx(1,1)(j)) * tARLQ_ * intPointWeightFunction(index);




                // laplMatrix[2*i  ][2*j  ] += (LC + AM) * weight_ * djac_;
                // laplMatrix[2*i+1][2*j+1] += (LC + AM) * weight_ * djac_;

                // laplMatrix[2*i  ][12+j] += Lpx * weight_ * djac_;
                // laplMatrix[2*i+1][12+j] += Lpy * weight_ * djac_;
                // laplMatrix[12+j][2*i  ] += Lpx * weight_ * djac_;
                // laplMatrix[12+j][2*i+1] += Lpy * weight_ * djac_;

            };

            //ARLEQUIN STABILIZATION TERMS
            double Amx = 0.; double Amy = 0.;
            double LCx = 0.; double LCy = 0.;
            double LPx = 0.; double LPy = 0.;
             double LLy = 0.;

            // if (iTimeStep > 5){
                // Amx = -phi_[i] * axm_ * tARLQ_;
                // Amy = -phi_[i] * aym_ * tARLQ_;

            // }

            //if (iTimeStep < 10){
            for (int k = DIM; k--; ){
                double LLx = 0.;
                for (int m = DIM; m--; ) LLx -= dphi_dx[i][m] * dL_dx[k][m]/wna_ * tARLQ_ / dens_;

                arlequinStabVector[DIM*i+k] += (Amx + LLx) * weight_ * djac_ * wna_;
            }

                // LLx = -(dphi_dx[0][i] * (dLx_dx/wna_) + dphi_dx[1][i] * (dLx_dy/wna_)) * tARLQ_ / dens_;
                // LLy = -(dphi_dx[0][i] * (dLy_dx/wna_) + dphi_dx[1][i] * (dLy_dy/wna_)) * tARLQ_ / dens_;
           //}else{
                // LLx = -(dphi_dx[0][i] * (dLx_dx/intPointWeightFunction(index) - duna_dx*duna_dx - dvna_dx*duna_dy) + 
                //         dphi_dx[1][i] * (dLx_dy/intPointWeightFunction(index) - duna_dy*duna_dx - dvna_dy*duna_dy)) * tARLQ_ / dens_;
                // LLy = -(dphi_dx[0][i] * (dLy_dx/intPointWeightFunction(index) - duna_dx*dvna_dx - dvna_dx*dvna_dy) + 
                //         dphi_dx[1][i] * (dLy_dy/intPointWeightFunction(index) - duna_dy*dvna_dx - dvna_dy*dvna_dy)) * tARLQ_ / dens_;  
                // LLx = -(dphi_dx[0][i] * (dLx_dx/intPointWeightFunction(index) - dp_dxx) + 
                //         dphi_dx[1][i] * (dLx_dy/intPointWeightFunction(index) - dp_dxy)) * tARLQ_ / dens_;
                // LLy = -(dphi_dx[0][i] * (dLy_dx/intPointWeightFunction(index) - dp_dxy) + 
                //         dphi_dx[1][i] * (dLy_dy/intPointWeightFunction(index) - dp_dyy)) * tARLQ_ / dens_;  

                // if (iTimeStep > 10){
                //     LLx +=  -(dphi_dx[0][i]*dax_dx + dphi_dx[1][i]*dax_dy) * tARLQ_;
                //     LLx +=  -(dphi_dx[0][i]*day_dx + dphi_dx[1][i]*day_dy) * tARLQ_ ;
                // };

            //}
            
            // arlequinStabVector[2*i  ] += (Amx + LLx) * weight_ * djac_ * wna_;
            // arlequinStabVector[2*i+1] += (Amy + LLy) * weight_ * djac_ * wna_;
        };         
        
        index++;        
    };  

    for (int i = 0; i < nElNodes; i++){
        for (int k = DIM; k--; ){
            if ((*nodes_)[connect_[i]] -> getConstrains(k) == 1) {
                for (int j = 0; j < nLocDOF; j++){
                    arlequinStab[DIM*i+k][j] = 0.;
                    arlequinStab[j][DIM*i+k] = 0.;
                };
                arlequinStabVector[DIM*i+k] = 0.0;
            };
        };
    };

    for (int i = nElNodes; i--; ) delete [] dphi_dx[i];
    delete [] dphi_dx;
    for (int i = DIM; i--; ) delete [] ainv_[i];
    delete [] ainv_;

    return;
};

//------------------------------------------------------------------------------
//----------------------------STEADY LAPLACE PROBEM-----------------------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Element<DIM,DEG>::getLagrangeMultipliersDifferentMesh(int &ielem, double &tPSPG2_, double* press, 
                                                     double* velx, double* vely, double* velxPrev, double* velyPrev,
                                                     double **lagrMultMatrix, double *rhsVectorLM, double *rhsVector){

    double xsi[DIM]= {};
    double xsi_intp[DIM] = {};
    QuadShapeFunction<DIM,DEG> shapeQuad;
    int index = 0;
    SpecialQuad sQuad = SpecialQuad();

    //tARLQ_ = -tPSPG2_;

    double phiLM_[nElNodes] = {};
    double phi_[nElNodes] = {};
    
    double **dphi_dx;
    dphi_dx = new double*[nElNodes];
    for (int i = nElNodes; i--; ) dphi_dx[i] = new double[DIM];
    double **dphiL_dx;
    dphiL_dx = new double*[nElNodes];
    for (int i = nElNodes; i--; ) dphiL_dx[i] = new double[DIM];        
    
    double &dTime_ = parameters.getTimeStep();
    double &visc_ = parameters.getViscosity();
    double &dens_ = parameters.getDensity();
    double &alpha_f = parameters.getAlphaF();
    double &alpha_m = parameters.getAlphaM();
    double &gamma = parameters.getGamma();
    double &k1 = parameters.getArlequinK1();
    double &k2 = parameters.getArlequinK2();
    int &iTimeStep = parameters.getTimeInstant();

    double **ainv_;
    ainv_ = new double*[DIM];
    for (int i = DIM; i--; ) ainv_[i] = new double[DIM];

    // if (index_ == 4936) std::cout << "PSPG Fine " << ielem << " " << tPSPG_ << std::endl;
    for(double* it = sQuad.begin(); it != sQuad.end(); it++){
        
        if ((intPointCorrespElem[index] == ielem)){

            //Defines the integration points adimentional coordinates
            for (int k = 0; k < DIM; k++) xsi[k] = sQuad.PointList(index,k);
            
            //Computes the velocity shape functions
            shapeQuad.evaluate(xsi,phi_);
            
            for (int k = 0; k < DIM; k++) xsi_intp[k] = intPointCorrespXsi[index][k];

            //Computes the coarse mesh shape functions
            shapeQuad.evaluate(xsi_intp,phiLM_);
            
            //Returns the quadrature integration weight
            double weight_ = sQuad.WeightList(index);
            
            double djac_ = 0.;
            //Computes the jacobian matrix
            getJacobianMatrix(xsi_intp, ainv_, djac_);
                        
            getSpatialDerivatives(xsi_intp, ainv_, dphi_dx);

            for (int i = 0; i < nElNodes; ++i)
                for (int j = 0; j < DIM; ++j)
                    dphiL_dx[i][j] = dphi_dx[i][j];

            djac_ = 0.;
            getJacobianMatrix(xsi, ainv_, djac_);
            getSpatialDerivatives(xsi, ainv_, dphi_dx);

            //Lagrange Multiplier
            double lagM_[DIM], una_[DIM] = {}, duna_dx[DIM][DIM] = {};
            interpolateLagMultiplier(index, lagM_);

            //Lagrange Multiplier Derivatives
            double dL_dx[DIM][DIM];
            interpolateLagMultiplierDerivatives(dphi_dx, dL_dx);

            double wna_ = alpha_f * intPointWeightFunctionSpecial[index] + (1. - alpha_f) * intPointWeightFunctionSpecialPrev[index];

            double dalpha_dx = 0.;
            double dalpha_dy = 0.;
            double WJ = weight_ * djac_;

            for (int i = 0; i < nElNodes; ++i){
                // dalpha_dx += ((*nodes_)[connect_[i]] -> getWeightFunction()) * dphi_dx[0][i];
                // dalpha_dy += ((*nodes_)[connect_[i]] -> getWeightFunction()) * dphi_dx[1][i];
                
                una_[0] += alpha_f * velx[i] * phiLM_[i] + (1. - alpha_f) * velxPrev[i] * phiLM_[i];
                una_[1] += alpha_f * vely[i] * phiLM_[i] + (1. - alpha_f) * velyPrev[i] * phiLM_[i];
            }
            for (int i = 0; i < nElNodes; i++){
                for (int j = 0; j < nElNodes; j++){
                    double l2 = phi_[i] * phiLM_[j] * WJ * k1;
                    for (int k = 0; k < DIM; k++){
                        // L2 COUPLING OPERATOR
                        lagrMultMatrix[DIM*i+k][DIM*j+k] += l2;
                        for (int l = 0; l < DIM; l++){
                            //H1 COUPLING OPERATOR
                            double K = dphi_dx[i][l] * dphiL_dx[j][k];
                            if (k==l) for (int m = DIM; m--; ) K += dphi_dx[i][m] * dphiL_dx[j][m];

                            lagrMultMatrix[DIM*i+k][DIM*j+l] += K * WJ * k2;
                        };
                    };
                };
                // Lagrange multipliers residual
                for (int k = 0; k < DIM; k++){
                    double L2 = lagM_[k] * phiLM_[i] * k1;

                    double H1 = 0.;
                    for (int l=DIM; l--; ) H1 += dphiL_dx[i][l] * dL_dx[k][l] * k2;
                    for (int l=DIM; l--; ) H1 += dphiL_dx[i][l] * dL_dx[l][k] * k2;

                    rhsVectorLM[DIM*i+k] -= (L2 + H1) * WJ;

                    double L2u = una_[k] * phi_[i] * k1;

                    double H1u = 0.;
                    for (int l=DIM; l--; ) H1u += dphi_dx[i][l] * duna_dx[k][l] * k2;
                    for (int l=DIM; l--; ) H1u += dphi_dx[i][l] * duna_dx[l][k] * k2;

                    rhsVector[DIM*i+k] -= (L2u + H1u) * WJ;
                };
            };
        };
        index++;        
    };  
    
    for (int i = nElNodes; i--; ) delete [] dphi_dx[i];
    delete [] dphi_dx;
    for (int i = nElNodes; i--; ) delete [] dphiL_dx[i];
    delete [] dphiL_dx;
    for (int i = DIM; i--; ) delete [] ainv_[i];
    delete [] ainv_;

    return;
};


//------------------------------------------------------------------------------
//----------------------------STEADY LAPLACE PROBEM-----------------------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Element<DIM,DEG>::getLagrangeMultipliersSUPG_PSPG_DifferentMesh(int &ielem, double &tPSPG2_, double* press, double* velx, double* vely,
                                                               double **jacobianNRMatrix, double *rhsVector){

    // double xsi[2] = {};
    // double xsi_intp[2] = {};
    // QuadShapeFunction<2,2> shapeQuad;
    // int index = 0;
    // SpecialQuad sQuad = SpecialQuad();

    // //tARLQ_ = -tPSPG2_;

    // double phi_[6] = {};
    // double phiLM_[6] = {};
    
    // double **dphi_dx;
    // dphi_dx = new double*[DIM];
    // for (int i = DIM; i--; ) dphi_dx[i] = new double[nElNodes];
    // double **dphiL_dx;
    // dphiL_dx = new double*[DIM];
    // for (int i = DIM; i--; ) dphiL_dx[i] = new double[nElNodes];        

    // double &dens_ = parameters.getDensity();
    // double &alpha_f = parameters.getAlphaF();
    // double &k1 = parameters.getArlequinK1();
    // double &k2 = parameters.getArlequinK2();
    // int &iTimeStep = parameters.getTimeInstant();

    // // jacobianNRMatrix.clear();
    // // rhsVector.clear();

    // double tSUPG_; double tPSPG_; double tLSIC_;

    // double **ainv_;
    // ainv_ = new double*[DIM];
    // for (int i = DIM; i--; ) ainv_[i] = new double[DIM];

    //    //std::cout << "PSPG Fine " << ielem << " " << tPSPG_ << std::endl;
    // for(double *it = sQuad.begin(); it != sQuad.end(); it++){
        
    //     if ((intPointCorrespElem[index] == ielem)){

    //         //Defines the integration points adimentional coordinates
    //         xsi[0] = sQuad.PointList(index,0);
    //         xsi[1] = sQuad.PointList(index,1);
            
    //         //Computes the velocity shape functions
    //         shapeQuad.evaluate(xsi,phi_);
            
    //         xsi_intp[0] = intPointCorrespXsi[index][0];
    //         xsi_intp[1] = intPointCorrespXsi[index][1];

    //         //Computes the coarse mesh shape functions
    //         shapeQuad.evaluate(xsi_intp,phiLM_);
            
    //         //Returns the quadrature integration weight
    //         double weight_ = sQuad.WeightList(index);
            
    //         //Computes the jacobian matrix
    //         getJacobianMatrix(xsi_intp, ainv_);
                        
    //         getSpatialDerivatives(xsi_intp, ainv_, dphi_dx);

    //         for (int i = 0; i < 2; ++i)
    //             for (int j = 0; j < 6; ++j)
    //                 dphiL_dx[i][j] = dphi_dx[i][j];
    //         // dphiL_dx = dphi_dx;
    //         // ddphiL_dx = ddphi_dx;

    //         getJacobianMatrix(xsi, ainv_);
    //         getSpatialDerivatives(xsi, ainv_, dphi_dx);

    //         getParameterSUPG(tSUPG_, tPSPG_, tLSIC_, phi_, dphi_dx);

    //         double una_ = alpha_f * u_ + (1. - alpha_f) * uPrev_;
    //         double vna_ = alpha_f * v_ + (1. - alpha_f) * vPrev_;

    //         double duna_dx = alpha_f * du_dx + (1. - alpha_f) * duprev_dx;
    //         double duna_dy = alpha_f * du_dy + (1. - alpha_f) * duprev_dy;
    //         double dvna_dx = alpha_f * dv_dx + (1. - alpha_f) * dvprev_dx;
    //         double dvna_dy = alpha_f * dv_dy + (1. - alpha_f) * dvprev_dy;

    //         double lxna_ = lagMx_;
    //         double lyna_ = lagMy_;

    //         double umeshna_ = alpha_f * umesh_ + (1. - alpha_f) * umeshPrev_;
    //         double vmeshna_ = alpha_f * vmesh_ + (1. - alpha_f) * vmeshPrev_;

    //         double wna_ = alpha_f * intPointWeightFunctionSpecial[index] + (1. - alpha_f) * intPointWeightFunctionSpecialPrev[index];

    //         u_ = 0.;
    //         v_ = 0.;

    //         du_dx = 0.;
    //         du_dy = 0.;
    //         dv_dx = 0.;
    //         dv_dy = 0.;

    //         //Interpolates the velocity components and its spatial derivatives
    //         for (int i = 0; i < 6; i++){
    //             u_ += velx[i] * phiLM_[i];
    //             v_ += vely[i] * phiLM_[i];

    //             du_dx += velx[i] * dphiL_dx[0][i];
    //             du_dy += velx[i] * dphiL_dx[1][i];
    //             dv_dx += vely[i] * dphiL_dx[0][i];
    //             dv_dy += vely[i] * dphiL_dx[1][i];
    //         };  
            
    //         for (int i = 0; i < 6; i++){
    //             for (int j = 0; j < 6; j++){
    //                 //SUPG STABILIZATION TERM
    //                 double LM = 0.;
    //                 // LM = ((una_ - umeshna_) * dphi_dx[0][j] + (vna_ - vmeshna_) * dphi_dx[0][j])
    //                 //     * phiLM_(i) * tSUPG_;
                    
    //                 jacobianNRMatrix[2*j  ][2*i  ] += LM * weight_ * djac_;
    //                 jacobianNRMatrix[2*j+1][2*i+1] += LM * weight_ * djac_;
                    
    //                 //PSPG STABILIZATION TERM
    //                 double Lx = 0.; double Ly = 0.;
                    
    //                 // Lx = -dphiL_dx(0,i) * phi_[j] * tPSPG_ / dens_ * k1 * wna_;
    //                 // Ly = -dphiL_dx(1,i) * phi_[j] * tPSPG_ / dens_ * k1 * wna_;
                    
    //                 jacobianNRMatrix[2*j  ][12+i] += Lx * weight_ * djac_;
    //                 jacobianNRMatrix[2*j+1][12+i] += Ly * weight_ * djac_;
    //             };

    //             //SUPG AND PSPG STABILIZATION TERMS
    //             double LMx = 0.; double LMy = 0.; double LMp = 0.;
             
    //             // LMx = ((una_ - umeshna_) * dphiL_dx(0,i) + 
    //             //        (vna_ - vmeshna_) * dphiL_dx(1,i)) * lxna_ * tSUPG_;
    //             // LMy = ((una_ - umeshna_) * dphiL_dx(0,i) + 
    //             //        (vna_ - vmeshna_) * dphiL_dx(1,i)) * lyna_ * tSUPG_;

    //             // LMp = (dphiL_dx(0,i) * lagMx_ + dphiL_dx(1,i) * lagMy_) * tPSPG_/ dens_;
                
    //             // if (iTimeStep > 5)
    //             //     LMp = (dphiL_dx(0,i) * lxna_ + dphiL_dx(1,i) * lyna_) * tPSPG_ / dens_ * k1 +
    //             //           (dphiL_dx(0,i) * (2. * dLx_dxx + dLx_dxy + dLx_dyy) +
    //             //            dphiL_dx(1,i) * (2. * dLy_dyy + dLy_dxy + dLy_dxx)) * tPSPG_ / dens_ * k2;
                
    //             rhsVector[2*i  ] += LMx * weight_ * djac_ * k1;
    //             rhsVector[2*i+1] += LMy * weight_ * djac_ * k1;
    //             rhsVector[12+i] += LMp * weight_ * djac_ * k1;
    //         };
    //     };
    //     index++;        
    // };  
    
    // for (int i = DIM; i--; ) delete [] dphi_dx[i];
    // delete [] dphi_dx;
    // for (int i = DIM; i--; ) delete [] dphiL_dx[i];
    // delete [] dphiL_dx;
    // for (int i = DIM; i--; ) delete [] ainv_[i];
    // delete [] ainv_;

    return;
};


//------------------------------------------------------------------------------
//----------------------------STEADY LAPLACE PROBEM-----------------------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Element<DIM,DEG>::getLagrangeMultipliersArlequinDifferentMesh(int &ielem, double &tPSPG2_,double* press, double* velx, double* vely,
                                                             double **arlequinStab, double **laplMatrix, double *arlequinStabVector){

    double xsi[DIM] = {};
    double xsi_intp[DIM] = {};
    QuadShapeFunction<DIM,DEG> shapeQuad;
    int index = 0;
    SpecialQuad sQuad = SpecialQuad();

    //tARLQ_ = -tPSPG2_;

    double phi_[nElNodes] = {};
    double phiLM_[nElNodes] = {};

    double **dphi_dx;
    dphi_dx = new double*[nElNodes];
    for (int i = nElNodes; i--; ) dphi_dx[i] = new double[DIM];
    double **dphiL_dx;
    dphiL_dx = new double*[nElNodes];
    for (int i = nElNodes; i--; ) dphiL_dx[i] = new double[DIM];        
    
    double &dens_ = parameters.getDensity();
    double &alpha_f = parameters.getAlphaF();
    double &k1 = parameters.getArlequinK1();
    double &k2 = parameters.getArlequinK2();

    // arlequinStab.clear();
    // arlequinStabVector.clear();
    // laplMatrix.clear();

    double **ainv_;
    ainv_ = new double*[DIM];
    for (int i = DIM; i--; ) ainv_[i] = new double[DIM];

    double tSUPG_, tPSPG_, tLSIC_, tARLQ_;

    for(double* it = sQuad.begin(); it != sQuad.end(); it++){
        
        if ((intPointCorrespElem[index] == ielem)){

            //Defines the integration points adimentional coordinates
            for (int k = 0; k < DIM; k++) xsi[k] = sQuad.PointList(index,k);
            
            //Computes the velocity shape functions
            shapeQuad.evaluate(xsi,phi_);
            
            for (int k = 0; k < DIM; k++) xsi_intp[k] = intPointCorrespXsi[index][k];

            //Computes the coarse mesh shape functions
            shapeQuad.evaluate(xsi_intp,phiLM_);
            
            //Returns the quadrature integration weight
            double weight_ = sQuad.WeightList(index);
            
            double djac_ = 0.;
            //Computes the jacobian matrix
            getJacobianMatrix(xsi_intp, ainv_, djac_);
                        
            getSpatialDerivatives(xsi_intp, ainv_, dphi_dx);

            for (int i = 0; i < nElNodes; ++i)
                for (int j = 0; j < DIM; ++j)
                    dphiL_dx[i][j] = dphi_dx[i][j];
            // dphiL_dx = dphi_dx;
            // ddphiL_dx = ddphi_dx;

            djac_ = 0.;
            getJacobianMatrix(xsi, ainv_, djac_);
            getSpatialDerivatives(xsi, ainv_, dphi_dx);

            getParameterArlequin(index, tARLQ_, tSUPG_, tPSPG_, tLSIC_, dphi_dx);


            double wna_ = alpha_f * intPointWeightFunctionSpecial[index] + (1. - alpha_f) * intPointWeightFunctionSpecialPrev[index];

            // u_ = 0.;
            // v_ = 0.;

            // du_dx = 0.;
            // du_dy = 0.;
            // dv_dx = 0.;
            // dv_dy = 0.;

            // //Interpolates the velocity components and its spatial derivatives
            // for (int i = 0; i < 6; i++){
            //     u_ += velx[i] * phiLM_[i];
            //     v_ += vely[i] * phiLM_[i];

            //     du_dx += velx[i] * dphiL_dx[0][i];
            //     du_dy += velx[i] * dphiL_dx[1][i];
            //     dv_dx += vely[i] * dphiL_dx[0][i];
            //     dv_dy += vely[i] * dphiL_dx[1][i];
            // };  

            //Lagrange Multiplier Derivatives
            double dL_dx[DIM][DIM];
            interpolateLagMultiplierDerivatives(dphi_dx, dL_dx);

            for (int i = 0; i < nElNodes; i++){
                for (int j = 0; j < nElNodes; j++){     
                    double AM = 0.;
                    double Lpx = 0.; double Lpy = 0.;
                    double LC = 0.; double LL = 0.;

                    // AM = -phiLM_(i) * phi_[j] * 
                    // intPointWeightFunction(index) * tARLQ_;

                    // LL = phiLM_(i) * phi_[j] * tARLQ_ / dens_;

                    for (int m = DIM; m--; ) LL += dphi_dx[i][m] * dphi_dx[j][m] * tARLQ_ / dens_;

                    for (int k = DIM; k--; )
                        arlequinStab[DIM*i+k][DIM*j+k] += LL * weight_ * djac_;


                    // LL = (dphi_dx[0][i] * dphi_dx[0][j] + dphi_dx[1][i] * dphi_dx[1][j]) * tARLQ_ / dens_;

                    // Lpx = phi_[i] * dphi_dx[0][i] * intPointWeightFunction(index)
                    //     * tARLQ_ / dens_;
                    // Lpy = phi_[i] * dphi_dx[1][i] * intPointWeightFunction(index)
                    //     * tARLQ_ / dens_;

                    // LC = -phi_[j] * ((una_ - umesh_) * dphi_dx[0][j] + (vna_ - vmesh_) * dphi_dx[1][j]) * phiLM_(i) * tARLQ_ * intPointWeightFunction(index);

                    // arlequinStab[2*j  ][2*i  ] += LL * weight_ * djac_;
                    // arlequinStab[2*j+1][2*i+1] += LL * weight_ * djac_;
                    // arlequinStab(12+i,12+j) += 0 * weight_ * djac_;

                    // LC = (dphi_dx[0][i]*(du_dx*dphiL_dx(0,j) + dv_dx*dphiL_dx(1,j)) +
                    //       dphi_dx[1][i]*(du_dy*dphiL_dx(0,j) + dv_dy*dphiL_dx(1,j)) + 
                    //       dphi_dx[0][i]*(u_*ddphiL_dx(0,0)(j) + v_*ddphiL_dx(0,1)(j)) + 
                    //       dphi_dx[1][i]*(u_*ddphiL_dx(1,0)(j) + v_*ddphiL_dx(1,1)(j))) * tARLQ_ * (1-intPointWeightFunction(index));

                    // laplMatrix(2*i  ,2*j  ) += LC * weight_ * djac_;
                    // laplMatrix(2*i+1,2*j+1) += LC * weight_ * djac_;

                     Lpx = 0.; Lpy = 0.;

                    // Lpx = (dphi_dx[0][i] * ddphi_dx(0,0)(j) + 
                    //        dphi_dx[1][i] * ddphi_dx(0,1)(j)) * tARLQ_ * intPointWeightFunction(index);
                    // Lpy = (dphi_dx[0][i] * ddphi_dx(1,0)(j) + 
                    //        dphi_dx[1][i] * ddphi_dx(1,1)(j)) * tARLQ_ * intPointWeightFunction(index);

                    // laplMatrix[2*j  ][2*i  ] += LC * weight_ * djac_;
                    // laplMatrix[2*j+1][2*i+1] += LC * weight_ * djac_;

                    // laplMatrix[2*j  ][12+i] += Lpx * weight_ * djac_;
                    // laplMatrix[2*j+1][12+i] += Lpy * weight_ * djac_;
                    // laplMatrix[12+i][2*j  ] += Lpx * weight_ * djac_;
                    // laplMatrix[12+i][2*j+1] += Lpy * weight_ * djac_;

                };

                //ARLEQUIN STABILIZATION TERMS
                double Amx = 0.; double Amy = 0.;
                double LCx = 0.; double LCy = 0.;
                double LPx = 0.; double LPy = 0.;
                double LLy = 0.;
                

                for (int k = DIM; k--; ){
                    double LLx = 0.;
                    for (int m = DIM; m--; ) LLx -= dphiL_dx[i][m] * dL_dx[k][m]/wna_ * tARLQ_ / dens_;

                    arlequinStabVector[DIM*i+k] += (Amx + LLx) * weight_ * djac_ * wna_;
                }


                // if (iTimeStep > 5){
                //     Amx = -phi_[i] * axm_ * intPointWeightFunctionSpecial(index) * tARLQ_;
                //     Amy = -phi_[i] * aym_ * intPointWeightFunctionSpecial(index) * tARLQ_;
                // }
                //if (iTimeStep < 10){
                    // LLx = -(dphi_dx[0][i] * dLx_dx/(1-wna_) + dphi_dx[1][i] * dLx_dy/(1-wna_)) * tARLQ_;
                    // LLy = -(dphi_dx[0][i] * dLy_dx/(1-wna_) + dphi_dx[1][i] * dLy_dy/(1-wna_)) * tARLQ_;
                //}else{
                    // LLx = -(dphi_dx[0][i] * (dLx_dx/(1-intPointWeightFunction(index)) + duna_dx*duna_dx + dvna_dx*duna_dy) + 
                    //         dphi_dx[1][i] * (dLx_dy/(1-intPointWeightFunction(index)) + duna_dy*duna_dx + dvna_dy*duna_dy)) * tARLQ_ / dens_;
                    // LLy = -(dphi_dx[0][i] * (dLy_dx/(1-intPointWeightFunction(index)) + duna_dx*dvna_dx + dvna_dx*dvna_dy) + 
                    //         dphi_dx[1][i] * (dLy_dy/(1-intPointWeightFunction(index)) + duna_dy*dvna_dx + dvna_dy*dvna_dy)) * tARLQ_ / dens_;
                    // LLx = -(dphi_dx[0][i] * (dLx_dx/(1-intPointWeightFunction(index)) + dp_dxx) + 
                    //         dphi_dx[1][i] * (dLx_dy/(1-intPointWeightFunction(index)) + dp_dxy)) * tARLQ_ / dens_;
                    // LLy = -(dphi_dx[0][i] * (dLy_dx/(1-intPointWeightFunction(index)) + dp_dxy) + 
                    //         dphi_dx[1][i] * (dLy_dy/(1-intPointWeightFunction(index)) + dp_dyy)) * tARLQ_ / dens_;
                //};
     
                // if (iTimeStep > 5){
                //     LLx +=  (dphi_dx[0][i]*dax_dx + dphi_dx[1][i]*dax_dy) * tARLQ_;
                //     LLx +=  (dphi_dx[0][i]*day_dx + dphi_dx[1][i]*day_dy) * tARLQ_ ;
                // };



                // arlequinStabVector[2*i  ] += (Amx + LLx) * weight_ * djac_ * (1-wna_) *0;
                // arlequinStabVector[2*i+1] += (Amy + LLy) * weight_ * djac_ * (1-wna_) *0;
            };
        };
        index++;        
    };  
    
    for (int i = 0; i < nElNodes; i++){
        for (int k = DIM; k--; ){
            if ((*nodes_)[connect_[i]] -> getConstrains(k) == 1) {
                for (int j = 0; j < nLocDOF; j++){
                    arlequinStab[DIM*i+k][j] = 0.;
                    arlequinStab[j][DIM*i+k] = 0.;
                };
                //lagrMultMatrix(12+i,12+i) = 1.;
                arlequinStabVector[DIM*i+k] = 0.0;
                std::cout << "AQUI2 elem different mesh" << std::endl;
            };
        }
    };

    for (int i = nElNodes; i--; ) delete [] dphi_dx[i];
    delete [] dphi_dx;
    for (int i = nElNodes; i--; ) delete [] dphiL_dx[i];
    delete [] dphiL_dx;
    for (int i = DIM; i--; ) delete [] ainv_[i];
    delete [] ainv_;

    return;
};



#endif

