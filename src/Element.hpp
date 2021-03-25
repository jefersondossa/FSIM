
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
    int           *connect_;           //Velocity mesh connectivity 
    int           index_;             //Element index
    double        djac_;              //Jacobian determinant
    double        x_, y_;
    double        u_, v_, p_;     //Interpolated velocity and pressure
    double        uPrev_, vPrev_;
    double        ax_, ay_, axprev_, ayprev_;      //Interpolated acceleration
    double        dax_dx, dax_dy, day_dx, day_dy;      //Interpolated acceleration
    double        umesh_, vmesh_; //Interpolated mesh velocity
    double        umeshPrev_, vmeshPrev_; //Interpolated mesh velocity
    double        du_dx, du_dy, dv_dx, dv_dy;//Interpolated fluid spatial derivatives
    double        duprev_dx, duprev_dy, dvprev_dx, dvprev_dy;
    double        dp_dx, dp_dy;//Interpolated pressure spatial derivative
    double        dumesh_dx, dumesh_dy, dvmesh_dx, dvmesh_dy;//Interpolated mesh velocity derivatives
    double        dumeshPrev_dx, dumeshPrev_dy, dvmeshPrev_dx, dvmeshPrev_dy;
    double        lagMx_,lagMy_;
    double        dLx_dx, dLx_dy, dLy_dx, dLy_dy;
    double        tARLQ_;

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
        connect_ = connect;
        nodes_ = &nodes;
        parameters = param;

        // std::cout <<"AASASD 1"<<std::endl;
        DI = dInt;
        // std::cout <<"AASASD 2"<<std::endl;

        djac_ = 0.; 
        u_ = 0.;          v_ = 0.;                p_ = 0.; 
        umesh_ = 0.;      vmesh_ = 0.;      
        du_dx = 0.;       du_dy = 0.;
        dv_dx = 0.;       dv_dy = 0.;
        dumesh_dx = 0.;       dumesh_dy = 0.;    
        dvmesh_dx = 0.;       dvmesh_dy = 0.;    
        glueZone = false;        FSIInterface = false;
        sideBoundary_ = 0;

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
    /// @param Connectivity element connectivity
    void setConnectivity(int *connect){connect_ = connect;};

    /// Gets the element connectivity
    /// @return element connectivity
    int* getConnectivity(){return connect_;};

    /// Compute and store the spatial jacobian matrix
    /// @param bounded_vector integration point adimensional coordinates
    void getJacobianMatrix(double *xsi, double **ainv_);

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

    /// Interpolate velocity derivatives
    /// @param double** shape function spatial derivatives @param double** velocity derivatives
    /// @param double** previous time step velocity derivatives
    void interpolateVelDerivatives(double **dphi_dx, double du_dx[DIM][DIM], double duprev_dx[DIM][DIM]);


    /// Compute and stores the interpolated velocities, mesh velocities, 
    /// previous mesh velocity, acceleration, among others
    void getVelAndDerivatives(double *phi_, double **dphi_dx);

    /// Compute and store the SUPG, PSPG and LSIC stabilization parameters
    void getParameterSUPG(double &tSUPG_, double &tPSPG_, double &tLSIC_, double *phi_, double **dphi_dx);
    void getParameterArlequin(double &tSUPG_, double &tPSPG_, double &tLSIC_, double *phi_, double **dphi_dx);
    void getParameterArlequin2();

    /// Gets the element jacobian determinant
    /// @return element jacobinan determinant
    double getJacobian(){
        double xsi[DIM] = {};
        double phi_[nElNodes] = {};

        double **ainv_;
        ainv_ = new double*[DIM];
        for (int i = DIM; i--; ) ainv_[i] = new double[DIM];
        
        double **dphi_dx;
        dphi_dx = new double*[DIM];
        for (int i = DIM; i--; ) dphi_dx[i] = new double[nElNodes];

        QuadShapeFunction<DIM,DEG>    shapeQuad;

        xsi[0] = 0.5;
        xsi[1] = 0.5;
                
        shapeQuad.evaluate(xsi,phi_);
        
        //Computes the jacobian matrix
        getJacobianMatrix(xsi, ainv_);
        
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
    void getElemMatrix(int &index, double *phi_, double **dphi_dx, double &tSUPG_, double &tPSPG_, double &tLSIC_, double &weight_, double **jacobianNRMatrix);
    void getElemMatrixInitial(int &index, double *phi_, double **dphi_dx, double &tSUPG_, double &tPSPG_, double &tLSIC_, double &weight_, double **jacobianNRMatrix);

    /// Compute and store the element matrix for the Laplace/Poisson problem
    void getElemLaplMatrix(double &weight_, double **dphi_dx, double** jacobianNRMatrix);

    /// Sets the boundary conditions for the incompressible flow problem
    void setBoundaryConditions(double **jacobianNRMatrix, double *rhsVector);
    void setBoundaryConditions2(double **jacobianNRMatrix, double *rhsVector);

    /// Sets the boundary conditions for the Laplace/Poisson problem
    void setBoundaryConditionsLaplace(double** jacobianNRMatrix, double* rhsVector);

    ///Compute and store the residual vector for the incompressible flow problem
    /// @param int integration point index
    void getResidualVector(int &index, double *phi_, double **dphi_dx, double &tSUPG_, double &tPSPG_, double &tLSIC_, double &weight_, double *rhsVector);
    void getResidualVectorInitial(int &index, double *phi_, double **dphi_dx, double &tSUPG_, double &tPSPG_, double &tLSIC_, double &weight_, double *rhsVector);

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
    void getTransientNavierStokesInitial(double *rhsVector, double **jacobianNRMatrix);

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
//----------------------SET ELEMENT INTERSECTION PARAMETERS---------------------
//------------------------------------------------------------------------------
template<>
void Element<2,2>::setIntersectionParameters(double *x, double *X) {
    xK[0] = x[0]; xK[1] = x[1]; 
    XK[0] = X[0]; XK[1] = X[1];
    return;
};

//------------------------------------------------------------------------------
//----------------------SET ELEMENT INTERSECTION PARAMETERS---------------------
//------------------------------------------------------------------------------
template<>
void Element<2,2>::setIntegPointWeightFunction() {
    
    double xsi[2] = {};
    ShapeFunction shapeQuad;
    double phi_[6] = {};
    
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

       for (int j=0; j<6; j++){
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

       for (int j=0; j<6; j++){
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
template<>
void Element<2,2>::getIntegPointCoordinates(){

    SpecialQuad sQuad = SpecialQuad();
    double xsi[2] = {};
    ShapeFunction   shapeQuad;
    double phi_[6] = {};

    for (int i = 0; i < sQuad.end() - sQuad.begin(); i++){
        double x[2] = {};

        xsi[0] = sQuad.PointList(i,0);
        xsi[1] = sQuad.PointList(i,1);

        shapeQuad.evaluate(xsi,phi_);

        intPointCoordinates[i][0] = 0.;
        intPointCoordinates[i][1] = 0.;

        for (int j = 0; j < 6; j++){
            intPointCoordinates[i][0] += (*nodes_)[connect_[j]] -> getCoordinateValue(0) * phi_[j];
            intPointCoordinates[i][1] += (*nodes_)[connect_[j]] -> getCoordinateValue(1) * phi_[j];
        }
    };


    return;
};

//------------------------------------------------------------------------------
//---------------------------CLEAR ELEMENT VARIABLES----------------------------
//------------------------------------------------------------------------------
template<>
void Element<2,2>::clearVariables(){

    djac_ = 0.;
        
    u_ = 0.;          v_ = 0.;        p_ = 0.; 
    umesh_ = 0.;      vmesh_ = 0.; 
    du_dx = 0.;       du_dy = 0.;
    dv_dx = 0.;       dv_dy = 0.;
    
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
void Element<DIM,DEG>::getJacobianMatrix(double *xsi, double **ainv_) {

    //Computes the spatial Jacobian matrix and its inverse
    double **dphi;
    dphi = new double*[DIM];
    for (int i = DIM; i--; ) dphi[i] = new double[nElNodes];
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
                dx_dxsi[j][k] += xna_[j] * dphi[k][i];
            };
        };
    };

    //Computing the jacobian determinant and Inverse
    computeInverseAndDeterminant(djac_, ainv_, dx_dxsi);

    for (int i = DIM; i--; ) delete [] dphi[i];
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
    dphi = new double*[DIM];
    for (int i = DIM; i--; ) dphi[i] = new double[nElNodes];

    ShapeFunction shapeQuad;
    
    shapeQuad.evaluateGradient(xsi,dphi);
    // shapeQuad.evaluateHessian(xsi,ddphi);
    
    for (int i = nElNodes; i--; )
        for (int j = DIM; j--; )
            dphi_dx[j][i] = 0.;

    //Quadratic shape functions spatial first derivatives
    // noalias(dphi_dx) = prod(ainv_,dphi);
    for (int i = nElNodes; i--; )
        for (int k = DIM; k--; )
            for (int j = DIM; j--; )
                dphi_dx[k][i] += ainv_[k][j] * dphi[j][i];

    for (int i = DIM; i--; ) delete [] dphi[i];
    delete [] dphi;

    return;
};

//------------------------------------------------------------------------------
//-------------INTERPOLATES VELOCITY, PRESSURE AND ITS DERIVATIVES--------------
//------------------------------------------------------------------------------
template<>
void Element<2,2>::getVelAndDerivatives(double *phi_, double **dphi_dx) {


    u_ = 0.; v_ = 0.;
    x_ = 0.; y_ = 0.;
    ax_ = 0.; ay_ = 0.; axprev_ = 0.; ayprev_ = 0.;
    uPrev_ = 0.; vPrev_ = 0.;
    umesh_ = 0.; vmesh_ = 0.;
    p_ = 0.;
    du_dx = 0.; du_dy = 0.; dv_dx = 0.; dv_dy = 0.;
    dp_dx = 0.; dp_dy = 0.;
    lagMx_ = 0.; lagMy_ = 0.;
    dLx_dx = 0.; dLx_dy = 0.; dLy_dx = 0.; dLy_dy = 0.;
    duprev_dx = 0.; duprev_dy = 0.; dvprev_dx = 0.; dvprev_dy = 0.;
    dax_dx = 0.; dax_dy = 0.; day_dx = 0.; day_dy = 0.;
    dumesh_dx = 0.; dumesh_dy = 0.; dvmesh_dx = 0.; dvmesh_dy = 0.;
    umeshPrev_ = 0.;    vmeshPrev_ = 0.;
    dumeshPrev_dx = 0.; dumeshPrev_dy = 0.; dvmeshPrev_dx = 0.; dvmeshPrev_dy = 0.;

    double &alpha_f = parameters.getAlphaF();
    double &alpha_m = parameters.getAlphaM();

    //Interpolates the velocity components and its spatial derivatives
    for (int i = 0; i < 6; i++){
        u_ += (*nodes_)[connect_[i]] -> getVelocity(0) * phi_[i];
        v_ += (*nodes_)[connect_[i]] -> getVelocity(1) * phi_[i];

        double xna_ = alpha_f * (*nodes_)[connect_[i]] -> getCoordinateValue(0) + 
                      (1. - alpha_f) * (*nodes_)[connect_[i]] -> getPreviousCoordinateValue(0);
        double yna_ = alpha_f * (*nodes_)[connect_[i]] -> getCoordinateValue(1) + 
                      (1. - alpha_f) * (*nodes_)[connect_[i]] -> getPreviousCoordinateValue(1);
        x_ += xna_ * phi_[i];
        y_ += yna_ * phi_[i];

        uPrev_ += (*nodes_)[connect_[i]] -> getPreviousVelocity(0) * phi_[i];
        vPrev_ += (*nodes_)[connect_[i]] -> getPreviousVelocity(1) * phi_[i];

        ax_ += (*nodes_)[connect_[i]] -> getAcceleration(0) * phi_[i];
        ay_ += (*nodes_)[connect_[i]] -> getAcceleration(1) * phi_[i];
        
        axprev_ += (*nodes_)[connect_[i]] -> getPreviousAcceleration(0) * phi_[i];
        ayprev_ += (*nodes_)[connect_[i]] -> getPreviousAcceleration(1) * phi_[i];

        dax_dx += ((*nodes_)[connect_[i]] -> getAcceleration(0) + (1.-alpha_m) * (*nodes_)[connect_[i]] -> getPreviousAcceleration(0)) * dphi_dx[0][i];
        dax_dy += ((*nodes_)[connect_[i]] -> getAcceleration(0) + (1.-alpha_m) * (*nodes_)[connect_[i]] -> getPreviousAcceleration(0)) * dphi_dx[1][i];
        day_dx += ((*nodes_)[connect_[i]] -> getAcceleration(1) + (1.-alpha_m) * (*nodes_)[connect_[i]] -> getPreviousAcceleration(1)) * dphi_dx[0][i];
        day_dy += ((*nodes_)[connect_[i]] -> getAcceleration(1) + (1.-alpha_m) * (*nodes_)[connect_[i]] -> getPreviousAcceleration(1)) * dphi_dx[1][i];

        du_dx += (*nodes_)[connect_[i]] -> getVelocity(0) * dphi_dx[0][i];
        du_dy += (*nodes_)[connect_[i]] -> getVelocity(0) * dphi_dx[1][i];
        dv_dx += (*nodes_)[connect_[i]] -> getVelocity(1) * dphi_dx[0][i];
        dv_dy += (*nodes_)[connect_[i]] -> getVelocity(1) * dphi_dx[1][i];

        duprev_dx += (*nodes_)[connect_[i]] -> getPreviousVelocity(0) * dphi_dx[0][i];
        duprev_dy += (*nodes_)[connect_[i]] -> getPreviousVelocity(0) * dphi_dx[1][i];
        dvprev_dx += (*nodes_)[connect_[i]] -> getPreviousVelocity(1) * dphi_dx[0][i];
        dvprev_dy += (*nodes_)[connect_[i]] -> getPreviousVelocity(1) * dphi_dx[1][i];

        // du_dxx += (*nodes_)[connect_[i]] -> getVelocity(0) * ddphi_dx(0,0)(i);
        // du_dyy += (*nodes_)[connect_[i]] -> getVelocity(0) * ddphi_dx(1,1)(i);
        // du_dxy += (*nodes_)[connect_[i]] -> getVelocity(0) * ddphi_dx(0,1)(i);
        // dv_dxx += (*nodes_)[connect_[i]] -> getVelocity(1) * ddphi_dx(0,0)(i);
        // dv_dyy += (*nodes_)[connect_[i]] -> getVelocity(1) * ddphi_dx(1,1)(i);
        // dv_dxy += (*nodes_)[connect_[i]] -> getVelocity(1) * ddphi_dx(0,1)(i);

        p_ += (*nodes_)[connect_[i]] -> getPressure() * phi_[i];
        
        dp_dx += (*nodes_)[connect_[i]] -> getPressure() * dphi_dx[0][i];
        dp_dy += (*nodes_)[connect_[i]] -> getPressure() * dphi_dx[1][i];

        umesh_ += (*nodes_)[connect_[i]] -> getMeshVelocity(0) * phi_[i];
        vmesh_ += (*nodes_)[connect_[i]] -> getMeshVelocity(1) * phi_[i];

        umeshPrev_ += (*nodes_)[connect_[i]] -> getPreviousMeshVelocity(0) * phi_[i];
        vmeshPrev_ += (*nodes_)[connect_[i]] -> getPreviousMeshVelocity(1) * phi_[i];
        
        dumesh_dx += (*nodes_)[connect_[i]] -> getMeshVelocity(0) * dphi_dx[0][i];
        dumesh_dy += (*nodes_)[connect_[i]] -> getMeshVelocity(0) * dphi_dx[1][i];
        dvmesh_dx += (*nodes_)[connect_[i]] -> getMeshVelocity(1) * dphi_dx[0][i];
        dvmesh_dy += (*nodes_)[connect_[i]] -> getMeshVelocity(1) * dphi_dx[1][i];

        dumeshPrev_dx += (*nodes_)[connect_[i]] -> getPreviousMeshVelocity(0) * dphi_dx[0][i];
        dumeshPrev_dy += (*nodes_)[connect_[i]] -> getPreviousMeshVelocity(0) * dphi_dx[1][i];
        dvmeshPrev_dx += (*nodes_)[connect_[i]] -> getPreviousMeshVelocity(1) * dphi_dx[0][i];
        dvmeshPrev_dy += (*nodes_)[connect_[i]] -> getPreviousMeshVelocity(1) * dphi_dx[1][i];
        // std::cout << "velocity " << phi_[i] << " " << (*nodes)[connect_[i]] -> getVelocity(0) << std::endl;

        lagMx_ += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(0) * phi_[i];
        lagMy_ += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(1) * phi_[i];

        // dLx_dxx += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(0) * ddphi_dx(0,0)(i);
        // dLx_dyy += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(0) * ddphi_dx(1,1)(i);
        // dLx_dxy += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(0) * ddphi_dx(0,1)(i);
        // dLy_dxx += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(1) * ddphi_dx(0,0)(i);
        // dLy_dyy += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(1) * ddphi_dx(1,1)(i);
        // dLy_dxy += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(1) * ddphi_dx(0,1)(i);        

        dLx_dx += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(0) * dphi_dx[0][i];
        dLx_dy += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(0) * dphi_dx[1][i];
        dLy_dx += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(1) * dphi_dx[0][i];
        dLy_dy += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(1) * dphi_dx[1][i];
    };  





    //std::cout << "VEL " << umesh_ << " " << vmesh_ << std::endl;

    // ax_ = (u_ - uPrev_) / dTime_;
    // ay_ = (v_ - vPrev_) / dTime_;

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
        double shapeFi = DI->phi_[i][index];

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
        double shapeFi = DI->phi_[i][index];
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
        double shapeFi = DI->phi_[i][index];
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
void Element<DIM,DEG>::interpolateVelocity(int &index, double *u_, double *uPrev_) {

    for (int i = DIM; i--; ) {
        u_[i] = 0.;
        uPrev_[i] = 0.;
    }

    for (int i = nElNodes; i--; ){
        double shapeFi = DI->phi_[i][index];
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
//-------------INTERPOLATES VELOCITY, PRESSURE AND ITS DERIVATIVES--------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Element<DIM,DEG>::getBoundaryLoad(double* xsi, double* load) {

    double &visc_ = parameters.getViscosity();

    double phi_[nElNodes] = {};

    double **dphi_dx;
    dphi_dx = new double*[DIM];
    for (int i = DIM; i--; ) dphi_dx[i] = new double[nElNodes];

    double shearStress[DIM][DIM] = {};
    
    double **dphi;
    dphi = new double*[DIM];
    for (int i = DIM; i--; ) dphi[i] = new double[nElNodes];
    double **ainv_;
    ainv_ = new double*[DIM];
    for (int i = DIM; i--; ) ainv_[i] = new double[DIM];

    double t_vector[2] = {};
    double n_vector[2] = {};
    double ident[2][2] = {}; ident[0][0] = 1.; ident[1][1] = 1.;
    ShapeFunction shapeQuad;
    

    //Computes the shape functions        
    shapeQuad.evaluate(xsi,phi_);

    //Computes the jacobian matrix
    getJacobianMatrix(xsi, ainv_);
    

    //Computes spatial derivatives
    getSpatialDerivatives(xsi, ainv_, dphi_dx);

    shapeQuad.evaluateGradient(xsi, dphi); 

    //Interpolates velocity and its derivatives values
    getVelAndDerivatives(phi_, dphi_dx);
   
    shearStress[0][0] = 2. * visc_ * du_dx;
    shearStress[0][1] = visc_ * (du_dy + dv_dx);
    shearStress[1][0] = visc_ * (du_dy + dv_dx);
    shearStress[1][1] = 2. * visc_ * dv_dy;

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

            dx_dxsiB[j][0] += xna_[j] * dphi[j][aux[i]];    
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
    

    for (int i = DIM; i--; ) delete [] dphi_dx[i];
    delete [] dphi_dx;
    for (int i = DIM; i--; ) delete [] dphi[i];
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
    double phi_[6] = {};
    
    double **dphi_dx;
    dphi_dx = new double*[DIM];
    for (int i = DIM; i--; ) dphi_dx[i] = new double[nElNodes];



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
        
        //Computes the jacobian matrix
        getJacobianMatrix(xsi, ainv_);

        //Computes spatial derivatives
        getSpatialDerivatives(xsi, ainv_, dphi_dx);

        //Interpolates velocity and its derivatives values
        getVelAndDerivatives(phi_, dphi_dx);        

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

        shearStress[0][0] = 2. * visc_ * du_dx;
        shearStress[0][1] = visc_ * (du_dy + dv_dx);
        shearStress[1][0] = visc_ * (du_dy + dv_dx);
        shearStress[1][1] = 2. * visc_ * dv_dy;

        for (int i = 0; i < 2; i++){
            for (int j = 0; j < 2; j++){
                load_pressure[i] += -p_ * ident[i][j] * n_vector[j] * jacb_ * weightB;
                load_friction[i] += shearStress[i][j] * n_vector[j] * jacb_ * weightB;
            }
        }

        moment += ((-p_ + shearStress[0][0] + shearStress[1][0]) * y_ 
                 -(-p_ + shearStress[0][1] + shearStress[1][1]) * (x_ - 0.248792267683901))
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

    for (int i = DIM; i--; ) delete [] dphi_dx[i];
    delete [] dphi_dx;
    for (int i = DIM; i--; ) delete [] ainv_[i];
    delete [] ainv_;

    return;
};




//------------------------------------------------------------------------------
//------------------COMPUTES THE SUPG STABILIZATION PARAMETER-------------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Element<DIM,DEG>::getParameterSUPG(double &tSUPG_, double &tPSPG_, double &tLSIC_, double *phi_, double **dphi_dx) {

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
            u__[j] += ua * phi_[i];
            a1 += ua*ua;
        }

        // a1 = std::sqrt(a1);
        for (int j = DIM; j--; ) r[j] += std::sqrt(a1) * dphi_dx[j][i];
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
            hRGN_ += r[j] * dphi_dx[j][i];
            hUGN_ += s[j] * dphi_dx[j][i];
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
template<>
void Element<2,2>::getParameterArlequin(double &tSUPG_, double &tPSPG_, double &tLSIC_, double *phi_, double **dphi_dx) {

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

    for (int i = 0; i < 6; i++){
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
        
        u__ += ua * phi_[i];
        v__ += va * phi_[i];

        lx__ += lxa * phi_[i];
        ly__ += lya * phi_[i];
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


    for (int i = 0; i < 6; i++){
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
        
        r[0] += std::sqrt(ua * ua + va * va) * dphi_dx[0][i];
        r[1] += std::sqrt(ua * ua + va * va) * dphi_dx[1][i];

        rl[0] += std::sqrt(lxa * lxa + lya * lya) * dphi_dx[0][i];
        rl[1] += std::sqrt(lxa * lxa + lya * lya) * dphi_dx[1][i];
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
    for (int i = 0; i < 6; i++){
        hRGN_ += std::fabs(r[0] * dphi_dx[0][i] + r[1] * dphi_dx[1][i]);
        hUGN_ += std::fabs(s[0] * dphi_dx[0][i] + s[1] * dphi_dx[1][i]);        

        hRGNL_ += std::fabs(rl[0] * dphi_dx[0][i] + rl[1] * dphi_dx[1][i]);
        hUGNL_ += std::fabs(sl[0] * dphi_dx[0][i] + sl[1] * dphi_dx[1][i]);        
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
    tARLQ_ = -1. * k1 * tSUPG_ * 1.e0;
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

//------------------------------------------------------------------------------
//------------------COMPUTES THE SUPG STABILIZATION PARAMETER-------------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Element<DIM,DEG>::getParameterArlequin2() {

    double xsi[2] = {};
    double phi_[6] = {};
    
    double **dphi_dx;
    dphi_dx = new double*[DIM];
    for (int i = DIM; i--; ) dphi_dx[i] = new double[nElNodes];


    ShapeFunction           shapeQuad;
    double **ainv_;
    ainv_ = new double*[DIM];
    for (int i = DIM; i--; ) ainv_[i] = new double[DIM];
    NormalQuad nQuad = NormalQuad();

    double &visc_ = parameters.getViscosity();
    double &dens_ = parameters.getDensity();
    double &dTime_ = parameters.getTimeStep();
    double &k1 = parameters.getArlequinK1();

    double lambda[18][18] = {};
    double inercia[18][18] = {};
    double conveccao[18][18] = {};
    double visccc[18][18] = {};

    double le = 0.;
    double li = 0.;
    double lc = 0.;
    double lv = 0.;

    tARLQ_ = 0.;
    
    int index = 0;

    for(double* it = nQuad.begin(); it != nQuad.end(); it++){

        //Defines the integration points adimentional coordinates
        xsi[0] = nQuad.PointList(index,0);
        xsi[1] = nQuad.PointList(index,1);

        //Computes the velocity shape functions
        shapeQuad.evaluate(xsi,phi_);

        //Returns the quadrature integration weight
        double weight_ = nQuad.WeightList(index);

        //Computes the jacobian matrix
        getJacobianMatrix(xsi, ainv_);

        getSpatialDerivatives(xsi, ainv_, dphi_dx);
        
        getVelAndDerivatives(phi_, dphi_dx);

        for (int i = 0; i < 6; i++){
            for (int j = 0; j < 6; j++){
                
                lambda[2*i  ][2*j  ] += phi_[i] * phi_[j] * weight_ * djac_;
                lambda[2*i+1][2*j+1] += phi_[i] * phi_[j] * weight_ * djac_;

                conveccao[2*i  ][2*j  ] += dphi_dx[0][i] * (du_dx*dphi_dx[0][j] + dv_dx*dphi_dx[1][j]) * weight_ * djac_;
                conveccao[2*i  ][2*j+1] += dphi_dx[1][i] * (du_dx*dphi_dx[0][j] + dv_dx*dphi_dx[1][j]) * weight_ * djac_;
                conveccao[2*i+1][2*j  ] += dphi_dx[0][i] * (du_dy*dphi_dx[0][j] + dv_dy*dphi_dx[1][j]) * weight_ * djac_;
                conveccao[2*i+1][2*j+1] += dphi_dx[1][i] * (du_dy*dphi_dx[0][j] + dv_dy*dphi_dx[1][j]) * weight_ * djac_;

                inercia[2*i  ][2*j  ] += dphi_dx[0][i] * dphi_dx[0][j] * weight_ * djac_;
                inercia[2*i  ][2*j+1] += dphi_dx[1][i] * dphi_dx[1][j] * weight_ * djac_;
                inercia[2*i+1][2*j  ] += dphi_dx[0][i] * dphi_dx[0][j] * weight_ * djac_;
                inercia[2*i+1][2*j+1] += dphi_dx[1][i] * dphi_dx[1][j] * weight_ * djac_;

                // visccc(2*i  ,2*j  ) += ((ddphi_dx(0,0)(i)+ddphi_dx(0,1)(i)) * (2. * ddphi_dx(0,0)(j) + ddphi_dx(1,1)(j))) * weight_ * djac_;
                // visccc(2*i  ,2*j+1) += ((ddphi_dx(0,0)(i)+ddphi_dx(0,1)(i)) * (ddphi_dx(0,1)(j))) * weight_ * djac_;
                // visccc(2*i+1,2*j  ) += ((ddphi_dx(1,1)(i)+ddphi_dx(0,1)(i)) * (ddphi_dx(0,1)(j))) * weight_ * djac_;
                // visccc(2*i+1,2*j+1) += ((ddphi_dx(1,1)(i)+ddphi_dx(0,1)(i)) * (2. * ddphi_dx(1,1)(j) + ddphi_dx(0,0)(j))) * weight_ * djac_;
            }
        }
        index++;
    }

    double U_[18] = {};
    double DU_[18] = {};

    for (int i=0; i<6; i++){
        U_[2*i  ] = (*nodes_)[connect_[i]] -> getVelocity(0);
        U_[2*i+1] = (*nodes_)[connect_[i]] -> getVelocity(1);
        DU_[2*i  ] = (*nodes_)[connect_[i]] -> getAcceleration(0);
        DU_[2*i+1] = (*nodes_)[connect_[i]] -> getAcceleration(1);
    };

    double a1[18] = {};
    double a2[18] = {};
    double a3[18] = {};
    double a4[18] = {};

    for (int i = 0; i < 18; i++){
        for (int j = 0; j < 18; j++){
            a1[i] += lambda[i][j] * U_[j];
            a2[i] += conveccao[i][j] * U_[j];
            a3[i] += inercia[i][j] * DU_[j];
            a4[i] += visccc[i][j] * U_[j];
        }
    }

    le = 0.; //norm_2(prod(lambda,U_));
    lc = 0.; //norm_2(prod(conveccao,U_));
    li = 0.; //norm_2(prod(inercia,DU_));
    lv = 0.; //norm_2(prod(visccc,U_));
    
    for (int i = 0; i < 18; i++){
        le += a1[i] * a1[i];
        lc += a2[i] * a2[i];
        li += a3[i] * a3[i];
        lv += a4[i] * a4[i];
    }


    // //std::cout << "asd " << index_ << " " << le << " " << lc << " " << li << " " << tARLQ_ << std::endl;
    //  // std::cout << "antes " << index_ << " " << tARLQ_ << std::endl;

    double tp1 = std::fabs(le) / std::fabs(lc);
    double tp2 = std::fabs(le) / std::fabs(li);
    double tp3 = tp1 / visc_;
    tp3 = le / (visc_*lv);

    if (std::fabs(tp1) <= 1.e-10) tp1 = 1.e-10;
    if (std::fabs(tp2) <= 1.e-10) tp2 = 1.e-10;
    if (std::fabs(tp3) <= 1.e-10) tp3 = 1.e-10;

    tARLQ_ = -1000.0 * k1 / std::sqrt(1. / (tp1*tp1) +
                            1. / (tp2*tp2) +
                            1. / (tp3*tp3));

    // std::cout << "depois " << index_ << " " << le << " " << lv << " " << tARLQ_ << std::endl;

    for (int i = DIM; i--; ) delete [] dphi_dx[i];
    delete [] dphi_dx;
    for (int i = DIM; i--; ) delete [] ainv_[i];
    delete [] ainv_;

    return;
};

//------------------------------------------------------------------------------
//----------------------ELEMENT DIFFUSION/VISCOSITY MATRIX----------------------
//------------------------------------------------------------------------------
template<>
void Element<2,2>::getElemMatrix(int &index, double *phi_, double **dphi_dx, double &tSUPG_, double &tPSPG_, double &tLSIC_, double &weight_, double **jacobianNRMatrix){

    double &dTime_ = parameters.getTimeStep();
    double &visc_ = parameters.getViscosity();
    double &dens_ = parameters.getDensity();
    double &alpha_f = parameters.getAlphaF();
    double &alpha_m = parameters.getAlphaM();
    double &gamma = parameters.getGamma();
    int &iTimeStep = parameters.getTimeInstant();
    double &pi = parameters.getPi();


    double una_ = alpha_f * u_ + (1. - alpha_f) * uPrev_;
    double vna_ = alpha_f * v_ + (1. - alpha_f) * vPrev_;

    double duna_dx = alpha_f * du_dx + (1. - alpha_f) * duprev_dx;
    double duna_dy = alpha_f * du_dy + (1. - alpha_f) * duprev_dy;
    double dvna_dx = alpha_f * dv_dx + (1. - alpha_f) * dvprev_dx;
    double dvna_dy = alpha_f * dv_dy + (1. - alpha_f) * dvprev_dy;

    double axm_ = alpha_m * ax_ + (1. - alpha_m) * axprev_;
    double aym_ = alpha_m * ay_ + (1. - alpha_m) * ayprev_;

    double umeshna_ = alpha_f * umesh_ + (1. - alpha_f) * umeshPrev_;
    double vmeshna_ = alpha_f * vmesh_ + (1. - alpha_f) * vmeshPrev_;

    double wna_ = alpha_f * intPointWeightFunction[index] + (1. - alpha_f) * intPointWeightFunctionPrev[index];

    // if ((wna_ < 0.99) && (model == false) && (index_ == 2144)) std::cout << "IP weight " << " " << index_ << " " << wna_ << " " << intPointWeightFunction(index) << " " << intPointWeightFunctionPrev(index) << std::endl;
     
    for (int i = 0; i < 6; i++){
        for (int j = 0; j < 6; j++){

            double wSUPGi = (una_ - umeshna_) * dphi_dx[0][i] + (vna_ - vmeshna_) * dphi_dx[1][i];
            double wSUPGj = (una_ - umeshna_) * dphi_dx[0][j] + (vna_ - vmeshna_) * dphi_dx[1][j];
            
            //Mass matrix (use for both directions)
            double mM = phi_[i] * phi_[j] * dens_ * alpha_m + 
                        wSUPGi * phi_[j] * tSUPG_ * dens_ * alpha_m;
                       
            //Difusion matrix (viscosity)
            double Kxx = (2. * dphi_dx[0][i] * dphi_dx[0][j] + dphi_dx[1][i] * dphi_dx[1][j]) 
                        * visc_ * alpha_f * gamma * dTime_;
            double Kxy = dphi_dx[1][i] * dphi_dx[0][j] 
                        * visc_ * alpha_f * gamma * dTime_;
            double Kyx = dphi_dx[0][i] * dphi_dx[1][j] 
                        * visc_ * alpha_f * gamma * dTime_;
            double Kyy = (2. * dphi_dx[1][i] * dphi_dx[1][j] + dphi_dx[0][i] * dphi_dx[0][j]) 
                        * visc_ * alpha_f * gamma * dTime_;
            
            //Convection matrix
            double Cxx = (dphi_dx[0][j] * (una_ - umeshna_) + dphi_dx[1][j] * (vna_ - vmeshna_)) * phi_[i] 
                        * dens_ * alpha_f * gamma * dTime_
                        + wSUPGi * wSUPGj * tSUPG_ * dens_ * alpha_f * gamma * dTime_;
            double Cyy  = Cxx;
            
            double Cuu = (phi_[i] * duna_dx * phi_[j] +
                          tSUPG_ * wSUPGi * duna_dx * phi_[j] + 
                          tSUPG_ * dphi_dx[0][i] * phi_[j] * 
                          (una_ * duna_dx + vna_ * duna_dy)) * dens_ * alpha_f * gamma * dTime_;
            
            double Cuv = (phi_[i] * duna_dy * phi_[j] +
                          tSUPG_ * wSUPGi * duna_dy * phi_[j] + 
                          tSUPG_ * dphi_dx[1][i] * phi_[j] * 
                          (una_ * duna_dx + vna_ * duna_dy)) * dens_ * alpha_f * gamma * dTime_;
            
            double Cvu = (phi_[i] * dvna_dx * phi_[j] +
                          tSUPG_ * wSUPGi * dvna_dx * phi_[j] + 
                          tSUPG_ * dphi_dx[0][i] * phi_[j] * 
                          (una_ * dvna_dx + vna_ * dvna_dy)) * dens_ * alpha_f * gamma * dTime_;
            
            double Cvv = (phi_[i] * dvna_dy * phi_[j] +
                          tSUPG_ * wSUPGi * dvna_dy * phi_[j] + 
                          tSUPG_ * dphi_dx[1][i] * phi_[j] * 
                          (una_ * dvna_dx + vna_ * dvna_dy)) * dens_ * alpha_f * gamma * dTime_;
            
            double KLSxx = dphi_dx[0][i] * dphi_dx[0][j] * tLSIC_ * dens_ * alpha_f * gamma * dTime_;
            double KLSxy = dphi_dx[0][i] * dphi_dx[1][j] * tLSIC_ * dens_ * alpha_f * gamma * dTime_;
            double KLSyx = dphi_dx[1][i] * dphi_dx[0][j] * tLSIC_ * dens_ * alpha_f * gamma * dTime_;
            double KLSyy = dphi_dx[1][i] * dphi_dx[1][j] * tLSIC_ * dens_ * alpha_f * gamma * dTime_;

            jacobianNRMatrix[2*i  ][2*j  ] += (mM + Kxx + KLSxx + Cxx + Cuu) * weight_ * djac_ * wna_;
            jacobianNRMatrix[2*i+1][2*j+1] += (mM + Kyy + KLSyy + Cyy + Cvv) * weight_ * djac_ * wna_;
            jacobianNRMatrix[2*i  ][2*j+1] += (Kxy + Cuv + KLSxy) * weight_ * djac_ * wna_;
            jacobianNRMatrix[2*i+1][2*j  ] += (Kyx + Cvu + KLSyx) * weight_ * djac_ * wna_; 

            
            //SINAL DA PARCELA QUE MULTIPLICA O TSUPG ESTA
            //COM SINAL TROCADO NA FORMULAÇAO DO TEZDUYAR
            //multipy pressure direction x
            double QSUPGx = - (dphi_dx[0][i] * phi_[j]) + 
                             ((dphi_dx[0][i] * (una_ - umeshna_) + 
                               dphi_dx[1][i] * (vna_ - vmeshna_)) * dphi_dx[0][j] * tSUPG_);
            //multiply pressure direction y
            double QSUPGy = - (dphi_dx[1][i] * phi_[j]) +
                             ((dphi_dx[0][i] * (una_ - umeshna_) + 
                               dphi_dx[1][i] * (vna_ - vmeshna_)) * dphi_dx[1][j] * tSUPG_);
            //multiply velocity direction x
            double Qx = dphi_dx[0][i] * phi_[j] * alpha_f * gamma * dTime_;
            //multiply velocity direction y
            double Qy = dphi_dx[1][i] * phi_[j] * alpha_f * gamma * dTime_;
            
            jacobianNRMatrix[12+j][2*i  ] += Qx * weight_ * djac_ * wna_;
            jacobianNRMatrix[12+j][2*i+1] += Qy * weight_ * djac_ * wna_;
            jacobianNRMatrix[2*i  ][12+j] += QSUPGx * weight_ * djac_ * wna_;
            jacobianNRMatrix[2*i+1][12+j] += QSUPGy * weight_ * djac_ * wna_;


            double Hx = dphi_dx[0][i] * phi_[j] * tPSPG_ * alpha_m;
            double Hy = dphi_dx[1][i] * phi_[j] * tPSPG_ * alpha_m;
            
            double Gx = dphi_dx[0][i] * wSUPGj * tPSPG_ * alpha_f * gamma * dTime_;
            double Gy = dphi_dx[1][i] * wSUPGj * tPSPG_ * alpha_f * gamma * dTime_;

            double Guu = (dphi_dx[0][i] * duna_dx * phi_[j] + 
                          dphi_dx[1][i] * dvna_dx * phi_[j]) * tPSPG_ * alpha_f * gamma * dTime_ * 0.;
            double Gvv = (dphi_dx[0][i] * duna_dy * phi_[j] +  
                          dphi_dx[1][i] * dvna_dy * phi_[j]) * tPSPG_ * alpha_f * gamma * dTime_ * 0.;
            
            double Q = (dphi_dx[0][i] * dphi_dx[0][j] + 
                        dphi_dx[1][i] * dphi_dx[1][j]) * tPSPG_ / (dens_);

            jacobianNRMatrix[12+j][2*i  ] += (Hx + Gx + Guu) * weight_ * djac_ * wna_;
            jacobianNRMatrix[12+j][2*i+1] += (Hy + Gy + Gvv) * weight_ * djac_ * wna_;
            jacobianNRMatrix[12+j][12+i] += Q * weight_ * djac_ * wna_;
        };
    };

    //std::cout << intPointWeightFunction(index) << std::endl;

    if (model){
        //        std::cout << " asas " << dens_ << " " << visc_ << " " << intPointWeightFunction(index) << std::endl;
    // std::cout << "matrix " << std::endl;
    // for (int i=0; i<18; i++){
    //     for (int j=0; j<18; j++){        
    //         std::cout << jacobianNRMatrix(i,j)<< " " ;
    //     };
    //     std::cout << std::endl;
    // };
    };
    return;
};
//------------------------------------------------------------------------------
//----------------------ELEMENT DIFFUSION/VISCOSITY MATRIX----------------------
//------------------------------------------------------------------------------
template<>
void Element<2,2>::getElemMatrixInitial(int &index, double *phi_, double **dphi_dx, double &tSUPG_, double &tPSPG_, double &tLSIC_, double &weight_, double **jacobianNRMatrix){
    
    tSUPG_ = 0.;
    // tPSPG_ = 0.;
    tLSIC_ = 0.;

    double &dTime_ = parameters.getTimeStep();
    double &visc_ = parameters.getViscosity();
    double &dens_ = parameters.getDensity();
    double &alpha_f = parameters.getAlphaF();
    double &alpha_m = parameters.getAlphaM();
    double &gamma = parameters.getGamma();

    double una_ = alpha_f * u_ + (1. - alpha_f) * uPrev_;
    double vna_ = alpha_f * v_ + (1. - alpha_f) * vPrev_;

    double duna_dx = alpha_f * du_dx + (1. - alpha_f) * duprev_dx;
    double duna_dy = alpha_f * du_dy + (1. - alpha_f) * duprev_dy;
    double dvna_dx = alpha_f * dv_dx + (1. - alpha_f) * dvprev_dx;
    double dvna_dy = alpha_f * dv_dy + (1. - alpha_f) * dvprev_dy;

    double umeshna_ = alpha_f * umesh_ + (1. - alpha_f) * umeshPrev_;
    double vmeshna_ = alpha_f * vmesh_ + (1. - alpha_f) * vmeshPrev_;
    
    for (int i = 0; i < 6; i++){
        for (int j = 0; j < 6; j++){

            double wSUPGi = (una_ - umeshna_) * dphi_dx[0][i] + (vna_ - vmeshna_) * dphi_dx[1][i];
            double wSUPGj = (una_ - umeshna_) * dphi_dx[0][j] + (vna_ - vmeshna_) * dphi_dx[1][j];
            
            //Mass matrix (use for both directions)
            double mM = phi_[i] * phi_[j] * dens_ * alpha_m + 
                        wSUPGi * phi_[j] * tSUPG_ * dens_ * alpha_m;

            jacobianNRMatrix[2*i  ][2*j  ] += mM * weight_ * djac_ * intPointWeightFunction[index];
            jacobianNRMatrix[2*i+1][2*j+1] += mM * weight_ * djac_ * intPointWeightFunction[index];
            
            //SINAL DA PARCELA QUE MULTIPLICA O TSUPG ESTA
            //COM SINAL TROCADO NA FORMULAÇAO DO TEZDUYAR
            //multipy pressure direction x
            double QSUPGx = - (dphi_dx[0][i] * phi_[j]) + 
                             ((dphi_dx[0][i] * (una_ - umeshna_) + 
                               dphi_dx[1][i] * (vna_ - vmeshna_)) * dphi_dx[0][j] * tSUPG_);
            //multiply pressure direction y
            double QSUPGy = - (dphi_dx[1][i] * phi_[j]) +
                             ((dphi_dx[0][i] * (una_ - umeshna_) + 
                               dphi_dx[1][i] * (vna_ - vmeshna_)) * dphi_dx[1][j] * tSUPG_);
            
            jacobianNRMatrix[2*i  ][12+j] += QSUPGx * weight_ * djac_ * intPointWeightFunction[index];
            jacobianNRMatrix[2*i+1][12+j] += QSUPGy * weight_ * djac_ * intPointWeightFunction[index];


            double Hx = dphi_dx[0][i] * phi_[j] * tPSPG_ * alpha_m*0;
            double Hy = dphi_dx[1][i] * phi_[j] * tPSPG_ * alpha_m*0;
            
            double Gx = dphi_dx[0][i] * wSUPGj * tPSPG_ * alpha_f * gamma * dTime_*0;
            double Gy = dphi_dx[1][i] * wSUPGj * tPSPG_ * alpha_f * gamma * dTime_*0;

            double Guu = (dphi_dx[0][i] * duna_dx * phi_[j] + 
                          dphi_dx[1][i] * dvna_dx * phi_[j]) * tPSPG_ * alpha_f * gamma * dTime_ * 0.;
            double Gvv = (dphi_dx[0][i] * duna_dy * phi_[j] +  
                          dphi_dx[1][i] * dvna_dy * phi_[j]) * tPSPG_ * alpha_f * gamma * dTime_ * 0.;
            
            double Q = (dphi_dx[0][i] * dphi_dx[0][j] + 
                        dphi_dx[1][i] * dphi_dx[1][j]) * tPSPG_ / (dens_);

            jacobianNRMatrix[12+j][2*i  ] += (Hx + Gx + Guu) * weight_ * djac_ * intPointWeightFunction[index];
            jacobianNRMatrix[12+j][2*i+1] += (Hy + Gy + Gvv) * weight_ * djac_ * intPointWeightFunction[index];
            jacobianNRMatrix[12+j][12+i] += Q * weight_ * djac_ * intPointWeightFunction[index];
        };
    };

    //std::cout << intPointWeightFunction(index) << std::endl;

    if (model){
        //        std::cout << " asas " << dens_ << " " << visc_ << " " << intPointWeightFunction(index) << std::endl;
    // std::cout << "matrix " << std::endl;
    // for (int i=0; i<18; i++){
    //     for (int j=0; j<18; j++){        
    //         std::cout << jacobianNRMatrix(i,j)<< " " ;
    //     };
    //     std::cout << std::endl;
    // };
    };
    return;
};


//------------------------------------------------------------------------------
//--------------------APPLY THE DIRICHLET BOUNDARY CONDITIONS-------------------
//------------------------------------------------------------------------------
template<>
void Element<2,2>::setBoundaryConditions(double **jacobianNRMatrix, double *rhsVector){

    for (int i = 0; i < 6; i++){


        if (((*nodes_)[connect_[i]] -> getConstrains(0) == 1) ||
            ((*nodes_)[connect_[i]] -> getConstrains(0) == 3)){//} ||
            for (int j = 0; j < 18; j++){
                jacobianNRMatrix[2*i  ][j] = 0.;
                jacobianNRMatrix[j][2*i  ] = 0.;
            };
            jacobianNRMatrix[2*i  ][2*i  ] = 1.;
            rhsVector[2*i  ] = 0.;
        };

        if (((*nodes_)[connect_[i]] -> getConstrains(1) == 1) ||
            ((*nodes_)[connect_[i]] -> getConstrains(1) == 3)){//} ||
            for (int j = 0; j < 18; j++){
                jacobianNRMatrix[2*i+1][j] = 0.;
                jacobianNRMatrix[j][2*i+1] = 0.;
            };
            jacobianNRMatrix[2*i+1][2*i+1] = 1.;
            rhsVector[2*i+1] =  0.;
        };

    
        // if ((model == false) && ((*nodes_)[connect_[i]] -> getDistFunction() > -1.2))  {

        //    // std::cout << "ELEMN " << index_ << " " << connect_[i] << std::endl;
        //     for (int j = 0; j < 18; j++){
        //         jacobianNRMatrix(2*i  ,j) *= 1000.e0;
        //         jacobianNRMatrix(j,2*i  ) *= 1000.e0;
        //         jacobianNRMatrix(2*i+1,j) *= 1000.e0;
        //         jacobianNRMatrix(j,2*i+1) *= 1000.e0;
        //     };
        //     rhsVector(2*i  ) *= 1000.e0;
        //     rhsVector(2*i+1) *= 1000.e0;
        // };

        

        // if ((model == false) && ((*nodes_)[connect_[i]] -> getDistFunction() > 0.3)){
        //     for (int j = 0; j < 18; j++){
        //         jacobianNRMatrix(12+i,j) = 0.;
        //         jacobianNRMatrix(j,12+i) = 0.;
        //     };
        //     jacobianNRMatrix(12+i,12+i) = 1.;
        //     rhsVector(12+i) =  0.;

        //     // for (int j = 0; j < 18; j++){
        //     //     jacobianNRMatrix(2*i  ,j) = 0.;
        //     //     jacobianNRMatrix(j,2*i  ) = 0.;
        //     //     jacobianNRMatrix(2*i+1,j) = 0.;
        //     //     jacobianNRMatrix(j,2*i+1) = 0.;                
        //     // };
        //     // jacobianNRMatrix(2*i  ,2*i  ) = 1.;
        //     // rhsVector(2*i  ) =  0.;
        //     // jacobianNRMatrix(2*i+1,2*i+1) = 1.;
        //     // rhsVector(2*i+1) =  0.;


        //     // std::cout<< "AUX mass " << index_ << " " << djac_ << std::endl;
        // }


        // double w = (*nodes_)[connect_[i]] -> getWeightFunction();
        // if (w < 0.012){
        //     rhsVector(2*i  ) = 0.;
        //     rhsVector(2*i+1) = 0.;
        // }


    };

    // typename Nodes::VecLocD x;

    // for (int i = 0; i < 6; i++){
    //     x = (*nodes_)[connect_[i]] -> getCoordinates();
    //     if((x(0) > 0.999) && (x(1) > 0.999)){
    //         for (int j = 0; j < 18; j++){
    //             jacobianNRMatrix(12+i,j) = 0.;
    //             jacobianNRMatrix(j,12+i) = 0.;
    //             jacobianNRMatrix(12+i,j) = 0.;
    //             jacobianNRMatrix(j,12+i) = 0.;
    //         };
    //         jacobianNRMatrix(12+i,12+i) = 1.;
    //         rhsVector(12+i) =  0.;
    //     };
    // };

    // for (int i = 0; i < 6; i++){
    //     x = (*nodes_)[connect_[i]] -> getCoordinates();
    //     if((x(0) > 0.999) || (x(0) < 0.001)){
    //         for (int j = 0; j < 18; j++){
    //             jacobianNRMatrix(2*i  ,j) = 0.;
    //             jacobianNRMatrix(j,2*i  ) = 0.;
    //             jacobianNRMatrix(2*i+1,j) = 0.;
    //             jacobianNRMatrix(j,2*i+1) = 0.;
    //         };
    //         jacobianNRMatrix(2*i+1,2*i+1) = 1.;
    //         jacobianNRMatrix(2*i  ,2*i  ) = 1.;
    //         rhsVector(2*i+1) =  0.;
    //         rhsVector(2*i  ) =  0.;
    //     };
    // };

    return;
};



//------------------------------------------------------------------------------
//--------------------APPLY THE DIRICHLET BOUNDARY CONDITIONS-------------------
//------------------------------------------------------------------------------
template<>
void Element<2,2>::setBoundaryConditions2(double **jacobianNRMatrix, double *rhsVector){

    for (int i = 0; i < 6; i++){


        if (((*nodes_)[connect_[i]] -> getConstrains(0) == 1) ||
            ((*nodes_)[connect_[i]] -> getConstrains(0) == 3)){//} ||
            for (int j = 0; j < 18; j++){
                jacobianNRMatrix[2*i  ][j] = 0.;
                jacobianNRMatrix[j][2*i  ] = 0.;
            };
            jacobianNRMatrix[2*i  ][2*i  ] = 1.;
            rhsVector[2*i  ] = 0.;
        };

        if (((*nodes_)[connect_[i]] -> getConstrains(1) == 1) ||
            ((*nodes_)[connect_[i]] -> getConstrains(1) == 3)){//} ||
            for (int j = 0; j < 18; j++){
                jacobianNRMatrix[2*i+1][j] = 0.;
                jacobianNRMatrix[j][2*i+1] = 0.;
            };
            jacobianNRMatrix[2*i+1][2*i+1] = 1.;
            rhsVector[2*i+1]=  0.;
        };
    };

    return;
};

//------------------------------------------------------------------------------
//--------------------APPLY THE DIRICHLET BOUNDARY CONDITIONS-------------------
//------------------------------------------------------------------------------
template<>
void Element<2,2>::setBoundaryConditionsLaplace(double **jacobianNRMatrix, double *rhsVector){

    for (int i = 0; i < 6; i++){

        if ((*nodes_)[connect_[i]] -> getConstrainsLaplace(0) == 1) {
            for (int j = 0; j < 18; j++){
                jacobianNRMatrix[2*i  ][j] = 0.;
                jacobianNRMatrix[j][2*i  ] = 0.;
            };
            jacobianNRMatrix[2*i  ][2*i  ] = 1.;
            rhsVector[2*i  ] = 0.;
        };

        if ((*nodes_)[connect_[i]] -> getConstrainsLaplace(1) == 1) {
            for (int j = 0; j < 18; j++){
                jacobianNRMatrix[2*i+1][j] = 0.;
                jacobianNRMatrix[j][2*i+1] = 0.;
            };
            jacobianNRMatrix[2*i+1][2*i+1] = 1.;
            rhsVector[2*i+1] =  0.;
        };
    };
    
    return;
};


//------------------------------------------------------------------------------
//---------------APPLY BOUNDARY CONDITIONS TO LAGRANGE MULTIPLIERS--------------
//------------------------------------------------------------------------------
template<>
void Element<2,2>::setBoundaryConditionsLagrangeMultipliers(double** jacobianNRMatrix, double* rhsVector){

    for (int i = 0; i < 18; i++) rhsVector[i] = 0.;
    
    double U_[18] = {};

    double &alpha_f = parameters.getAlphaF();

    for (int i=0; i<6; i++){
        // U_(2*i  ) = (*nodes_)[connect_[i]] -> getVelocity(0);
        // U_(2*i+1) = (*nodes_)[connect_[i]] -> getVelocity(1);
        U_[2*i  ] = alpha_f * (*nodes_)[connect_[i]] -> getVelocity(0) + (1. - alpha_f) * (*nodes_)[connect_[i]] -> getPreviousVelocity(0);
        U_[2*i+1] = alpha_f * (*nodes_)[connect_[i]] -> getVelocity(1) + (1. - alpha_f) * (*nodes_)[connect_[i]] -> getPreviousVelocity(1);
    };

    for (int i = 0; i < 18; i++)
        for (int j = 0; j < 18; j++)
            rhsVector[i] -= jacobianNRMatrix[i][j] * U_[j];
    
    for (int i = 0; i < 6; i++){
        if ((*nodes_)[connect_[i]] -> getConstrains(0) == 1) {
            for (int j = 0; j < 18; j++){
                //jacobianNRMatrix(2*i  ,j) = 0.;
                jacobianNRMatrix[j][2*i  ] = 0.;
            };
            //jacobianNRMatrix(12+i,12+i) = 1.;
            rhsVector[2*i  ] = 0.0;
        };

        if ((*nodes_)[connect_[i]] -> getConstrains(1) == 1) {
            // std::cout << "aqui" << std::endl;
            for (int j = 0; j < 18; j++){
                //jacobianNRMatrix(2*i+1,j) = 0.;
                jacobianNRMatrix[j][2*i+1] = 0.;
            };
            jacobianNRMatrix[12+i][12+i] = 1.;
            rhsVector[2*i+1] = 0.0;
        };
    };

    return;
};

//------------------------------------------------------------------------------
//-----------------------------RESIDUAL - RHS VECTOR----------------------------
//------------------------------------------------------------------------------
template<>
void Element<2,2>::getResidualVector(int &index, double *phi_, double **dphi_dx, double &tSUPG_, double &tPSPG_, double &tLSIC_, double &weight_, double *rhsVector){

    double &dTime_ = parameters.getTimeStep();
    double &visc_ = parameters.getViscosity();
    double &dens_ = parameters.getDensity();
    double &alpha_f = parameters.getAlphaF();
    double &alpha_m = parameters.getAlphaM();
    double &gamma = parameters.getGamma();
    double* &fieldForce = parameters.getFieldForce();

    double una_ = alpha_f * u_ + (1. - alpha_f) * uPrev_;
    double vna_ = alpha_f * v_ + (1. - alpha_f) * vPrev_;

    double axm_ = alpha_m * ax_ + (1. - alpha_m) * axprev_;
    double aym_ = alpha_m * ay_ + (1. - alpha_m) * ayprev_;

    double duna_dx = alpha_f * du_dx + (1. - alpha_f) * duprev_dx;
    double duna_dy = alpha_f * du_dy + (1. - alpha_f) * duprev_dy;
    double dvna_dx = alpha_f * dv_dx + (1. - alpha_f) * dvprev_dx;
    double dvna_dy = alpha_f * dv_dy + (1. - alpha_f) * dvprev_dy;

    double umeshna_ = alpha_f * umesh_ + (1. - alpha_f) * umeshPrev_;
    double vmeshna_ = alpha_f * vmesh_ + (1. - alpha_f) * vmeshPrev_;

    double wna_ = alpha_f * intPointWeightFunction[index] + (1. - alpha_f) * intPointWeightFunctionPrev[index];
    
    for (int i = 0; i < 6; i++){

        double mx = phi_[i] * axm_ * dens_ + 
                    ((una_ - umeshna_) * dphi_dx[0][i] + (vna_ - vmeshna_) * dphi_dx[1][i]) * 
                     axm_ * tSUPG_ * dens_;
        double my = phi_[i] * aym_ * dens_ + 
                    ((una_ - umeshna_) * dphi_dx[0][i] + (vna_ - vmeshna_) * dphi_dx[1][i]) * 
                     aym_ * tSUPG_ * dens_;
        
        double Kx = (2. * dphi_dx[0][i] * duna_dx * visc_ + 
                     dphi_dx[1][i] * duna_dy * visc_ + 
                     dphi_dx[1][i] * dvna_dx * visc_);
        double Ky = (dphi_dx[0][i] * duna_dy * visc_ + 
                     2. * dphi_dx[1][i] * dvna_dy * visc_ + 
                     dphi_dx[0][i] * dvna_dx * visc_);

        double KLSx = dphi_dx[0][i] * (duna_dx + dvna_dy) * tLSIC_ * dens_;
        double KLSy = dphi_dx[1][i] * (duna_dx + dvna_dy) * tLSIC_ * dens_;

        double Cx = (duna_dx * (una_ - umeshna_) + duna_dy * (vna_ - vmeshna_)) * phi_[i] * dens_ +
                    ((una_ - umeshna_) * dphi_dx[0][i] + (vna_ - vmeshna_) * dphi_dx[1][i]) *
                    ((una_ - umeshna_) * duna_dx + (vna_ - vmeshna_) * duna_dy) * tSUPG_ * dens_;
        double Cy = (dvna_dx * (una_ - umeshna_) + dvna_dy * (vna_ - vmeshna_)) * phi_[i] * dens_ +
                    ((una_ - umeshna_) * dphi_dx[0][i] + (vna_ - vmeshna_) * dphi_dx[1][i]) *
                    ((una_ - umeshna_) * dvna_dx + (vna_ - vmeshna_) * dvna_dy) * tSUPG_ * dens_;
 
        double Px = - (dphi_dx[0][i] * p_)
                    + ((dphi_dx[0][i] * (una_ - umeshna_) + dphi_dx[1][i] * (vna_ - vmeshna_)) * dp_dx * tSUPG_);
        double Py = - (dphi_dx[1][i] * p_) 
                    + ((dphi_dx[0][i] * (una_ - umeshna_) + dphi_dx[1][i] * (vna_ - vmeshna_)) * dp_dy * tSUPG_);
           
        double Q = ((duna_dx + dvna_dy) * phi_[i]) +
                   (dphi_dx[0][i] * dp_dx + dphi_dx[1][i] * dp_dy) * tPSPG_ / dens_ +
                    dphi_dx[0][i] * axm_ * tPSPG_ +
                    dphi_dx[1][i] * aym_ * tPSPG_ +
                    dphi_dx[0][i] * ((una_ - umeshna_) * duna_dx +
                                    (vna_ - vmeshna_) * duna_dy) * tPSPG_ + 
                    dphi_dx[1][i] * ((una_ - umeshna_) * dvna_dx +
                                    (vna_ - vmeshna_) * dvna_dy) * tPSPG_;

        if (model == false){

            // mx -= una_ * (intPointWeightFunction(index) - 
            //               intPointWeightFunctionPrev(index)) / dTime_ * phi_[i] * dens_;
            // my -= vna_ * (intPointWeightFunction(index) - 
            //               intPointWeightFunctionPrev(index)) / dTime_ * phi_[i] * dens_;
            // if (index_ == 2117) std::cout << "Int point " << intPointWeightFunction(index) << " " << intPointWeightFunctionPrev(index) << std::endl;
            // if (std::fabs((intPointWeightFunction(index) - intPointWeightFunctionPrev(index)) / dTime_) > 1.e-8) std::cout << "AQUI " << (intPointWeightFunction(index) - 
            //                     intPointWeightFunctionPrev(index)) / dTime_ << " " << una_ << " " << vna_ << std::endl;
        } else {
            // mx += una_ * dens_ * (dumesh_dx + dvmesh_dy) * phi_[i] * dens_;
            // my += vna_ * dens_ * (dumesh_dx + dvmesh_dy) * phi_[i] * dens_;
        }
        
        
        double Fx = fieldForce[0] * phi_[i];
        double Fy = fieldForce[1] * phi_[i];

        rhsVector[2*i  ] += (-mx -Kx - Px - Cx - KLSx + Fx) * weight_ * djac_ * wna_;
        rhsVector[2*i+1] += (-my -Ky - Py - Cy - KLSy + Fx) * weight_ * djac_ * wna_;
        rhsVector[12+i] += -Q * weight_ * djac_ * wna_;
                             
    };

    // if(model) {
    //     if (index_ == 601){ 
    //         std::cout << "RHS " << index << std::endl;
    //         for (int i=0; i<18; i++){
    //             std::cout << rhsVector(i) << " ";
    //         };
    //         std::cout << std::endl;
    //     };
    // };


    return;
};

//------------------------------------------------------------------------------
//-----------------------------RESIDUAL - RHS VECTOR----------------------------
//------------------------------------------------------------------------------
template<>
void Element<2,2>::getResidualVectorInitial(int &index, double *phi_, double **dphi_dx, double &tSUPG_, double &tPSPG_, double &tLSIC_, double &weight_, double *rhsVector){

    double &dTime_ = parameters.getTimeStep();
    double &visc_ = parameters.getViscosity();
    double &dens_ = parameters.getDensity();
    double &alpha_f = parameters.getAlphaF();
    double &alpha_m = parameters.getAlphaM();
    double &gamma = parameters.getGamma();
    double* &fieldForce = parameters.getFieldForce();

    double una_ = alpha_f * u_ + (1. - alpha_f) * uPrev_;
    double vna_ = alpha_f * v_ + (1. - alpha_f) * vPrev_;

    double unm_ = alpha_m * u_ + (1. - alpha_m) * uPrev_;
    double vnm_ = alpha_m * v_ + (1. - alpha_m) * vPrev_;

    double duna_dx = alpha_f * du_dx + (1. - alpha_f) * duprev_dx;
    double duna_dy = alpha_f * du_dy + (1. - alpha_f) * duprev_dy;
    double dvna_dx = alpha_f * dv_dx + (1. - alpha_f) * dvprev_dx;
    double dvna_dy = alpha_f * dv_dy + (1. - alpha_f) * dvprev_dy;
    
    double umeshna_ = alpha_f * umesh_ + (1. - alpha_f) * umeshPrev_;
    double vmeshna_ = alpha_f * vmesh_ + (1. - alpha_f) * vmeshPrev_;

    for (int i = 0; i < 6; i++){
        
        double Kx = (2. * dphi_dx[0][i] * duna_dx * visc_ + 
                     dphi_dx[1][i] * duna_dy * visc_ + 
                     dphi_dx[1][i] * dvna_dx * visc_);
        double Ky = (dphi_dx[0][i] * duna_dy * visc_ + 
                     2. * dphi_dx[1][i] * dvna_dy * visc_ + 
                     dphi_dx[0][i] * dvna_dx * visc_);

        double KLSx = dphi_dx[0][i] * (duna_dx + dvna_dy) * tLSIC_ * dens_;
        double KLSy = dphi_dx[1][i] * (duna_dx + dvna_dy) * tLSIC_ * dens_;

        double Cx = (duna_dx * (una_ - umeshna_) + duna_dy * (vna_ - vmeshna_)) * phi_[i] * dens_ +
                    ((una_ - umeshna_) * dphi_dx[0][i] + (vna_ - vmeshna_) * dphi_dx[1][i]) *
                    ((una_ - umeshna_) * duna_dx + (vna_ - vmeshna_) * duna_dy) * tSUPG_ * dens_;
        double Cy = (dvna_dx * (una_ - umeshna_) + dvna_dy * (vna_ - vmeshna_)) * phi_[i] * dens_ +
                    ((una_ - umeshna_) * dphi_dx[0][i] + (vna_ - vmeshna_) * dphi_dx[1][i]) *
                    ((una_ - umeshna_) * dvna_dx + (vna_ - vmeshna_) * dvna_dy) * tSUPG_ * dens_;
 

        double Q = ((duna_dx + dvna_dy) * phi_[i]);// +
                   // (dphi_dx[0][i] * dp_dx + dphi_dx[1][i] * dp_dy) * tPSPG_ / dens_ +
                   //  dphi_dx[0][i] * (unm_ - uPrev_) * tPSPG_ +
                   //  dphi_dx[1][i] * (vnm_ - vPrev_) * tPSPG_ +
                   //  dphi_dx[0][i] * ((una_ - umesh_) * duna_dx +
                   //                  (vna_ - vmesh_) * duna_dy) * tPSPG_ + 
                   //  dphi_dx[1][i] * ((una_ - umesh_) * dvna_dx +
                   //                  (vna_ - vmesh_) * dvna_dy) * tPSPG_;        
        
        double Fx = fieldForce[0] * phi_[i];
        double Fy = fieldForce[1] * phi_[i];

        rhsVector[2*i  ] += (-Kx - Cx - KLSx + Fx) * weight_ * djac_ * alpha_f * gamma * dTime_ * intPointWeightFunction[index];
        rhsVector[2*i+1] += (-Ky - Cy - KLSy + Fx) * weight_ * djac_ * alpha_f * gamma * dTime_ * intPointWeightFunction[index];
        rhsVector[12+i] += -Q * weight_ * djac_ * alpha_f * gamma * dTime_ * intPointWeightFunction[index];
                             
    };

    // if(model) {
    //     if (index_ == 601){ 
    //         std::cout << "RHS " << index << std::endl;
    //         for (int i=0; i<18; i++){
    //             std::cout << rhsVector(i) << " ";
    //         };
    //         std::cout << std::endl;
    //     };
    // };


    return;
};

//------------------------------------------------------------------------------
//-----------------------------RESIDUAL - RHS VECTOR----------------------------
//------------------------------------------------------------------------------

template<>
void Element<2,2>::getResidualVectorLaplace(double *rhsVector){

    // rhsVector.clear();
    
    double U_[18] = {};

    for (int i=0; i<6; i++){

        double *x_up = (*nodes_)[connect_[i]] -> getUpdatedCoordinates();        
        
        if ((*nodes_)[connect_[i]] -> getConstrainsLaplace(0) == 1){
            U_[2*i  ] = x_up[0] - (*nodes_)[connect_[i]] -> getCoordinateValue(0);
        };
        if ((*nodes_)[connect_[i]] -> getConstrainsLaplace(1) == 1){
            U_[2*i+1] = x_up[1] - (*nodes_)[connect_[i]] -> getCoordinateValue(1);
        };

        //if(connect_[i] == 29) std::cout << "Element 29 " << x_up(0) << " " << x_up(1) <<std::endl;
    };
    
    // noalias(rhsVector) += U_;//- prod(laplMatrix,U_); 

    return;
};

//------------------------------------------------------------------------------
//---------------------------ELEMENT LAPLACIAN MATRIX---------------------------
//------------------------------------------------------------------------------
template<>
void Element<2,2>::getElemLaplMatrix(double &weight_, double **dphi_dx, double **jacobianNRMatrix){

     
    for (int i = 0; i < 6; i++){
        for (int j = 0; j < 6; j++){        
            jacobianNRMatrix[2*i  ][2*j  ] += (dphi_dx[0][i] * dphi_dx[0][j] +
                                              dphi_dx[1][i] * dphi_dx[1][j]) 
                                        * weight_ * djac_ * meshMovingParameter;
            jacobianNRMatrix[2*i+1][2*j+1] += (dphi_dx[0][i] * dphi_dx[0][j] +
                                              dphi_dx[1][i] * dphi_dx[1][j]) 
                                        * weight_ * djac_ * meshMovingParameter;
        };
    };

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

    double xsi[2] = {};
    double phi_[6] = {};
    
    double **dphi_dx;
    dphi_dx = new double*[DIM];
    for (int i = DIM; i--; ) dphi_dx[i] = new double[nElNodes];

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
        xsi[0] = nQuad.PointList(index,0);
        xsi[1] = nQuad.PointList(index,1);

        //Computes the velocity shape functions
        shapeQuad.evaluate(xsi,phi_);

        //Returns the quadrature integration weight
        double weight_ = nQuad.WeightList(index);

        //Computes the jacobian matrix
        getJacobianMatrix(xsi, ainv_);

        //Computes spatial derivatives
        getSpatialDerivatives(xsi, ainv_, dphi_dx);

        //Interpolates velocity and its derivatives values
        getVelAndDerivatives(phi_, dphi_dx);

        //Compute Stabilization Parameters
        getParameterSUPG(tSUPG_, tPSPG_, tLSIC_, phi_, dphi_dx);

        //Computes the element diffusion/viscosity matrix
        getElemMatrix(index, phi_, dphi_dx, tSUPG_, tPSPG_, tLSIC_, weight_, jacobianNRMatrix);

        //Computes the RHS vector
        getResidualVector(index, phi_, dphi_dx, tSUPG_, tPSPG_, tLSIC_, weight_, rhsVector); 

        index++;        
    };  
    
    //Apply boundary conditions
    setBoundaryConditions(jacobianNRMatrix, rhsVector);

    for (int i = DIM; i--; ) delete [] dphi_dx[i];
    delete [] dphi_dx;
    for (int i = DIM; i--; ) delete [] ainv_[i];
    delete [] ainv_;

    return;
};

//------------------------------------------------------------------------------
//-----------------------TRANSIENT NAVIER-STOKES PROBEM-------------------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Element<DIM,DEG>::getTransientNavierStokesInitial(double *rhsVector, double **jacobianNRMatrix){

    double xsi[2] = {};
    double phi_[6] = {};  
    ShapeFunction           shapeQuad;

    double **dphi_dx;
    dphi_dx = new double*[DIM];
    for (int i = DIM; i--; ) dphi_dx[i] = new double[nElNodes];

    int index = 0;
    NormalQuad nQuad = NormalQuad();
    
    for (int i = 0; i < 18; i++){
        rhsVector[i] = 0.;
        for (int j = 0; j < 18; j++){
            jacobianNRMatrix[i][j] = 0.;
        }
    }
    

    double tSUPG_;
    double tPSPG_;
    double tLSIC_;

    double **ainv_;
    ainv_ = new double*[DIM];
    for (int i = DIM; i--; ) ainv_[i] = new double[DIM];

    for(double* it = nQuad.begin(); it != nQuad.end(); it++){

        //Defines the integration points adimentional coordinates
        xsi[0] = nQuad.PointList(index,0);
        xsi[1] = nQuad.PointList(index,1);

        //Computes the velocity shape functions
        shapeQuad.evaluate(xsi,phi_);

        //Returns the quadrature integration weight
        double weight_ = nQuad.WeightList(index);

        //Computes the jacobian matrix
        getJacobianMatrix(xsi, ainv_);

        //Computes spatial derivatives
        getSpatialDerivatives(xsi, ainv_, dphi_dx);

        //Interpolates velocity and its derivatives values
        getVelAndDerivatives(phi_, dphi_dx);

        //Compute Stabilization Parameters
        getParameterSUPG(tSUPG_, tPSPG_, tLSIC_, phi_, dphi_dx);

        //Computes the element diffusion/viscosity matrix
        getElemMatrixInitial(index, phi_, dphi_dx, tSUPG_, tPSPG_, tLSIC_, weight_, jacobianNRMatrix);

        //Computes the RHS vector
        getResidualVectorInitial(index, phi_, dphi_dx, tSUPG_, tPSPG_, tLSIC_, weight_, rhsVector); 

        index++;        
    };  
    
    //Apply boundary conditions
    setBoundaryConditions2(jacobianNRMatrix, rhsVector);

    for (int i = DIM; i--; ) delete [] dphi_dx[i];
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

    double xsi[2] = {};
    double phi_[6] = {}; 
    ShapeFunction           shapeQuad;
    
    double **dphi_dx;
    dphi_dx = new double*[DIM];
    for (int i = DIM; i--; ) dphi_dx[i] = new double[nElNodes];

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
        xsi[0] = nQuad.PointList(index,0);
        xsi[1] = nQuad.PointList(index,1);

        //Computes the velocity shape functions
        shapeQuad.evaluate(xsi,phi_);

        //Returns the quadrature integration weight
        double weight_ = nQuad.WeightList(index);

        //Computes the jacobian matrix
        getJacobianMatrix(xsi, ainv_);

        //Computes spatial derivatives
        getSpatialDerivatives(xsi, ainv_, dphi_dx);

        // ublas::bounded_matrix<double, 3,18> bMatrix,aux;
        // bMatrix.clear(); aux.clear();

        // for (int i = 0; i < 6; i++){                
        //     bMatrix(0,2*i  ) = dphi_dx[0][i];
        //     bMatrix(2,2*i  ) = dphi_dx[1][i];

        //     bMatrix(1,2*i+1) = dphi_dx[1][i];
        //     bMatrix(2,2*i+1) = dphi_dx[0][i];
        // };

        // aux = prod(hooke,bMatrix);        

        // laplMatrix += prod(trans(bMatrix),aux) * djac_ * weight_;

        getElemLaplMatrix(weight_, dphi_dx, jacobianNRMatrix);

        index++;        
    };  
    
    //Computes the RHS vector
    getResidualVectorLaplace(rhsVector);

    //Apply boundary conditions
    setBoundaryConditionsLaplace(jacobianNRMatrix, rhsVector);

    for (int i = DIM; i--; ) delete [] dphi_dx[i];
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

    double xsi[2] = {};
    double phi_[6] = {}; 
    double **dphi;
    dphi = new double*[DIM];
    for (int i = DIM; i--; ) dphi[i] = new double[nElNodes];
    
    ShapeFunction           shapeQuad;
    // ShapeFunctionDerivative dphi_dx;
    int index = 0;
    NormalQuad nQuad = NormalQuad();

    double &dTime_ = parameters.getTimeStep();
        
    for(double* it = nQuad.begin(); it != nQuad.end(); it++){
        
        //Defines the integration points adimentional coordinates
        xsi[0] = nQuad.PointList(index,0);
        xsi[1] = nQuad.PointList(index,1);

        //Computes the velocity shape functions
        shapeQuad.evaluate(xsi,phi_);
        shapeQuad.evaluateGradient(xsi,dphi);

        //Returns the quadrature integration weight
        double weight_ = nQuad.WeightList(index);

        //COMPUTE A0
        double dx_dxsi[2][2];
        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++)
                dx_dxsi[i][j] = 0.0;

        for (int i = 0; i < 6; ++i){
            double* initialCoord = (*nodes_)[connect_[i]] -> getInitialCoordinates();
            dx_dxsi[0][0] += initialCoord[0] * dphi[0][i];
            dx_dxsi[0][1] += initialCoord[0] * dphi[1][i];

            dx_dxsi[1][0] += initialCoord[1] * dphi[0][i];
            dx_dxsi[1][1] += initialCoord[1] * dphi[1][i];
        }
        double j0 = dx_dxsi[0][0] * dx_dxsi[1][1] - dx_dxsi[0][1] * dx_dxsi[1][0];
        
        //dxsi_dx
        double dxsi_dx[2][2];

        dxsi_dx[0][0] = dx_dxsi[1][1] / j0;
        dxsi_dx[0][1] = -dx_dxsi[0][1] / j0;
        dxsi_dx[1][0] = -dx_dxsi[1][0] / j0;
        dxsi_dx[1][1] = dx_dxsi[0][0] / j0;

        //dphi_dx
        double dphi_dx[6][2]; 
        for (int i = 0; i < 6; i++)
        {
            dphi_dx[i][0] = dxsi_dx[0][0] * dphi[0][i] + dxsi_dx[1][0] * dphi[1][i];
            dphi_dx[i][1] = dxsi_dx[0][1] * dphi[0][i] + dxsi_dx[1][1] * dphi[1][i];
        }


        //COMPUTE A1
        double dy_dxsi[2][2]; //row = cartesian, column = parametric
        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++)
                dy_dxsi[i][j] = 0.0;
        
        for (int i = 0; i < 6; i++){
            double* currentCoord = (*nodes_)[connect_[i]] -> getCoordinates();
            dy_dxsi[0][0] += currentCoord[0] * dphi[0][i];
            dy_dxsi[0][1] += currentCoord[0] * dphi[1][i];

            dy_dxsi[1][0] += currentCoord[1] * dphi[0][i];
            dy_dxsi[1][1] += currentCoord[1] * dphi[1][i];
        }

        //dy_dx
        double dy_dx[2][2];
        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++)
                dy_dx[i][j] = dy_dxsi[i][0] * dxsi_dx[0][j] + dy_dxsi[i][1] * dxsi_dx[1][j];

        //jacobian
        double jac = dy_dx[0][0] * dy_dx[1][1] - dy_dx[0][1] * dy_dx[1][0];

        //Green-Lagrange strain tensor
        double E[2][2];
        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++)
                E[i][j] = 0.5 * (dy_dx[0][i] * dy_dx[0][j] + dy_dx[1][i] * dy_dx[1][j]);
        E[0][0] -= 0.5; E[1][1] -= 0.5;

        //Second Piola-Kirchhoff stress tensor
        double S[2][2];
        double young = 1. / j0;
        double poisson = 0.3;
        S[0][0] = young / (1.0-(poisson*poisson)) * (E[0][0] + poisson * E[1][1]);
        S[0][1] = 2.0 * (young / (2.0 * (1.0+poisson))) * E[0][1];
        S[1][0] = S[0][1];
        S[1][1] = young / (1.0-(poisson*poisson)) * (E[1][1] + poisson * E[0][0]);
        //     S[0][0] = (young / ( (1.0+poisson) * (1.0-2.0*poisson))) * ((1.0-poisson) * E[0][0] + poisson * E[1][1]);
        //     S[1][1] = (young / ( (1.0+poisson) * (1.0-2.0*poisson))) * ((1.0-poisson) * E[1][1] + poisson * E[0][0]);


        //element rhs vector
        for (int a = 0; a < 6; a++){
            for (int k = 0; k < 2; k++){
                double dE_dyak[2][2];
                for (int i = 0; i < 2; i++)
                    for (int j = 0; j < 2; j++)
                        dE_dyak[i][j] = 0.5 * (dphi_dx[a][i] * dy_dx[k][j] + dy_dx[k][i] * dphi_dx[a][j]);
                
                //internal force
                double f = 0.0;
                for (int i = 0; i < 2; i++)
                    for (int j = 0; j < 2; j++)
                        f += S[i][j] * dE_dyak[i][j];

                // double accel = 0.0;
                // for (int i = 0; i < 6; i++)
                //     accel += phi(i) * nodes_[i]->getCurrentAcceleration()(k);
                
                // //inertial force
                // double m = density * phi(a) * accel;

                double vel = 0.0;
                for (int i = 0; i < 6; i++)
                    vel += phi_[i] * (*nodes_)[connect_[i]]->getMeshVelocity(k);

                double c =  phi_[a] * vel*0;

                //domain force
                // double b;
                // (k==1) ? b = phi(a) * density * gravity : 0.0;

                rhsVector[2 * a + k] -= (f+c) * weight_ * j0;

                //element tangent matrix
                for (int b = 0; b < 6; b++){
                    for (int l = 0; l < 2; l++){
                        double dE_dybl[2][2];
                        for (int i = 0; i < 2; i++)
                            for (int j = 0; j < 2; j++)
                                dE_dybl[i][j] = 0.5 * (dphi_dx[b][i] * dy_dx[l][j] + dy_dx[l][i] * dphi_dx[b][j]);

                        double d2E_dyakbl[2][2];
                        if (k==l)
                            for (int i = 0; i < 2; i++)
                                for (int j = 0; j < 2; j++)
                                    d2E_dyakbl[i][j] = 0.5 * (dphi_dx[a][i] * dphi_dx[b][j] + dphi_dx[b][i] * dphi_dx[a][j]);
                        else
                            for (int i = 0; i < 2; i++)
                                for (int j = 0; j < 2; j++)
                                    d2E_dyakbl[i][j] = 0.0; 

                        double dS_dybl[2][2];
                        dS_dybl[0][0] = young / (1.0-(poisson*poisson)) * (dE_dybl[0][0] + poisson * dE_dybl[1][1]);
                        dS_dybl[0][1] = 2.0 * (young / (2.0 * (1.0+poisson))) * dE_dybl[0][1];
                        dS_dybl[1][0] = 2.0 * (young / (2.0 * (1.0+poisson))) * dE_dybl[1][0];
                        dS_dybl[1][1] = young / (1.0-(poisson*poisson)) * (dE_dybl[1][1] + poisson * dE_dybl[0][0]);
                        // if (analysisType == "PLANE_STRAIN")
                        // {
                        //     dS_dybl[0][0] = (young / ( (1.0+poisson) * (1.0-2.0*poisson))) * ((1.0-poisson) * dE_dybl[0][0] + poisson * dE_dybl[1][1]);
                        //     dS_dybl[1][1] = (young / ( (1.0+poisson) * (1.0-2.0*poisson))) * ((1.0-poisson) * dE_dybl[1][1] + poisson * dE_dybl[0][0]);
                        // }

                        //elastic and geometric componentes of tangent matrix
                        double e = 0.0;
                        for (int i = 0; i < 2; i++)
                            for (int j = 0; j < 2; j++)
                                e += dS_dybl[i][j] * dE_dyak[i][j] + S[i][j] * d2E_dyakbl[i][j];

                        //mass matrix
                        double m;
                        // (k==l)? m = (1.0 / (beta * deltat * deltat)) * density * phi(a) * phi(b) : m = 0.0;
                        (k==l)? m = (1.0 / (0.25 * dTime_)) *0* phi_[a] * phi_[b] : m = 0.0;

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
    

    // if (index_ == 0) std::cout << "asdas  " << jacobianNRMatrix << " " << rhsVector << std::endl;

    //Computes the RHS vector
    // getResidualVectorLaplace(rhsVector);

    //Apply boundary conditions
    setBoundaryConditionsLaplace(jacobianNRMatrix, rhsVector);

    for (int i = DIM; i--; ) delete [] dphi[i];
    delete [] dphi;

    return;
};

//------------------------------------------------------------------------------
//----------------------------STEADY LAPLACE PROBEM-----------------------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Element<DIM,DEG>::getLagrangeMultipliersSameMesh(double **lagrMultMatrix, double *lagrMultVector, double *rhsVector){

    double xsi[2] = {};
    double phi_[6] = {};
    
    double **dphi_dx;
    dphi_dx = new double*[DIM];
    for (int i = DIM; i--; ) dphi_dx[i] = new double[nElNodes];

    double **ainv_;
    ainv_ = new double*[DIM];
    for (int i = DIM; i--; ) ainv_[i] = new double[DIM];

    QuadShapeFunction<2,2> shapeQuad;
    int index = 0;
    NormalQuad nQuad = NormalQuad();

    double diffMatrix[18][18] = {};

    double tSUPG_; double tPSPG_; double tLSIC_;

    double &k1 = parameters.getArlequinK1();
    double &k2 = parameters.getArlequinK2();
    double &alpha_f = parameters.getAlphaF();
    double &gamma = parameters.getGamma();
    double &dTime_ = parameters.getTimeStep();
    
    for(double* it = nQuad.begin(); it != nQuad.end(); it++){
        
        //Defines the integration points adimentional coordinates
        xsi[0] = nQuad.PointList(index,0);
        xsi[1] = nQuad.PointList(index,1);

        //Computes the velocity shape functions
        shapeQuad.evaluate(xsi,phi_);

        //Returns the quadrature integration weight
        double weight_ = nQuad.WeightList(index);

        //Computes the jacobian matrix
        getJacobianMatrix(xsi, ainv_);

        getSpatialDerivatives(xsi, ainv_, dphi_dx);
        
        getVelAndDerivatives(phi_, dphi_dx);

        // H1 COUPLING OPERATOR  
        if (k2 > 0){      
            double Bmatrix[3][18] = {};
            double Baux[3][18] = {};
            double BmatrixT[18][3] = {}; 
            double Maux[3][3] = {};

            for (int i = 0; i < 6; i++){
                Bmatrix[0][2*i  ] = dphi_dx[0][i];
                Bmatrix[2][2*i  ] = dphi_dx[1][i];       
                
                BmatrixT[2*i  ][0] = dphi_dx[0][i];
                BmatrixT[2*i  ][2] = dphi_dx[1][i];        
                
                Bmatrix[1][2*i+1] = dphi_dx[1][i];
                Bmatrix[2][2*i+1] = dphi_dx[0][i];
                
                BmatrixT[2*i+1][1] = dphi_dx[1][i];
                BmatrixT[2*i+1][2] = dphi_dx[0][i];
            };
            
            Maux[0][0] = 2.; 
            Maux[1][1] = 2.;
            Maux[2][2] = 1.;
            
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 18; j++)
                    for (int k = 0; k < 3; k++)
                        Baux[i][j] += Maux[i][k] * Bmatrix[k][j];

            for (int i = 0; i < 18; i++)
                for (int j = 0; j < 18; j++)
                    for (int k = 0; k < 3; k++)
                        diffMatrix[i][j] += BmatrixT[i][k] * Baux[k][j] * weight_ * djac_ * k2;
        }

        double una_ = alpha_f * u_ + (1. - alpha_f) * uPrev_;
        double vna_ = alpha_f * v_ + (1. - alpha_f) * vPrev_;
        // double una_ = uPrev_ + alpha_f * dTime_ * ((1. - gamma) * axprev_ + gamma * ax_);
        // double vna_ = vPrev_ + alpha_f * dTime_ * ((1. - gamma) * ayprev_ + gamma * ay_);

        double duna_dx = alpha_f * du_dx + (1. - alpha_f) * duprev_dx;
        double duna_dy = alpha_f * du_dy + (1. - alpha_f) * duprev_dy;
        double dvna_dx = alpha_f * dv_dx + (1. - alpha_f) * dvprev_dx;
        double dvna_dy = alpha_f * dv_dy + (1. - alpha_f) * dvprev_dy;
        
        for (int i = 0; i < 6; i++){
            for (int j = 0; j < 6; j++){
                // L2 COUPLING OPERATOR
                lagrMultMatrix[2*i  ][2*j  ] -= phi_[i] * phi_[j] * weight_ * djac_ * k1;
                lagrMultMatrix[2*i+1][2*j+1] -= phi_[i] * phi_[j] * weight_ * djac_ * k1;
            };
            // Lagrange multipliers residual
            lagrMultVector[2*i  ] += (lagMx_ * phi_[i] * k1 +
                                      (2. * dphi_dx[0][i] * dLx_dx + 
                                       dphi_dx[1][i] * dLx_dy + 
                                       dphi_dx[1][i] * dLy_dx) * k2) * weight_ * djac_;
            lagrMultVector[2*i+1] += (lagMy_ * phi_[i] * k1 +
                                      (dphi_dx[0][i] * dLx_dy + 
                                       2. * dphi_dx[1][i] * dLy_dy + 
                                       dphi_dx[0][i] * dLy_dx) * k2) * weight_ * djac_;

            rhsVector[2*i  ] += (una_ * phi_[i] * k1 +
                                (2. * dphi_dx[0][i] * duna_dx + 
                                 dphi_dx[1][i] * duna_dy + 
                                 dphi_dx[1][i] * dvna_dx) * k2) * weight_ * djac_;
            rhsVector[2*i+1] += (vna_ * phi_[i] * k1 +
                                (dphi_dx[0][i] * duna_dy + 
                                 2. * dphi_dx[1][i] * dvna_dy + 
                                 dphi_dx[0][i] * dvna_dx) * k2) * weight_ * djac_;
        };         
        
        index++; 
    }; 

    if (k2 > 0)
        for (int i = 0; i < 18; i++)
            for (int j = 0; j < 18; j++)
                lagrMultMatrix[i][j] -= diffMatrix[i][j];


    for (int i = DIM; i--; ) delete [] dphi_dx[i];
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

    double xsi[2] = {};
    double phi_[6] = {};
    
    double **dphi_dx;
    dphi_dx = new double*[DIM];
    for (int i = DIM; i--; ) dphi_dx[i] = new double[nElNodes];

    double **ainv_;
    ainv_ = new double*[DIM];
    for (int i = DIM; i--; ) ainv_[i] = new double[DIM];  

    QuadShapeFunction<2,2> shapeQuad;
    int index = 0;
    NormalQuad nQuad = NormalQuad();

    double tSUPG_; double tPSPG_; double tLSIC_;

    double &dens_ = parameters.getDensity();
    double &alpha_f = parameters.getAlphaF();
    double &k1 = parameters.getArlequinK1();
    
    for(double* it = nQuad.begin(); it != nQuad.end(); it++){
        
        //Defines the integration points adimentional coordinates
        xsi[0] = nQuad.PointList(index,0);
        xsi[1] = nQuad.PointList(index,1);

        //Computes the velocity shape functions
        shapeQuad.evaluate(xsi,phi_);

        //Returns the quadrature integration weight
        double weight_ = nQuad.WeightList(index);

        //Computes the jacobian matrix
        getJacobianMatrix(xsi, ainv_);

        getSpatialDerivatives(xsi, ainv_, dphi_dx);
        
        getVelAndDerivatives(phi_, dphi_dx);

        getParameterSUPG(tSUPG_, tPSPG_, tLSIC_, phi_, dphi_dx);

        double una_ = alpha_f * u_ + (1. - alpha_f) * uPrev_;
        double vna_ = alpha_f * v_ + (1. - alpha_f) * vPrev_;

        double umeshna_ = alpha_f * umesh_ + (1. - alpha_f) * umeshPrev_;
        double vmeshna_ = alpha_f * vmesh_ + (1. - alpha_f) * vmeshPrev_;

        double wna_ = alpha_f * intPointWeightFunction[index] + (1. - alpha_f) * intPointWeightFunctionPrev[index];
        
        for (int i = 0; i < 6; i++){
            for (int j = 0; j < 6; j++){        

                //SUPG STABILIZATION TERM
                double LM = 0.;
                // LM = - ((una_ - umeshna_) * dphi_dx[0][i] + (vna_ - vmeshna_) * dphi_dx[1][i])
                //     * phi_[j] * tSUPG_;
            
                jacobianNRMatrix[2*i  ][2*j  ] += LM * weight_ * djac_ * k1;
                jacobianNRMatrix[2*i+1][2*j+1] += LM * weight_ * djac_ * k1;

                //PSPG STABILIZATION TERM
                double Lx = 0.;
                double Ly = 0.;
                double Tx = 0.;
                double Ty = 0.;

                // Lx = dphi_dx[0][i] * phi_[j] * tPSPG_ / dens_ * k1;
                // Ly = dphi_dx[1][i] * phi_[j] * tPSPG_ / dens_ * k1;

                jacobianNRMatrix[2*i  ][12+j] += Lx * weight_ * djac_;
                jacobianNRMatrix[2*i+1][12+j] += Ly * weight_ * djac_;
            };

            //SUPG AND PSPG STABILIZATION TERMS
            double LMx = 0.;
            double LMy = 0.;
            double LMp = 0.;
            
            // LMx = (((una_ - umeshna_) * dphi_dx[0][i] + (vna_ - vmeshna_) * dphi_dx[1][i])
            //        * lagMx_) * tSUPG_;
            // LMy = (((una_ - umeshna_) * dphi_dx[0][i] + (vna_ - vmeshna_) * dphi_dx[1][i])
            //        * lagMy_) * tSUPG_;
            // LMp = (dphi_dx[0][i] * lagMx_ + dphi_dx[1][i] * lagMy_) * tPSPG_/ dens_;

            rhsVector[2*i  ] += (-LMx) * weight_ * djac_ * k1;
            rhsVector[2*i+1] += (-LMy) * weight_ * djac_ * k1;
            rhsVector[12+i] += (-LMp) * weight_ * djac_ * k1;

        };         
        
        index++;        
    }; 

    for (int i = DIM; i--; ) delete [] dphi_dx[i];
    delete [] dphi_dx;
    for (int i = DIM; i--; ) delete [] ainv_[i];
    delete [] ainv_;

    return;
};

//------------------------------------------------------------------------------
//----------------------------STEADY LAPLACE PROBEM-----------------------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Element<DIM,DEG>::getLagrangeMultipliersArlequinSameMesh(double **arlequinStab, double **laplMatrix, double *arlequinStabVector){

    double xsi[2] = {};
    double phi_[6] = {};
    
    double **dphi_dx;
    dphi_dx = new double*[DIM];
    for (int i = DIM; i--; ) dphi_dx[i] = new double[nElNodes];

    double **ainv_;
    ainv_ = new double*[DIM];
    for (int i = DIM; i--; ) ainv_[i] = new double[DIM];  

    QuadShapeFunction<2,2> shapeQuad;
    int index = 0;
    NormalQuad nQuad = NormalQuad();

    double tSUPG_; double tPSPG_; double tLSIC_;

    double &dens_ = parameters.getDensity();
    double &alpha_f = parameters.getAlphaF();
    double &k1 = parameters.getArlequinK1();

    getParameterArlequin2();

    for(double* it = nQuad.begin(); it != nQuad.end(); it++){
        //Defines the integration points adimentional coordinates
        xsi[0] = nQuad.PointList(index,0);
        xsi[1] = nQuad.PointList(index,1);

        //Computes the velocity shape functions
        shapeQuad.evaluate(xsi,phi_);

        //Returns the quadrature integration weight
        double weight_ = nQuad.WeightList(index);

        //Computes the jacobian matrix
        getJacobianMatrix(xsi, ainv_);

        getSpatialDerivatives(xsi, ainv_, dphi_dx);
        
        getVelAndDerivatives(phi_, dphi_dx);

        getParameterArlequin(tSUPG_, tPSPG_, tLSIC_, phi_, dphi_dx);

        double wna_ = alpha_f * intPointWeightFunction[index] + (1. - alpha_f) * intPointWeightFunctionPrev[index];
        
        for (int i = 0; i < 6; i++){
            for (int j = 0; j < 6; j++){        

                //ARLEQUIN STABILIZATION TERMS
                double AM = 0.;
                double Lpx = 0.; double Lpy = 0.;
                double LC = 0.; double LL = 0.;

                // AM = phi_[i] * phi_[j] * tARLQ_ * intPointWeightFunction(index) * alpha_m;

                // LL = -2 * phi_[i] * phi_[j] * tARLQ_ / dens_;
                LL = (dphi_dx[0][i] * dphi_dx[0][j] + dphi_dx[1][i] * dphi_dx[1][j]) * tARLQ_ / dens_;

                // Lpx = -(dphi_dx[0][i] * dp_dxx + dphi_dx[1][i] * dp_dxy) * intPointWeightFunction(index)
                //      * tARLQ_ / dens_;
                // Lpy = -(dphi_dx[0][i] * dp_dxy + dphi_dx[1][i] * dp_dyy) * intPointWeightFunction(index)
                //      * tARLQ_ / dens_;

                // LC = phi_[i] * ((una_ - umesh_) * dphi_dx[0][i] + (vna_ - vmesh_) * dphi_dx[1][i]) * phi_[j] * tARLQ_ * intPointWeightFunction(index);

                arlequinStab[2*i  ][2*j  ] += LL * weight_ * djac_;
                arlequinStab[2*i+1][2*j+1] += LL * weight_ * djac_;
                arlequinStab[12+i][12+j] += 0 * weight_ * djac_;

                // LC = -(dphi_dx[0][i]*(du_dx*dphi_dx[0][j] + dv_dx*dphi_dx[1][j]) +
                //        dphi_dx[1][i]*(du_dy*dphi_dx[0][j] + dv_dy*dphi_dx[1][j]) + 
                //        dphi_dx[0][i]*(u_*ddphi_dx(0,0)(j) + v_*ddphi_dx(0,1)(j)) + 
                //        dphi_dx[1][i]*(u_*ddphi_dx(1,0)(j) + v_*ddphi_dx(1,1)(j))) * tARLQ_ * intPointWeightFunction(index);

                // Lpx = 0.; Lpy = 0.;

                // Lpx = (dphi_dx[0][i] * ddphi_dx(0,0)(j) + 
                //        dphi_dx[1][i] * ddphi_dx(0,1)(j)) * tARLQ_ * intPointWeightFunction(index);
                // Lpy = (dphi_dx[0][i] * ddphi_dx(1,0)(j) + 
                //        dphi_dx[1][i] * ddphi_dx(1,1)(j)) * tARLQ_ * intPointWeightFunction(index);




                laplMatrix[2*i  ][2*j  ] += (LC + AM) * weight_ * djac_;
                laplMatrix[2*i+1][2*j+1] += (LC + AM) * weight_ * djac_;

                laplMatrix[2*i  ][12+j] += Lpx * weight_ * djac_;
                laplMatrix[2*i+1][12+j] += Lpy * weight_ * djac_;
                laplMatrix[12+j][2*i  ] += Lpx * weight_ * djac_;
                laplMatrix[12+j][2*i+1] += Lpy * weight_ * djac_;

            };

            //ARLEQUIN STABILIZATION TERMS
            double Amx = 0.; double Amy = 0.;
            double LCx = 0.; double LCy = 0.;
            double LPx = 0.; double LPy = 0.;
            double LLx = 0.; double LLy = 0.;

            // if (iTimeStep > 5){
                // Amx = -phi_[i] * axm_ * tARLQ_;
                // Amy = -phi_[i] * aym_ * tARLQ_;

            // }

// Iteration = 5 (1)  Du Norm = 3.465085e+01 2.740808e+00 1.351008e+00 3.367368e+01
            //if (iTimeStep < 10){
                LLx = -(dphi_dx[0][i] * (dLx_dx/wna_) + dphi_dx[1][i] * (dLx_dy/wna_)) * tARLQ_ / dens_;
                LLy = -(dphi_dx[0][i] * (dLy_dx/wna_) + dphi_dx[1][i] * (dLy_dy/wna_)) * tARLQ_ / dens_;
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
            
            arlequinStabVector[2*i  ] += (Amx + LLx) * weight_ * djac_ * wna_;
            arlequinStabVector[2*i+1] += (Amy + LLy) * weight_ * djac_ * wna_;
        };         
        
        index++;        
    };  

    for (int i = 0; i < 6; i++){
        if ((*nodes_)[connect_[i]] -> getConstrains(0) == 1) {
            for (int j = 0; j < 18; j++){
                arlequinStab[2*i  ][j] = 0.;
                arlequinStab[j][2*i  ] = 0.;
            };
            //lagrMultMatrix(12+i,12+i) = 1.;
            arlequinStabVector[2*i  ] = 0.0;
        };

        if ((*nodes_)[connect_[i]] -> getConstrains(1) == 1) {
            // std::cout << "aqui" << std::endl;
            for (int j = 0; j < 18; j++){
                arlequinStab[2*i+1][j] = 0.;
                arlequinStab[j][2*i+1] = 0.;
            };
            //lagrMultMatrix(12+i,12+i) = 1.;
            arlequinStabVector[2*i+1] = 0.0;
        };
    };

    for (int i = DIM; i--; ) delete [] dphi_dx[i];
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

    double xsi[2]= {};
    double xsi_intp[2] = {};
    QuadShapeFunction<2,2> shapeQuad;
    int index = 0;
    SpecialQuad sQuad = SpecialQuad();

    //tARLQ_ = -tPSPG2_;

    double phiLM_[6] = {};
    double phi_[6] = {};
    
    double **dphi_dx;
    dphi_dx = new double*[DIM];
    for (int i = DIM; i--; ) dphi_dx[i] = new double[nElNodes];
    double **dphiL_dx;
    dphiL_dx = new double*[DIM];
    for (int i = DIM; i--; ) dphiL_dx[i] = new double[nElNodes];        
    

    double &dTime_ = parameters.getTimeStep();
    double &visc_ = parameters.getViscosity();
    double &dens_ = parameters.getDensity();
    double &alpha_f = parameters.getAlphaF();
    double &alpha_m = parameters.getAlphaM();
    double &gamma = parameters.getGamma();
    double &k1 = parameters.getArlequinK1();
    double &k2 = parameters.getArlequinK2();
    int &iTimeStep = parameters.getTimeInstant();

    double diffMatrix[18][18] = {};

    double **ainv_;
    ainv_ = new double*[DIM];
    for (int i = DIM; i--; ) ainv_[i] = new double[DIM];

    // if (index_ == 4936) std::cout << "PSPG Fine " << ielem << " " << tPSPG_ << std::endl;
    for(double* it = sQuad.begin(); it != sQuad.end(); it++){
        
        if ((intPointCorrespElem[index] == ielem)){

            //Defines the integration points adimentional coordinates
            xsi[0] = sQuad.PointList(index,0);
            xsi[1] = sQuad.PointList(index,1);
            
            //Computes the velocity shape functions
            shapeQuad.evaluate(xsi,phi_);
            
            xsi_intp[0] = intPointCorrespXsi[index][0];
            xsi_intp[1] = intPointCorrespXsi[index][1];

            //Computes the coarse mesh shape functions
            shapeQuad.evaluate(xsi_intp,phiLM_);
            
            //Returns the quadrature integration weight
            double weight_ = sQuad.WeightList(index);
            
            //Computes the jacobian matrix
            getJacobianMatrix(xsi_intp, ainv_);
                        
            getSpatialDerivatives(xsi_intp, ainv_, dphi_dx);

            for (int i = 0; i < 2; ++i)
                for (int j = 0; j < 6; ++j)
                    dphiL_dx[i][j] = dphi_dx[i][j];
            // ddphiL_dx = ddphi_dx;

            getJacobianMatrix(xsi, ainv_);
            getSpatialDerivatives(xsi, ainv_, dphi_dx);
            getVelAndDerivatives(phi_, dphi_dx);

            double umeshna_ = alpha_f * umesh_ + (1. - alpha_f) * umeshPrev_;
            double vmeshna_ = alpha_f * vmesh_ + (1. - alpha_f) * vmeshPrev_;

            double dumeshna_dx = alpha_f * dumesh_dx + (1. - alpha_f) * dumeshPrev_dx;
            double dumeshna_dy = alpha_f * dumesh_dy + (1. - alpha_f) * dumeshPrev_dy;
            double dvmeshna_dx = alpha_f * dvmesh_dx + (1. - alpha_f) * dvmeshPrev_dx;
            double dvmeshna_dy = alpha_f * dvmesh_dy + (1. - alpha_f) * dvmeshPrev_dy;

            double una_ = alpha_f * u_ + (1. - alpha_f) * uPrev_;
            double vna_ = alpha_f * v_ + (1. - alpha_f) * vPrev_;

            double duna_dx = alpha_f * du_dx + (1. - alpha_f) * duprev_dx;
            double duna_dy = alpha_f * du_dy + (1. - alpha_f) * duprev_dy;
            double dvna_dx = alpha_f * dv_dx + (1. - alpha_f) * dvprev_dx;
            double dvna_dy = alpha_f * dv_dy + (1. - alpha_f) * dvprev_dy;

            double wna_ = alpha_f * intPointWeightFunctionSpecial[index] + (1. - alpha_f) * intPointWeightFunctionSpecialPrev[index];

            double unaL_ = 0.;
            double vnaL_ = 0.;
            double dalpha_dx = 0.;
            double dalpha_dy = 0.;
            for (int i = 0; i < 6; ++i){
                // dalpha_dx += (alpha_f * (*nodes_)[connect_[i]] -> getWeightFunction() + (1. - alpha_f) * (*nodes_)[connect_[i]] -> getPreviousWeightFunction()) * dphi_dx[0][i];
                // dalpha_dy += (alpha_f * (*nodes_)[connect_[i]] -> getWeightFunction() + (1. - alpha_f) * (*nodes_)[connect_[i]] -> getPreviousWeightFunction()) * dphi_dx[1][i];
                dalpha_dx += ((*nodes_)[connect_[i]] -> getWeightFunction()) * dphi_dx[0][i];
                dalpha_dy += ((*nodes_)[connect_[i]] -> getWeightFunction()) * dphi_dx[1][i];
                unaL_ += alpha_f * velx[i] * phiLM_[i] + (1. - alpha_f) * velxPrev[i] * phiLM_[i];
                vnaL_ += alpha_f * vely[i] * phiLM_[i] + (1. - alpha_f) * velyPrev[i] * phiLM_[i];
            }

            // if (index_ == 3056 ) std::cout << "DADX " << dalpha_dx << " " << dalpha_dy << " " << index << " " << ielem << std::endl;

            if (k2 > 0){
                //H1 COUPLING OPERATOR
                double Bmatrix[3][18] = {};
                double Baux[3][18] = {};
                double BmatrixT[18][3] = {}; 
                double Maux[3][3] = {};        
                
                for (int i = 0; i < 6; i++){
                    Bmatrix[0][2*i  ] = dphi_dx[0][i];
                    Bmatrix[2][2*i  ] = dphi_dx[1][i];       
                    
                    BmatrixT[2*i  ][0] = dphiL_dx[0][i];
                    BmatrixT[2*i  ][2] = dphiL_dx[1][i];        
                    
                    Bmatrix[1][2*i+1] = dphi_dx[1][i];
                    Bmatrix[2][2*i+1] = dphi_dx[0][i];
                    
                    BmatrixT[2*i+1][1] = dphiL_dx[1][i];
                    BmatrixT[2*i+1][2] = dphiL_dx[0][i];
                };

                Maux[0][0] = 2.; 
                Maux[1][1] = 2.;
                Maux[2][2] = 1.;
                
                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 18; j++)
                        for (int k = 0; k < 3; k++)
                            Baux[i][j] += Maux[i][k] * Bmatrix[k][j];

                for (int i = 0; i < 18; i++)
                    for (int j = 0; j < 18; j++)
                        for (int k = 0; k < 3; k++)
                            diffMatrix[i][j] += BmatrixT[i][k] * Baux[k][j] * weight_ * djac_ * k2;

            }

            for (int i = 0; i < 6; i++){
                for (int j = 0; j < 6; j++){     
                    //L2 COUPLING OPERATOR   
                    lagrMultMatrix[2*j  ][2*i  ] += (phiLM_[i] * phi_[j]) * weight_ * djac_ * k1;
                    lagrMultMatrix[2*j+1][2*i+1] += (phiLM_[i] * phi_[j]) * weight_ * djac_ * k1;
                };
                
                rhsVectorLM[2*i  ] += -(lagMx_ * phiLM_[i] * k1 + 
                                        (2. * dphiL_dx[0][i] * dLx_dx + 
                                         dphiL_dx[1][i] * dLx_dy + 
                                         dphiL_dx[1][i] * dLy_dx) * k2) * weight_ * djac_;
                rhsVectorLM[2*i+1] += -(lagMy_ * phiLM_[i] * k1 + 
                                        (dphiL_dx[0][i] * dLx_dy + 
                                         2. * dphiL_dx[1][i] * dLy_dy + 
                                         dphiL_dx[0][i] * dLy_dx) * k2) * weight_ * djac_;

                rhsVector[2*i  ] += -(unaL_ * phi_[i] * k1 +
                                      (2. * dphi_dx[0][i] * duna_dx + 
                                       dphi_dx[1][i] * duna_dy + 
                                       dphi_dx[1][i] * dvna_dx) * k2) * weight_ * djac_;
                rhsVector[2*i+1] += -(vnaL_ * phi_[i] * k1 +
                                      (dphi_dx[0][i] * duna_dy + 
                                       2. * dphi_dx[1][i] * dvna_dy + 
                                       dphi_dx[0][i] * dvna_dx) * k2) * weight_ * djac_;
            };
        };
        index++;        
    };  
    
    // if (index_ == 3056 ) std::cout << "RHS " << rhsVectorLM << std::endl;

    if (k2 > 0)
        for (int i = 0; i < 18; i++)
            for (int j = 0; j < 18; j++)
                lagrMultMatrix[i][j] += diffMatrix[j][i];

    // LocalVector U_;
    // U_.clear();
    // for (int i = 0; i < 6; ++i)
    // {
    //     U_(2*i  ) = velx(i);
    //     U_(2*i+1) = vely(i);
    // }

    // noalias(rhsVectorLM) = -prod(lagrMultMatrix,U_);

    for (int i = DIM; i--; ) delete [] dphi_dx[i];
    delete [] dphi_dx;
    for (int i = DIM; i--; ) delete [] dphiL_dx[i];
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

    double xsi[2] = {};
    double xsi_intp[2] = {};
    QuadShapeFunction<2,2> shapeQuad;
    int index = 0;
    SpecialQuad sQuad = SpecialQuad();

    //tARLQ_ = -tPSPG2_;

    double phi_[6] = {};
    double phiLM_[6] = {};
    
    double **dphi_dx;
    dphi_dx = new double*[DIM];
    for (int i = DIM; i--; ) dphi_dx[i] = new double[nElNodes];
    double **dphiL_dx;
    dphiL_dx = new double*[DIM];
    for (int i = DIM; i--; ) dphiL_dx[i] = new double[nElNodes];        

    double &dens_ = parameters.getDensity();
    double &alpha_f = parameters.getAlphaF();
    double &k1 = parameters.getArlequinK1();
    double &k2 = parameters.getArlequinK2();
    int &iTimeStep = parameters.getTimeInstant();

    // jacobianNRMatrix.clear();
    // rhsVector.clear();

    double tSUPG_; double tPSPG_; double tLSIC_;

    double **ainv_;
    ainv_ = new double*[DIM];
    for (int i = DIM; i--; ) ainv_[i] = new double[DIM];

       //std::cout << "PSPG Fine " << ielem << " " << tPSPG_ << std::endl;
    for(double *it = sQuad.begin(); it != sQuad.end(); it++){
        
        if ((intPointCorrespElem[index] == ielem)){

            //Defines the integration points adimentional coordinates
            xsi[0] = sQuad.PointList(index,0);
            xsi[1] = sQuad.PointList(index,1);
            
            //Computes the velocity shape functions
            shapeQuad.evaluate(xsi,phi_);
            
            xsi_intp[0] = intPointCorrespXsi[index][0];
            xsi_intp[1] = intPointCorrespXsi[index][1];

            //Computes the coarse mesh shape functions
            shapeQuad.evaluate(xsi_intp,phiLM_);
            
            //Returns the quadrature integration weight
            double weight_ = sQuad.WeightList(index);
            
            //Computes the jacobian matrix
            getJacobianMatrix(xsi_intp, ainv_);
                        
            getSpatialDerivatives(xsi_intp, ainv_, dphi_dx);

            for (int i = 0; i < 2; ++i)
                for (int j = 0; j < 6; ++j)
                    dphiL_dx[i][j] = dphi_dx[i][j];
            // dphiL_dx = dphi_dx;
            // ddphiL_dx = ddphi_dx;

            getJacobianMatrix(xsi, ainv_);
            getSpatialDerivatives(xsi, ainv_, dphi_dx);

            getParameterSUPG(tSUPG_, tPSPG_, tLSIC_, phi_, dphi_dx);

            double una_ = alpha_f * u_ + (1. - alpha_f) * uPrev_;
            double vna_ = alpha_f * v_ + (1. - alpha_f) * vPrev_;

            double duna_dx = alpha_f * du_dx + (1. - alpha_f) * duprev_dx;
            double duna_dy = alpha_f * du_dy + (1. - alpha_f) * duprev_dy;
            double dvna_dx = alpha_f * dv_dx + (1. - alpha_f) * dvprev_dx;
            double dvna_dy = alpha_f * dv_dy + (1. - alpha_f) * dvprev_dy;

            double lxna_ = lagMx_;
            double lyna_ = lagMy_;

            double umeshna_ = alpha_f * umesh_ + (1. - alpha_f) * umeshPrev_;
            double vmeshna_ = alpha_f * vmesh_ + (1. - alpha_f) * vmeshPrev_;

            double wna_ = alpha_f * intPointWeightFunctionSpecial[index] + (1. - alpha_f) * intPointWeightFunctionSpecialPrev[index];

            u_ = 0.;
            v_ = 0.;

            du_dx = 0.;
            du_dy = 0.;
            dv_dx = 0.;
            dv_dy = 0.;

            //Interpolates the velocity components and its spatial derivatives
            for (int i = 0; i < 6; i++){
                u_ += velx[i] * phiLM_[i];
                v_ += vely[i] * phiLM_[i];

                du_dx += velx[i] * dphiL_dx[0][i];
                du_dy += velx[i] * dphiL_dx[1][i];
                dv_dx += vely[i] * dphiL_dx[0][i];
                dv_dy += vely[i] * dphiL_dx[1][i];
            };  
            
            for (int i = 0; i < 6; i++){
                for (int j = 0; j < 6; j++){
                    //SUPG STABILIZATION TERM
                    double LM = 0.;
                    // LM = ((una_ - umeshna_) * dphi_dx[0][j] + (vna_ - vmeshna_) * dphi_dx[0][j])
                    //     * phiLM_(i) * tSUPG_;
                    
                    jacobianNRMatrix[2*j  ][2*i  ] += LM * weight_ * djac_;
                    jacobianNRMatrix[2*j+1][2*i+1] += LM * weight_ * djac_;
                    
                    //PSPG STABILIZATION TERM
                    double Lx = 0.; double Ly = 0.;
                    
                    // Lx = -dphiL_dx(0,i) * phi_[j] * tPSPG_ / dens_ * k1 * wna_;
                    // Ly = -dphiL_dx(1,i) * phi_[j] * tPSPG_ / dens_ * k1 * wna_;
                    
                    jacobianNRMatrix[2*j  ][12+i] += Lx * weight_ * djac_;
                    jacobianNRMatrix[2*j+1][12+i] += Ly * weight_ * djac_;
                };

                //SUPG AND PSPG STABILIZATION TERMS
                double LMx = 0.; double LMy = 0.; double LMp = 0.;
             
                // LMx = ((una_ - umeshna_) * dphiL_dx(0,i) + 
                //        (vna_ - vmeshna_) * dphiL_dx(1,i)) * lxna_ * tSUPG_;
                // LMy = ((una_ - umeshna_) * dphiL_dx(0,i) + 
                //        (vna_ - vmeshna_) * dphiL_dx(1,i)) * lyna_ * tSUPG_;

                // LMp = (dphiL_dx(0,i) * lagMx_ + dphiL_dx(1,i) * lagMy_) * tPSPG_/ dens_;
                
                // if (iTimeStep > 5)
                //     LMp = (dphiL_dx(0,i) * lxna_ + dphiL_dx(1,i) * lyna_) * tPSPG_ / dens_ * k1 +
                //           (dphiL_dx(0,i) * (2. * dLx_dxx + dLx_dxy + dLx_dyy) +
                //            dphiL_dx(1,i) * (2. * dLy_dyy + dLy_dxy + dLy_dxx)) * tPSPG_ / dens_ * k2;
                
                rhsVector[2*i  ] += LMx * weight_ * djac_ * k1;
                rhsVector[2*i+1] += LMy * weight_ * djac_ * k1;
                rhsVector[12+i] += LMp * weight_ * djac_ * k1;
            };
        };
        index++;        
    };  
    
    for (int i = DIM; i--; ) delete [] dphi_dx[i];
    delete [] dphi_dx;
    for (int i = DIM; i--; ) delete [] dphiL_dx[i];
    delete [] dphiL_dx;
    for (int i = DIM; i--; ) delete [] ainv_[i];
    delete [] ainv_;

    return;
};


//------------------------------------------------------------------------------
//----------------------------STEADY LAPLACE PROBEM-----------------------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Element<DIM,DEG>::getLagrangeMultipliersArlequinDifferentMesh(int &ielem, double &tPSPG2_,double* press, double* velx, double* vely,
                                                             double **arlequinStab, double **laplMatrix, double *arlequinStabVector){

    double xsi[2] = {};
    double xsi_intp[2] = {};
    QuadShapeFunction<2,2> shapeQuad;
    int index = 0;
    SpecialQuad sQuad = SpecialQuad();

    //tARLQ_ = -tPSPG2_;

    double phi_[6] = {};
    double phiLM_[6] = {};

    double **dphi_dx;
    dphi_dx = new double*[DIM];
    for (int i = DIM; i--; ) dphi_dx[i] = new double[nElNodes];
    double **dphiL_dx;
    dphiL_dx = new double*[DIM];
    for (int i = DIM; i--; ) dphiL_dx[i] = new double[nElNodes];        
    
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

    double tSUPG_;
    double tPSPG_;
    double tLSIC_;

    getParameterArlequin2();

    for(double* it = sQuad.begin(); it != sQuad.end(); it++){
        
        if ((intPointCorrespElem[index] == ielem)){

            //Defines the integration points adimentional coordinates
            xsi[0] = sQuad.PointList(index,0);
            xsi[1] = sQuad.PointList(index,1);
            
            //Computes the velocity shape functions
            shapeQuad.evaluate(xsi,phi_);
            
            xsi_intp[0] = intPointCorrespXsi[index][0];
            xsi_intp[1] = intPointCorrespXsi[index][1];

            //Computes the coarse mesh shape functions
            shapeQuad.evaluate(xsi_intp,phiLM_);
            
            //Returns the quadrature integration weight
            double weight_ = sQuad.WeightList(index);
            
            //Computes the jacobian matrix
            getJacobianMatrix(xsi_intp, ainv_);
                        
            getSpatialDerivatives(xsi_intp, ainv_, dphi_dx);

            for (int i = 0; i < 2; ++i)
                for (int j = 0; j < 6; ++j)
                    dphiL_dx[i][j] = dphi_dx[i][j];
            // dphiL_dx = dphi_dx;
            // ddphiL_dx = ddphi_dx;

            getJacobianMatrix(xsi, ainv_);
            getSpatialDerivatives(xsi, ainv_, dphi_dx);

            getParameterArlequin(tSUPG_, tPSPG_, tLSIC_, phi_, dphi_dx);

            double una_ = alpha_f * u_ + (1. - alpha_f) * uPrev_;
            double vna_ = alpha_f * v_ + (1. - alpha_f) * vPrev_;

            double duna_dx = alpha_f * du_dx + (1. - alpha_f) * duprev_dx;
            double duna_dy = alpha_f * du_dy + (1. - alpha_f) * duprev_dy;
            double dvna_dx = alpha_f * dv_dx + (1. - alpha_f) * dvprev_dx;
            double dvna_dy = alpha_f * dv_dy + (1. - alpha_f) * dvprev_dy;

            double lxna_ = lagMx_;
            double lyna_ = lagMy_;

            double umeshna_ = alpha_f * umesh_ + (1. - alpha_f) * umeshPrev_;
            double vmeshna_ = alpha_f * vmesh_ + (1. - alpha_f) * vmeshPrev_;

            double wna_ = alpha_f * intPointWeightFunctionSpecial[index] + (1. - alpha_f) * intPointWeightFunctionSpecialPrev[index];

            u_ = 0.;
            v_ = 0.;

            du_dx = 0.;
            du_dy = 0.;
            dv_dx = 0.;
            dv_dy = 0.;

            //Interpolates the velocity components and its spatial derivatives
            for (int i = 0; i < 6; i++){
                u_ += velx[i] * phiLM_[i];
                v_ += vely[i] * phiLM_[i];

                du_dx += velx[i] * dphiL_dx[0][i];
                du_dy += velx[i] * dphiL_dx[1][i];
                dv_dx += vely[i] * dphiL_dx[0][i];
                dv_dy += vely[i] * dphiL_dx[1][i];
            };  

            for (int i = 0; i < 6; i++){
                for (int j = 0; j < 6; j++){     
                    double AM = 0.;
                    double Lpx = 0.; double Lpy = 0.;
                    double LC = 0.; double LL = 0.;

                    // AM = -phiLM_(i) * phi_[j] * 
                    // intPointWeightFunction(index) * tARLQ_;

                    // LL = phiLM_(i) * phi_[j] * tARLQ_ / dens_;
                    LL = (dphi_dx[0][i] * dphi_dx[0][j] + dphi_dx[1][i] * dphi_dx[1][j]) * tARLQ_ / dens_;

                    // Lpx = phi_[i] * dphi_dx[0][i] * intPointWeightFunction(index)
                    //     * tARLQ_ / dens_;
                    // Lpy = phi_[i] * dphi_dx[1][i] * intPointWeightFunction(index)
                    //     * tARLQ_ / dens_;

                    // LC = -phi_[j] * ((una_ - umesh_) * dphi_dx[0][j] + (vna_ - vmesh_) * dphi_dx[1][j]) * phiLM_(i) * tARLQ_ * intPointWeightFunction(index);

                    arlequinStab[2*j  ][2*i  ] += LL * weight_ * djac_;
                    arlequinStab[2*j+1][2*i+1] += LL * weight_ * djac_;
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

                    laplMatrix[2*j  ][2*i  ] += LC * weight_ * djac_;
                    laplMatrix[2*j+1][2*i+1] += LC * weight_ * djac_;

                    laplMatrix[2*j  ][12+i] += Lpx * weight_ * djac_;
                    laplMatrix[2*j+1][12+i] += Lpy * weight_ * djac_;
                    laplMatrix[12+i][2*j  ] += Lpx * weight_ * djac_;
                    laplMatrix[12+i][2*j+1] += Lpy * weight_ * djac_;

                };

                //ARLEQUIN STABILIZATION TERMS
                double Amx = 0.; double Amy = 0.;
                double LCx = 0.; double LCy = 0.;
                double LPx = 0.; double LPy = 0.;
                double LLx = 0.; double LLy = 0.;
                
                // if (iTimeStep > 5){
                //     Amx = -phi_[i] * axm_ * intPointWeightFunctionSpecial(index) * tARLQ_;
                //     Amy = -phi_[i] * aym_ * intPointWeightFunctionSpecial(index) * tARLQ_;
                // }
                //if (iTimeStep < 10){
                    LLx = -(dphi_dx[0][i] * dLx_dx/(1-wna_) + dphi_dx[1][i] * dLx_dy/(1-wna_)) * tARLQ_;
                    LLy = -(dphi_dx[0][i] * dLy_dx/(1-wna_) + dphi_dx[1][i] * dLy_dy/(1-wna_)) * tARLQ_;
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



                arlequinStabVector[2*i  ] += (Amx + LLx) * weight_ * djac_ * (1-wna_) *0;
                arlequinStabVector[2*i+1] += (Amy + LLy) * weight_ * djac_ * (1-wna_) *0;
            };
        };
        index++;        
    };  
    
    for (int i = 0; i < 6; i++){
        if ((*nodes_)[connect_[i]] -> getConstrains(0) == 1) {
            for (int j = 0; j < 18; j++){
                arlequinStab[2*i  ][j] = 0.;
                arlequinStab[j][2*i  ] = 0.;
            };
            //lagrMultMatrix(12+i,12+i) = 1.;
            arlequinStabVector[2*i  ] = 0.0;
            std::cout << "AQUI2 elem different mesh" << std::endl;
        };

        if ((*nodes_)[connect_[i]] -> getConstrains(1) == 1) {
            // std::cout << "aqui" << std::endl;
            for (int j = 0; j < 18; j++){
                arlequinStab[2*i+1][j] = 0.;
                arlequinStab[j][2*i+1] = 0.;
            };
            //lagrMultMatrix(12+i,12+i) = 1.;
            arlequinStabVector[2*i+1] = 0.0;
            std::cout << "AQUI3 elem different mesh" << std::endl;
        };
    };

    for (int i = DIM; i--; ) delete [] dphi_dx[i];
    delete [] dphi_dx;
    for (int i = DIM; i--; ) delete [] dphiL_dx[i];
    delete [] dphiL_dx;
    for (int i = DIM; i--; ) delete [] ainv_[i];
    delete [] ainv_;

    return;
};



#endif

