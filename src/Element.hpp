//------------------------------------------------------------------------------
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
#include "IntegrationQuadrature.hpp"
//#include "IntegrationQuadrature11.hpp"
#include "PartitionedQuadrature.hpp"

/// Defines the fluid element object and all the element information

template<int DIM>
class Element{
 
public:
    /// Defines the class Node locally
    typedef Node<DIM>                                           Nodes;
    typedef IntegQuadratureSpecial<DIM> IntSpecial;
    typedef IntegQuadrature<DIM> IntNormal;


    /// Defines a type to store the element mesh connectivity
    typedef ublas::bounded_vector<int, 4*DIM-2>                 Connectivity;

    /// Defines a blas-type vector with dimension = DIM
    typedef ublas::bounded_vector<double, DIM>                  DimVector;

    /// Defines a blas-type matrix with dimension = DIM x DIM
    typedef ublas::bounded_matrix<double, DIM, DIM>             DimMatrix;

    /// Defines the vector which contains the element nodal coordinates
    typedef ublas::bounded_matrix<double, 4*DIM-2, DIM>         LocalNodes;

    /// Defines the local vector type with dimension 18 for DIM=2 
    /// and 40 for DIM=3
    typedef ublas::bounded_vector<double, 22*DIM-26>            LocalVector;

    /// Defines the local matrix type with dimension 18x18
    /// for DIM=2 and 40x40 for DIM=3
    typedef ublas::bounded_matrix<double, 22*DIM-26, 22*DIM-26> LocalMatrix;

    ///Defines the partitioned integration quadrature rule class locally
    //typedef PartQuadrature<DIM>                                 SpecialQuad;
    typedef IntegQuadratureSpecial<DIM>                         SpecialQuad;

    /// Defines the normal integration quadrature rule class locally
    typedef IntegQuadrature<DIM>                                NormalQuad;

    /// Defines the boundary integration quadrature rule class locally
    typedef BoundaryIntegQuadrature<DIM>                        BoundaryQuad;
    
    /// Define the type VecLocD from class Node locally
    typedef typename Nodes::VecLocD                             VecLoc;

    /// Define the type Integration Point interger variable
    /// from class Partitioned Quadrature
    typedef typename SpecialQuad::PointIntVar                   IntPointInteger;

    /// Define the type Integration Point coordinates variable
    /// from class Partitioned Quadrature
    typedef typename SpecialQuad::PointCoord                    IntPointCoord;

private:
    QuadShapeFunction<DIM> shapeQuad; //Quadratic shape function
    BoundShapeFunction<DIM>shapeBound;//Boundary shape function
    SpecialQuad            sQuad;     //Integration quadrature
    NormalQuad             nQuad;     //Integration quadrature
    std::vector<Nodes *>   nodes_;    //Velocity nodes
    Connectivity  connect_;           //Velocity mesh connectivity 
    int           index_;             //Element index
    LocalNodes    localNodes_;        //Nodal coordinates - velocity
    LocalNodes    localNodesBoundary_;//Nodal coordinates - velocity
    double        djac_;              //Jacobian determinant
    double        weight_;            //Integration point weight
    DimMatrix     ainv_;              //Inverse of Jacobian transform
    double        visc_;              //Fluid dynamic viscosity
    double        dens_;              //Fluid Density
    double        dTime_;             //Time Steps
    double        timeScheme_;        //Time Integration scheme
    LocalMatrix   diffMatrix;         //Diffusion / Viscosity matrix
    LocalMatrix   laplMatrix;         //Laplacian problem matrix
    LocalMatrix   lagrMultMatrix;     //Lagrange Multiplier matrix
    LocalMatrix   jacobianNRMatrix;   //Newton's method jacobian
    LocalMatrix   arlequinStab;       //Newton's method jacobian
    double        x_, y_;
    double        u_, v_, w_, p_;     //Interpolated velocity and pressure
    double        uPrev_, vPrev_, wPrev_, pPrev_;
    double        ax_, ay_, az_;      //Interpolated acceleration
    double        umesh_, vmesh_, wmesh_; //Interpolated mesh velocity
    double        du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, 
                  dw_dx, dw_dy, dw_dz;//Interpolated fluid spatial derivatives
    double        dax_dx, dax_dy, day_dx, day_dy;
    double        du_dxx, du_dyy, du_dxy, dv_dxx, dv_dyy, dv_dxy;
    double        dp_dx, dp_dy, dp_dz;//Interpolated pressure spatial derivative
    double        dp_dxx, dp_dxy, dp_dyx, dp_dyy;
    double        dumesh_dx, dumesh_dy, dumesh_dz,//Interpolated mesh 
                  dvmesh_dx, dvmesh_dy, dvmesh_dz,//velocity derivatives
                  dwmesh_dx, dwmesh_dy, dwmesh_dz;
    double        lagMx_,lagMy_;
    double        dLx_dx, dLx_dy, dLy_dx, dLy_dy;
    double        dLx_dxx, dLx_dxy, dLx_dyy, dLy_dxx, dLy_dxy, dLy_dyy;
    double        tSUPG_;             //SUPG parameter
    double        tPSPG_;             //PSPG parameter
    double        tLSIC_;             //LSIC parameter
    double        tARLQ_;             //Arlequin parameter
    double        tSUGN1_;
    double        tSUGN2_;
    double        tSUGN3_;
    double        hRGN_;
    DimVector     fieldForce_;        //Element field force
    LocalVector   externalForces;     //Field forces integrated over
                                      //the element domain
    LocalVector   rhsVector;          //RHS vector of Newton's method
    LocalVector   rhsVectorLM;        //Lag Mult RHS vector of Newton's method    
    LocalVector   arlequinStabVector;          //RHS vector of Newton's method

    DimVector     xK, XK;
    double        dCk[DIM+1], dck[DIM+1];
    LocalNodes    gradient_;
    int           sideBoundary_;
    double        meshMovingParameter;
    double        pressureDragForce;
    double        pressureLiftForce;
    double        frictionDragForce;
    double        frictionLiftForce;
    double        dragForce;
    double        liftForce;
    bool          arlequin_; //True if is an Arlequin problem
    double        k1;
    double        k2;

    double pi = M_PI;

    int iTimeStep;

    typename NormalQuad::PointWeight   intPointWeightFunction;
    typename SpecialQuad::PointWeight   intPointWeightFunctionSpecial;
    typename NormalQuad::PointWeight   intPointWeightFunctionPrev;
    typename SpecialQuad::PointWeight   intPointDistGlueZone;
    typename SpecialQuad::PointLogical  intPointGlueZone;
    bool          glueZone;
    bool          compressibility;   //True if the element is compressible
    bool          model; //true for local and false for global
    bool          FSIInterface;
    
    std::vector<ublas::bounded_vector<double,DIM> > di;

    //Second derivatives of velocity shape functions
    typename QuadShapeFunction<DIM>::ValueDDeriv ddphi_dx, ddphiL_dx; 
    //First derivatives of velocity shape functions
    typename QuadShapeFunction<DIM>::ValueDeriv  dphi_dx, dphiL_dx;
    //Values of velocity shape functins
    typename QuadShapeFunction<DIM>::Values      phi_;     
    //Values of velocity shape functins
    typename BoundShapeFunction<DIM>::Values     phib_;     
    //Values of velocity shape functins
    typename BoundShapeFunction<DIM>::ValueDeriv dphib_;     

    IntPointCoord                       intPointCorrespXsi;
    IntPointCoord                       intPointCoordinates;
    IntPointInteger                     intPointCorrespElem;

public:
    /// fluid element constructor
    /// @param int element index @param Connectivity element connectivity
    /// @param vector<Nodes> 
    Element(int index, Connectivity& connect, std::vector<Nodes *> nodes){
        index_ = index;
        connect_ = connect;
        nodes_ = nodes;

        localNodes_.clear();  
        djac_ = 0.;             weight_ = 0.;             ainv_.clear();
        visc_ = 0.;             dens_ = 0.;               dTime_ = 0.; 
        timeScheme_ = 0.;       laplMatrix.clear();
        jacobianNRMatrix.clear();lagrMultMatrix.clear();
        u_ = 0.;          v_ = 0.;          w_ = 0;          p_ = 0.; 
        umesh_ = 0.;      vmesh_ = 0.;      wmesh_ = 0.; 
        du_dx = 0.;       du_dy = 0.;       du_dz = 0.; 
        dv_dx = 0.;       dv_dy = 0.;       dv_dz = 0.;
        dw_dx = 0.;       dw_dy = 0.;       dw_dz = 0.;      tSUPG_ = 0.;
        dumesh_dx = 0.;       dumesh_dy = 0.;       dumesh_dz = 0.; 
        dvmesh_dx = 0.;       dvmesh_dy = 0.;       dvmesh_dz = 0.;
        dwmesh_dx = 0.;       dwmesh_dy = 0.;       dwmesh_dz = 0.;
        fieldForce_.clear();     externalForces.clear();     rhsVector.clear();
        rhsVectorLM.clear();
        xK.clear();              XK.clear();                 di.clear();
        glueZone = false;        FSIInterface = false;
        intPointGlueZone.clear();
        sideBoundary_ = 0;
        arlequin_ = false;

        for (int i=0; i < (nQuad.end() - nQuad.begin()); i++){
            intPointWeightFunction(i) = 1.;
            intPointWeightFunctionPrev(i) = 1.;
        };

        for (int i=0; i < (sQuad.end() - sQuad.begin()); i++){
            intPointWeightFunctionSpecial(i) = 1.;
        };

        setLocalNodes();
        getIntegPointCoordinates();

        DimVector xsi;        
        xsi.clear();        
        getJacobianMatrix(xsi);

        sQuad = IntSpecial();
        nQuad = IntNormal();

    };

    //........................Element basic information.........................
    /// Clear all element variables
    void clearVariables();

    /// Sets the element connectivity
    /// @param Connectivity element connectivity
    void setConnectivity(Connectivity& connect){connect_ = connect;};

    /// Gets the element connectivity
    /// @return element connectivity
    Connectivity getConnectivity(){return connect_;};
 
    int getIndex(){return index_;}

    /// Sets the element density
    /// @param double element density
    void setDensity(double& dens){dens_ = dens;}

    /// Sets the element viscosity
    /// @param double element viscosity
    void setViscosity(double& visc){visc_ = visc;}

    /// Sets the time step size
    /// @param double time step size
    void setTimeStep(double& dt){dTime_ = dt;}

    /// Sets the time integration scheme
    /// @param double time integration scheme: 0.0 - Explicit forward Euler;
    /// 1.0 - Implicit backward Euler;
    /// 0.5 - Implicit Trapezoidal Rule.
    void setTimeIntegrationScheme(double& b){timeScheme_ = b;}

    /// Sets the body forces
    /// @param double* body forces
    void setFieldForce(double* ff);

    /// Sets true if is an Arlequin problem
    void setArlequinProblem(){arlequin_ = true;};

    /// Sets the element vector of local nodes
    void setLocalNodes();

    void setTimeStepCounter(int i){iTimeStep = i;};

    /// Sets Arlequin coupling operator constants
    void setArlequinOperatorConstants(double& k11, double& k22){
        k1 = k11;
        k2 = k22;
    };

    /// Compute and store the spatial jacobian matrix
    /// @param bounded_vector integration point adimensional coordinates
    void getJacobianMatrix(ublas::bounded_vector<double, DIM>& xsi);

    /// Compute and store the shape function spatial derivatives
    /// @param bounded_vector integration point adimensional coordinates
    void getSpatialDerivatives(ublas::bounded_vector<double, DIM>& xsi);

    /// Compute and stores the interpolated velocities, mesh velocities, 
    /// previous mesh velocity, acceleration, among others
    void getVelAndDerivatives();

    /// Compute and store the SUPG, PSPG and LSIC stabilization parameters
    void getParameterSUPG();
    double getPSPG(){return tPSPG_;};

    /// Compute the vorticity field
    void computeVorticity();

    /// Gets the element jacobian determinant
    /// @return element jacobinan determinant
    double getJacobian(){return djac_;};

    /// Compute and store the velocity divergent
    void computeVelocityDivergent();

    /// Compute and store the drag and lift forces at the element boundary
    void computeDragAndLiftForces();

    /// Compute and store the drag and lift forces at the element boundary
    double computeSeparationAngle();

    /// Gets the element pressure drag force
    double getPressureDragForce(){return pressureDragForce;};

    /// Gets the element pressure lift force
    double getPressureLiftForce(){return pressureLiftForce;};

    /// Gets the element friction drag force
    double getFrictionDragForce(){return frictionDragForce;};

    /// Gets the element friction lift force
    double getFrictionLiftForce(){return frictionLiftForce;};

    /// Gets the element drag force
    double getDragForce(){return dragForce;};

    /// Gets the element lift force
    double getLiftForce(){return liftForce;};

    /// Compute and store the boundary forces
    ublas::bounded_vector<double,DIM> getBoundaryLoad(DimVector xsi);

    void getBoundaryIntegration();

    /// Gets the spatial jacobian matrix
    /// @param bounded_vector integration point coordinates
    /// @return spatial jacobian matrix
    DimMatrix getJacobianMatrixValues(ublas::bounded_vector<double, DIM>& xsi){
        getJacobianMatrix(xsi);
        return ainv_;};

    /// Gets the nodal gradient value for potential problem
    /// @param int direction @param int element node 
    /// @return nodal gradient value
    double getNodalGradientValue(int dir, int node){return gradient_(node,dir);}

    /// Sets the element side in boundary
    /// @param int side in boundary
    void setElemSideInBoundary(int side){sideBoundary_ = side;};

    /// Gets the element side in boundary
    /// @return side in boundary
    int getElemSideInBoundary(){return sideBoundary_;};

    /// Gets the results in a specific point for the potential problem
    /// @param VecLoc point adimensional coordinates
    /// @return potential results for the specific point
    std::pair<double, ublas::bounded_vector<double,DIM> >
    getPotentialResultsInPoint(VecLoc xsi);

    /// Sets the mesh moving weighting parameter for solving the Laplace problem
    /// @param double parameter value
    void setMeshMovingParameter(double value) {meshMovingParameter = value;};

    /// Gets the mesh moving weighting parameter
    /// @return mesh moving weighting parameter
    double getMeshMovingParameter(){return meshMovingParameter;};

    //.................Element intersection and correspondence..................
    /// Gets the element intersection parameters 
    /// @param minimum coordinates @param maximum coordinates 
    /// @param minimum inner product @param maximum inner product
    /// @param side lenght
    void setIntersectionParameters(DimVector x, DimVector X, 
                                   double *Dk, double *dk,
                           std::vector<ublas::bounded_vector<double,DIM> > dii);

    /// Gets the coordinates intersection parameters
    /// @return minimum and maximum coordinates
    std::pair<DimVector,DimVector> getXIntersectionParameter()
    {return std::make_pair(xK,XK);};

    /// Gets the inner product intersection parameters
    /// @return minimum and maximum inner products
    std::pair<double*,double*> getDIntersectionParameter()
    {return std::make_pair(dck,dCk);};

    /// Gets the side lenght intersection parameter
    /// @return side lenght intersection parameter
    std::vector<ublas::bounded_vector<double,DIM> >getDiIntersectionParameter()
    {return di;};

    //.............................Model functions..............................
    /// Sets if the element is in the gluing zone
    /// @param glueZone: if true is in the glue zone
    void setGlueZone(){glueZone = true;}

    /// Sets which model the fluid element belongs
    /// @param bool model: true = fine; false = coarse.
    void setModel(bool m){model = m;};

    /// Sets if the element belongs to the fluid structure interface
    void setFSIInterface(){FSIInterface = true;};

    /// Sets if the element must be incompressible
    /// If true, the element is compressible else it is incompressible
    void setCompressibility(){compressibility = true;};

    /// Clears the element compressibility condition setting false. 
    /// @see Element::setCompressibility()
    void clearCompressibility(){compressibility = false;};

    /// Gets the element compressibility constrain
    /// @return compressibility constrain @see Element::setCompressibility()
    bool getCompressibility(){return compressibility;};

    //......................Integration Points Information......................
    /// Gets the number of integration points of the special quadrature rule
    /// @retunr number of integration point of the special quadrature rule
    int getNumberOfIntegrationPoints(){return (sQuad.end() - sQuad.begin());};

    /// Sets the integration point correspondence to the overlapped mesh
    /// @param int element correspondent @param VecLoc Adimensional coordinates
    void setIntegrationPointCorrespondence(int ipoint, int elem, VecLoc x){
        intPointCorrespElem(ipoint) = elem;
        intPointCorrespXsi(ipoint,0) = x(0);
        intPointCorrespXsi(ipoint,1) = x(1); };

    /// Gets the integration point correspondence - element
    /// @return overlapped element correspondence
    int getIntegPointCorrespondenceElement(int index)
    {return intPointCorrespElem(index);};

    /// Compute and store the integration points global coordinates
    void getIntegPointCoordinates();

    /// Gets the integration point global coordinates
    /// @param int integration point index @return integration point coordinates
    VecLoc getIntegPointCoordinatesValue(int index){
        VecLoc x;
        x(0) = intPointCoordinates(index,0);
        x(1) = intPointCoordinates(index,1);
        return x;};

    /// Sets the integration point energy weight function
    /// @param int integration point index 
    /// @param double energy weight function value
    void setIntegPointWeightFunction();

    /// Gets the integration point energy weight function
    /// @param int integration point index @return energy weight function value
    double getIntegPointWeightFunction(int index)
    {return intPointWeightFunction(index);};

    /// Sets if the integration point is the gluing zone
    /// @param int integration point index 
    void setIntegPointInGlueZone(int index){intPointGlueZone(index) = true;};

    /// Gets true if the integration point is in the gluing zone
    /// @param int integration point index @return If is in the gluing zone
    bool getIntegPointInGlueZone(int index){return intPointGlueZone(index);};

    /// Sets the integration point signaled distance function
    /// @param int integration point index
    /// @param signaled distance function valye
    void setIntegPointDistFunction(int index, double val)
    {intPointDistGlueZone(index) = val;};

    /// Gets the integration point signaled distance function value
    /// @return integration point signaled distance function value
    double getIntegPointDistFunction(int index)
    {return intPointDistGlueZone(index);};

    //.......................Element vectors and matrices.......................
    /// Compute and store the element matrix for the incompressible flow problem
    /// @param int integration point index
    void getElemMatrix(int index);

    /// Compute and store the element matrix for the Laplace/Poisson problem
    void getElemLaplMatrix();

    /// Sets the boundary conditions for the incompressible flow problem
    void setBoundaryConditions();

    /// Sets the boundary conditions for the Laplace/Poisson problem
    void setBoundaryConditionsLaplace();

    ///Compute and store the residual vector for the incompressible flow problem
    /// @param int integration point index
    void getResidualVector(int index);

    /// Compute and store the residual vector for the Laplace/Poisson problem
    void getResidualVectorLaplace();

    /// Gets the Newton-Raphson's jacobian matrix
    /// @return Newton-Raphson's jacobian matrix
    LocalMatrix getJacNRMatrix(){return jacobianNRMatrix;};

    /// Gets the Newton-Raphson's jacobian matrix
    /// @return Newton-Raphson's jacobian matrix
    LocalMatrix getArlequinStabilizationMatrix(){return arlequinStab;};
    LocalMatrix getArlequinStabilizationMatrix2(){return laplMatrix;};


    /// Gets the Lagrange multiplier operator matrix
    /// @return Lagrange multiplier operator matrix
    LocalMatrix getLagrMultMatrix(){return lagrMultMatrix;};

    /// Gets the residual vector
    /// @return residual vecot
    LocalVector getRhsVector(){return rhsVector;};

    /// Gets the residual vector
    /// @return residual vecot
    LocalVector getArlequinStabilizationVector(){return arlequinStabVector;};   

    /// Gets the Lagrange multiplier residual vector
    /// @return Lagrange multiplier residual vector
    LocalVector getRhsVectorLagMult(){return rhsVectorLM;};

    /// Apply the boundary conditions and returns the matrix and residual vector
    /// for the Lagrange multiplier operator matrix, used when computing the 
    /// operator term from different meshes
    /// @param LocalMatrix Lagrange multiplier operator matrix
    /// @return residual vector and Lagrange multiplier operator matrix with
    /// boundary conditions applied
    std::pair<LocalVector,LocalMatrix> 
    getRhsVectorAndBoundaryConditions(LocalMatrix Ajac){
        jacobianNRMatrix = Ajac;
        setBoundaryConditionsLagrangeMultipliers();       
        return std::make_pair(rhsVector,jacobianNRMatrix);    };


    LocalVector getRhsLagrangeMultipliers(LocalMatrix Ajac){
        LocalVector U_,rhsVectorLM; U_.clear(); rhsVectorLM.clear();

        for (int i=0; i<6; i++){
            U_(2*i  ) = nodes_[connect_(i)] -> getLagrangeMultiplier(0);
            U_(2*i+1) = nodes_[connect_(i)] -> getLagrangeMultiplier(1);
        };
        noalias(rhsVectorLM) = - prod(Ajac,U_);
        return rhsVectorLM;
    }

    LocalVector getRhsVelocities(LocalMatrix Ajac){
        LocalVector U_,rhsVectorLM; U_.clear(); rhsVectorLM.clear();

        for (int i=0; i<6; i++){
            U_(2*i  ) = nodes_[connect_(i)] -> getVelocity(0);
            U_(2*i+1) = nodes_[connect_(i)] -> getVelocity(1);
        };
        noalias(rhsVectorLM) = - prod(Ajac,U_);
        return rhsVectorLM;
    }

    std::pair<LocalVector,LocalMatrix> getBoundaryConditionsLagMult(LocalVector vel, LocalMatrix Ajac){
        for (int i = 0; i < 6; i++){
            if (nodes_[connect_(i)] -> getConstrains(0) == 1) {
                //std::cout << "Aqui1" << std::endl;
                for (int j = 0; j < 18; j++){
                    //lagrMultMatrix(2*i  ,j) = 0.;
                    Ajac(2*i  ,j) = 0.;
                    //Ajac(j,2*i  ) = 0.;
                };
                //lagrMultMatrix(12+i,12+i) = 1.;
                vel(2*i  ) = 0.;
            };

            if (nodes_[connect_(i)] -> getConstrains(1) == 1) {
                //std::cout << "Aqui2" << std::endl;
                for (int j = 0; j < 18; j++){
                    //lagrMultMatrix(2*i+1,j) = 0.;
                    Ajac(2*i+1,j) = 0.;
                    //Ajac(j,2*i+1) = 0.;
                };
                Ajac(12+i,12+i) = 1.;
                vel(2*i+1) = 0.;
            };
            //std::cout << "asda " << connect_(i) << " " << nodes_[connect_(i)] -> getConstrains(1) << " " << vel(2*i+1) << std::endl;

        };
        return std::make_pair(vel,Ajac);
    }

    std::pair<LocalVector,LocalMatrix> getBoundaryConditionsVelocity(LocalVector vel, LocalMatrix Ajac){

        for (int i = 0; i < 6; i++){
            if (nodes_[connect_(i)] -> getConstrains(0) == 1) {
                //std::cout << "Aqui3" << std::endl;
                for (int j = 0; j < 18; j++){
                    //jacobianNRMatrix(2*i  ,j) = 0.;
                    Ajac(j,2*i  ) = 0.;
                    //Ajac(2*i  ,j) = 0.;
                };
                //jacobianNRMatrix(12+i,12+i) = 1.;
                vel(2*i  ) = 0.;
            };

            if (nodes_[connect_(i)] -> getConstrains(1) == 1) {
                //std::cout << "Aqui4" << std::endl;
                for (int j = 0; j < 18; j++){
                    //jacobianNRMatrix(2*i+1,j) = 0.;
                    Ajac(j,2*i+1) = 0.;
                    //Ajac(2*i+1,j) = 0.;                   
                };
                //Ajac(12+i,12+i) = 1.;
                vel(2*i+1) = 0.;
            };
        };


        return std::make_pair(vel,Ajac);
    };


    /// Sets the boundary conditions to the Lagrange multiplier operator
    void setBoundaryConditionsLagrangeMultipliers();

    /// Compute and store the nodal gradient value (for potential problems)
    void computeNodalGradient();

    /// Compute and store the Lagrange multiplier operator when integrating 
    /// the same mesh portion
    void getLagrangeMultipliersSameMesh();

    /// Compute and store the Lagrange multiplier operator when integrationg
    /// the different mesh portion
    /// @param int element of the coarse mesh (used to verify which integration
    /// point belongs to the coarse mesh element)
    void getLagrangeMultipliersDifferentMesh(int ielem,double tPSPG2_,ublas::bounded_vector<double, 6> press, ublas::bounded_vector<double, 6> velx, ublas::bounded_vector<double, 6> vely);

    /// Compute and store the stabilization for the Lagrange multiplier operator
    /// for the same mesh portion
    void getLMStabilizationSameMesh();

    /// Compute and store the stabilization for the Lagrange multiplier operator
    /// for the different mesh portion
    /// @param int element of the coarse mesh 
    /// @see Element::getLagrangeMultipliersDifferentMesh(int ielem)
    void getStabCoarse(int ielem);

    //...............................Problem type...............................
    /// Compute the Transient Navier-Stokes problem matrices and vectors
    void getTransientNavierStokes();

    /// Compute the Steady Laplace problem matrices and vectors 
    /// (usually for the mesh moving step)
    void getSteadyLaplace();

    /// Compute the Transient Laplace problem matrices and vectors
    /// (usually for the mesh moving step)
    void getTransientLaplace();

};


//------------------------------------------------------------------------------
//--------------------------------IMPLEMENTATION--------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//----------------------SET ELEMENT INTERSECTION PARAMETERS---------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::setIntersectionParameters(DimVector x, DimVector X, 
                                           double *Dk, double *dk,
                           std::vector<ublas::bounded_vector<double,2> > dii) {
    xK = x; XK = X; 
    dCk[0] = Dk[0];
    dCk[1] = Dk[1]; 
    dCk[2] = Dk[2];  
    
    dck[0] = dk[0];
    dck[1] = dk[1];
    dck[2] = dk[2];   
    
    di = dii;
    
    return;
};

//------------------------------------------------------------------------------
//----------------------SET ELEMENT INTERSECTION PARAMETERS---------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::setIntegPointWeightFunction() {
    
    typename QuadShapeFunction<2>::Coords xsi;
    
    intPointWeightFunctionPrev = intPointWeightFunction;
    intPointWeightFunction.clear();
    intPointWeightFunctionSpecial.clear();

    int index=0;

    for(typename NormalQuad::QuadratureListIt it = nQuad.begin(); 
        it != nQuad.end(); it++){
        
       xsi(0) = nQuad.PointList(index,0);
       xsi(1) = nQuad.PointList(index,1);
            
       //Computes the velocity shape functions
       shapeQuad.evaluate(xsi,phi_);

       for (int j=0; j<6; j++){
           intPointWeightFunction(index) += phi_(j) * nodes_[connect_(j)] -> 
               getWeightFunction();
       };
       // intPointWeightFunction(index) = 1.;
       index++;
    }; 

    index = 0;

    for(typename SpecialQuad::QuadratureListIt it = sQuad.begin(); 
        it != sQuad.end(); it++){
        
       xsi(0) = sQuad.PointList(index,0);
       xsi(1) = sQuad.PointList(index,1);
            
       //Computes the velocity shape functions
       shapeQuad.evaluate(xsi,phi_);

       for (int j=0; j<6; j++){
           intPointWeightFunctionSpecial(index) += phi_(j) * nodes_[connect_(j)] -> 
               getWeightFunction();
       };
       // intPointWeightFunction(index) = 1.;
       index++;
    }; 

    //     typename Nodes::VecLocD x;
    // int cont = 0;
    // if (model == false){
    //     for (int i = 0; i < 6; i++){
    //         x = nodes_[connect_(i)] -> getCoordinates();
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

    
    return;
};

//------------------------------------------------------------------------------
//----------------------------SET ELEMENT FIELD FORCE---------------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::setFieldForce(double *ff) {
    fieldForce_(0) = ff[0];
    fieldForce_(1) = ff[1];
    
    return;
};

template<>
void Element<3>::setFieldForce(double *ff) {
    fieldForce_(0) = ff[0];
    fieldForce_(1) = ff[1];
    fieldForce_(2) = ff[2];
    
    return;
};

//------------------------------------------------------------------------------
//------------------COMPUTES THE INTEGRATION POINT COORDINATE-------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::getIntegPointCoordinates(){

    typename SpecialQuad::NodalValuesQuad nodalCoordx, nodalCoordy;
    VecLoc x; 
    
    for (int i = 0; i < 6; i++){
        nodalCoordx(i) = localNodes_(i,0);
        nodalCoordy(i) = localNodes_(i,1);
    };

    for (int i = 0; i < sQuad.end() - sQuad.begin(); i++){
        x(0) = 0.; x(1) = 0.;

        x(0) = sQuad.interpolateQuadraticVariable(nodalCoordx, i);
        x(1) = sQuad.interpolateQuadraticVariable(nodalCoordy, i);

        intPointCoordinates(i,0) = x(0);
        intPointCoordinates(i,1) = x(1);
    };

    return;
};

template<>
void Element<3>::getIntegPointCoordinates(){

    typename SpecialQuad::NodalValuesQuad nodalCoordx, nodalCoordy, nodalCoordz;
    VecLoc x; 
    
    for (int i = 0; i < 10; i++){
        nodalCoordx(i) = localNodes_(i,0);
        nodalCoordy(i) = localNodes_(i,1);
    };

    for (int i = 0; i < sQuad.end() - sQuad.begin(); i++){
        x(0) = 0.; x(1) = 0.; x(2) = 0.;

        x(0) = sQuad.interpolateQuadraticVariable(nodalCoordx, i);
        x(1) = sQuad.interpolateQuadraticVariable(nodalCoordy, i);
        x(2) = sQuad.interpolateQuadraticVariable(nodalCoordz, i);

        intPointCoordinates(i,0) = x(0);
        intPointCoordinates(i,1) = x(1);
        intPointCoordinates(i,2) = x(2);
    };

    return;
};

//------------------------------------------------------------------------------
//------------------------SEARCH THE ELEMENT LOCAL NODES------------------------
//------------------------------------------------------------------------------
//Creates a vector with the local nodes coordinates
template<>
void Element<2>::setLocalNodes() {
    
    typename Nodes::VecLocD x;

    for (int i=0; i<6; i++){
        x = nodes_[connect_(i)] -> getCoordinates();
        localNodes_(i,0) = x(0);
        localNodes_(i,1) = x(1);
    };

    getIntegPointCoordinates();

    return;
};

template<>
void Element<3>::setLocalNodes() {
    
    typename Nodes::VecLocD x;

    for (int i=0; i<10; i++){
        x = nodes_[connect_(i)] -> getCoordinates();
        localNodes_(i,0) = x(0);
        localNodes_(i,1) = x(1);
        localNodes_(i,2) = x(2);
    };

    return;
};

//------------------------------------------------------------------------------
//---------------------------CLEAR ELEMENT VARIABLES----------------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::clearVariables(){

    djac_ = 0.;               weight_ = 0.;              ainv_.clear();
    
    diffMatrix.clear();        gradient_.clear();
    laplMatrix.clear();
    jacobianNRMatrix.clear();  lagrMultMatrix.clear();
    
    u_ = 0.;          v_ = 0.;          w_ = 0;          p_ = 0.; 
    umesh_ = 0.;      vmesh_ = 0.;      wmesh_ = 0.; 
    du_dx = 0.;       du_dy = 0.;       du_dz = 0.; 
    dv_dx = 0.;       dv_dy = 0.;       dv_dz = 0.;
    dw_dx = 0.;       dw_dy = 0.;       dw_dz = 0.; 
    tSUPG_ = 0.;      tPSPG_ = .0;      tLSIC_ = 0.;
    
    externalForces.clear();   rhsVector.clear();
    localNodes_.clear();      rhsVectorLM.clear();

    xK.clear();              XK.clear();                 di.clear();

    ddphi_dx.clear();        dphi_dx.clear();            phi_.clear();     

    intPointCorrespXsi.clear();    intPointCoordinates.clear();
    intPointCorrespElem.clear();

    glueZone = false;        compressibility = false;
    intPointGlueZone.clear();
    arlequin_ = false;

    for (int i=0; i < (sQuad.end() - sQuad.begin()); i++){
        intPointWeightFunction(i) = 1.;
    };
    
    setLocalNodes();
    getIntegPointCoordinates();

    return;
}; 

template<>
void Element<3>::clearVariables(){

    djac_ = 0.;               weight_ = 0.;              ainv_.clear();
    
    diffMatrix.clear();        gradient_.clear();
    laplMatrix.clear();
    jacobianNRMatrix.clear();  lagrMultMatrix.clear();
    
    u_ = 0.;          v_ = 0.;          w_ = 0;          p_ = 0.; 
    umesh_ = 0.;      vmesh_ = 0.;      wmesh_ = 0.; 
    du_dx = 0.;       du_dy = 0.;       du_dz = 0.; 
    dv_dx = 0.;       dv_dy = 0.;       dv_dz = 0.;
    dw_dx = 0.;       dw_dy = 0.;       dw_dz = 0.; 
    tSUPG_ = 0.;      tPSPG_ = 0.;      tLSIC_ = 0.;
    
    externalForces.clear();   rhsVector.clear();  rhsVectorLM.clear();
    
    glueZone = false;
    intPointGlueZone.clear();
    
    for (int i=0; i < (sQuad.end() - sQuad.begin()); i++){
        intPointWeightFunction(i) = 1.;
    };

    setLocalNodes();
    getIntegPointCoordinates();

    return;
}; 

//------------------------------------------------------------------------------
//-------------------------SPATIAL TRANSFORM - JACOBIAN-------------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::getJacobianMatrix(ublas::bounded_vector<double,2>& xsi) {

    //Computes the spatial Jacobian matrix and its inverse

    typename QuadShapeFunction<2>::ValueDeriv dphi;
    
    shapeQuad.evaluateGradient(xsi,dphi);

    double dx_dxsi1 = 0.;
    double dx_dxsi2 = 0.;
    double dy_dxsi1 = 0.;
    double dy_dxsi2 = 0.;

    for (int i=0; i<6; i++){
        dx_dxsi1 += localNodes_(i,0) * dphi(0,i);
        dx_dxsi2 += localNodes_(i,0) * dphi(1,i);
        dy_dxsi1 += localNodes_(i,1) * dphi(0,i);
        dy_dxsi2 += localNodes_(i,1) * dphi(1,i);        
    };

    //Computing the jacobian determinant
    djac_ = dx_dxsi1 * dy_dxsi2 - dx_dxsi2 * dy_dxsi1;
    //djac_ = fabs(djac_);

    //Computing Jacobian inverse
    ainv_(0,0) =  dy_dxsi2 / djac_;
    ainv_(0,1) = -dy_dxsi1 / djac_;
    ainv_(1,0) = -dx_dxsi2 / djac_;
    ainv_(1,1) =  dx_dxsi1 / djac_;

    return;
    
};

template<>
void Element<3>::getJacobianMatrix(ublas::bounded_vector<double,3>& xsi) {

    //Computes the spatial Jacobian matrix and its inverse

    typename QuadShapeFunction<3>::ValueDeriv dphi;
    
    shapeQuad.evaluateGradient(xsi,dphi);

    double dx_dxsi1 = 0.;
    double dx_dxsi2 = 0.;
    double dx_dxsi3 = 0.;
    double dy_dxsi1 = 0.;
    double dy_dxsi2 = 0.;
    double dy_dxsi3 = 0.;
    double dz_dxsi1 = 0.;
    double dz_dxsi2 = 0.;
    double dz_dxsi3 = 0.;

    for (int i=0; i<10; i++){
        dx_dxsi1 += localNodes_(i,0) * dphi(0,i);
        dx_dxsi2 += localNodes_(i,0) * dphi(1,i);
        dx_dxsi3 += localNodes_(i,0) * dphi(2,i);
        dy_dxsi1 += localNodes_(i,1) * dphi(0,i);
        dy_dxsi2 += localNodes_(i,1) * dphi(1,i);        
        dy_dxsi3 += localNodes_(i,1) * dphi(2,i);       
        dz_dxsi1 += localNodes_(i,2) * dphi(0,i);
        dz_dxsi2 += localNodes_(i,2) * dphi(1,i);        
        dz_dxsi3 += localNodes_(i,2) * dphi(2,i);        
    };
    
    //Computing the jacobian determinant
    djac_ = dx_dxsi1 * dy_dxsi2 * dz_dxsi3 - dx_dxsi1 * dy_dxsi3 * dz_dxsi2 
          + dx_dxsi2 * dy_dxsi3 * dz_dxsi1 - dx_dxsi2 * dy_dxsi1 * dz_dxsi3 
          + dx_dxsi3 * dy_dxsi1 * dz_dxsi2 - dx_dxsi3 * dy_dxsi2 * dz_dxsi1;

    //Computing Jacobian inverse
    ainv_(0,0) = (dy_dxsi2 * dz_dxsi3 - dy_dxsi3 * dz_dxsi2) / djac_;
    ainv_(0,1) = (dx_dxsi3 * dz_dxsi2 - dx_dxsi2 * dz_dxsi3) / djac_;
    ainv_(0,2) = (dx_dxsi2 * dy_dxsi3 - dx_dxsi3 * dy_dxsi2) / djac_;
    ainv_(1,0) = (dy_dxsi3 * dz_dxsi1 - dy_dxsi1 * dz_dxsi3) / djac_;
    ainv_(1,1) = (dx_dxsi1 * dz_dxsi3 - dx_dxsi3 * dz_dxsi1) / djac_;
    ainv_(1,2) = (dx_dxsi3 * dy_dxsi1 - dx_dxsi1 * dy_dxsi3) / djac_;
    ainv_(2,0) = (dy_dxsi1 * dz_dxsi2 - dy_dxsi2 * dz_dxsi1) / djac_;
    ainv_(2,1) = (dx_dxsi2 * dz_dxsi1 - dx_dxsi1 * dz_dxsi2) / djac_;
    ainv_(2,2) = (dx_dxsi1 * dy_dxsi2 - dx_dxsi2 * dy_dxsi1) / djac_;

    return;
    
};

//------------------------------------------------------------------------------
//-----------------------------SPATIAL DERIVATIVES------------------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::getSpatialDerivatives(ublas::bounded_vector<double,2>& xsi) {
    
    typename QuadShapeFunction<2>::ValueDDeriv ddphi;
    typename QuadShapeFunction<2>::ValueDeriv  dphi;
    
    shapeQuad.evaluateGradient(xsi,dphi);
    shapeQuad.evaluateHessian(xsi,ddphi);
    
    dphi_dx.clear();
    ddphi_dx.clear();

    //Quadratic shape functions spatial first derivatives
    noalias(dphi_dx) = prod(ainv_,dphi);

    //Quadratic shape functions spatial second derivatives
    // for (int i=0; i<6; i++){
    //     for (int j=0; i<2; i++){
    //         for (int k=0; i<2; i++){
    //             for (int l=0; i<2; i++){
    //                 for (int m=0; i<2; i++){                 
    //                     ddphi_dx(j,l)(i) += ainv_(j,k) * ainv_(l,m) * 
    //                                         ddphi(k,m)(i);
    //                 };
    //             };
    //         };
    //     };
    // };


    // Computing second derivatives
    
    double dx_dxsi1 = 0.;
    double dx_dxsi2 = 0.;
    double dy_dxsi1 = 0.;
    double dy_dxsi2 = 0.;

    double dx_dxsi11 = 0.;
    double dx_dxsi22 = 0.;
    double dx_dxsi12 = 0.;
    double dy_dxsi11 = 0.;
    double dy_dxsi22 = 0.;
    double dy_dxsi12 = 0.;

    for (int i=0; i<6; i++){
        dx_dxsi1 += localNodes_(i,0) * dphi(0,i);
        dx_dxsi2 += localNodes_(i,0) * dphi(1,i);
        dy_dxsi1 += localNodes_(i,1) * dphi(0,i);
        dy_dxsi2 += localNodes_(i,1) * dphi(1,i);

        dx_dxsi11 += localNodes_(i,0) * ddphi(0,0)(i);
        dx_dxsi22 += localNodes_(i,0) * ddphi(1,1)(i);
        dx_dxsi12 += localNodes_(i,0) * ddphi(0,1)(i);
        dy_dxsi11 += localNodes_(i,1) * ddphi(0,0)(i);
        dy_dxsi22 += localNodes_(i,1) * ddphi(1,1)(i);
        dy_dxsi12 += localNodes_(i,1) * ddphi(0,1)(i);        
    };        

    ublas::bounded_matrix<double, 3, 3> invM;
    ublas::bounded_vector<double, 3> vecM, resM;
    
    double a = dx_dxsi1 * dx_dxsi1;
    double b = dy_dxsi1 * dy_dxsi1;
    double c = 2. * dx_dxsi1 * dy_dxsi1;
    double d = dx_dxsi2 * dx_dxsi2;
    double e = dy_dxsi2 * dy_dxsi2;
    double f = 2. * dx_dxsi2 * dy_dxsi2;
    double g = dx_dxsi1 * dx_dxsi2;
    double h = dy_dxsi1 * dy_dxsi2;
    double i = dx_dxsi1 * dy_dxsi2 + dx_dxsi2 * dy_dxsi1;
    
    // matrix:
    // |a b c|
    // |d e f|
    // |g h i|
    
    double det = 1. / (a * e * i - a * f * h - b * d * i + 
                       b * f * g + c * d * h - c * e * g);
    
    invM(0,0) = (e * i - f * h) * det;
    invM(0,1) = (c * h - b * i) * det;
    invM(0,2) = (b * f - c * e) * det;
    invM(1,0) = (f * g - d * i) * det;
    invM(1,1) = (a * i - c * g) * det;
    invM(1,2) = (c * d - a * f) * det;
    invM(2,0) = (d * h - e * g) * det;
    invM(2,1) = (b * g - a * h) * det;
    invM(2,2) = (a * e - b * d) * det;
    
    for (int j = 0; j < 6; j++){
        
        vecM(0) = ddphi(0,0)(j)
            - dphi_dx(0,j) * dx_dxsi11 - dphi_dx(1,j) * dy_dxsi11;
        vecM(1) = ddphi(1,1)(j)
            - dphi_dx(0,j) * dx_dxsi22 - dphi_dx(1,j) * dy_dxsi22;
        vecM(2) = ddphi(0,1)(j)
            - dphi_dx(0,j) * dx_dxsi12 - dphi_dx(1,j) * dy_dxsi12;

        resM = prod(invM,vecM);
        
        ddphi_dx(0,0)(j) = resM(0);
        ddphi_dx(1,1)(j) = resM(1);
        ddphi_dx(0,1)(j) = resM(2);
        ddphi_dx(1,0)(j) = resM(2);
  
        // if(index_ == 0) std::cout << "SADASD " << vecM(0) << " " << vecM(1) << " " << vecM(2) << " " << det << " " << e << " " << f << " " << g << " " << h << " " << i << std::endl;
      
    };

  
    return;
};

template<>
void Element<3>::getSpatialDerivatives(ublas::bounded_vector<double,3>& xsi) {

    typename QuadShapeFunction<3>::ValueDDeriv ddphi;    
    typename QuadShapeFunction<3>::ValueDeriv  dphi;

    shapeQuad.evaluateGradient(xsi,dphi);
    shapeQuad.evaluateHessian(xsi,ddphi);

    dphi_dx.clear();
    ddphi_dx.clear();

    //Quadratic shape functions spatial first derivatives
    for (int i=0; i<10; i++){
        dphi_dx(0,i) += ainv_(0,0) * dphi(0,i) 
                      + ainv_(1,0) * dphi(1,i) 
                      + ainv_(2,0) * dphi(2,i);
 
        dphi_dx(1,i) += ainv_(0,1) * dphi(0,i) 
                      + ainv_(1,1) * dphi(1,i) 
                      + ainv_(2,1) * dphi(2,i);         
 
       dphi_dx(2,i) += ainv_(0,2) * dphi(0,i) 
                      + ainv_(1,2) * dphi(1,i) 
                      + ainv_(2,2) * dphi(2,i);         
    };

    //Quadratic shape functions spatial second derivatives
    for (int i=0; i<10; i++){
        for (int j=0; i<3; i++){
            for (int k=0; i<3; i++){
                for (int l=0; i<3; i++){
                    for (int m=0; i<3; i++){                 
                        ddphi_dx(j,l)(i) += ainv_(j,k) * ainv_(l,m) * 
                                            ddphi(k,m)(i);
                    };
                };
            };
        };
    };

    return;
};

//------------------------------------------------------------------------------
//-------------INTERPOLATES VELOCITY, PRESSURE AND ITS DERIVATIVES--------------
//------------------------------------------------------------------------------
template<>
void Element<2>::getVelAndDerivatives() {

    u_ = 0.; v_ = 0.;
    x_ = 0.; y_ = 0.;
    ax_ = 0.; ay_ = 0.;
    uPrev_ = 0.; vPrev_ = 0.; pPrev_ = 0.; 
    umesh_ = 0.; vmesh_ = 0.;
    p_ = 0.;
    du_dx = 0.; du_dy = 0.; dv_dx = 0.; dv_dy = 0.;
    du_dxx = 0.; du_dyy = 0.; du_dxy = 0.; dv_dxx = 0.; dv_dyy = 0.; dv_dxy = 0.;
    dp_dx = 0.; dp_dy = 0.;
    lagMx_ = 0.; lagMy_ = 0.;
    dLx_dx = 0.; dLx_dy = 0.; dLy_dx = 0.; dLy_dy = 0.;
    dLx_dxx = 0.; dLx_dxy = 0.; dLx_dyy = 0.; dLy_dxx = 0.; dLy_dxy = 0.; dLy_dyy = 0.;
    dp_dxx = 0.; dp_dxy = 0.; dp_dyx = 0.; dp_dyy = 0.;
    dax_dx = 0.; dax_dy = 0.; day_dx = 0.; day_dy = 0.;

    //Interpolates the velocity components and its spatial derivatives
    for (int i = 0; i < 6; i++){
        u_ += nodes_[connect_(i)] -> getVelocity(0) * phi_(i);
        v_ += nodes_[connect_(i)] -> getVelocity(1) * phi_(i);

        x_ += localNodes_(i,0) * phi_(i);
        y_ += localNodes_(i,1) * phi_(i);

        uPrev_ += nodes_[connect_(i)] -> getPreviousVelocity(0) * phi_(i);
        vPrev_ += nodes_[connect_(i)] -> getPreviousVelocity(1) * phi_(i);

        // ax_ += nodes_[connect_(i)] -> getAcceleration(0) * phi_(i);
        // ay_ += nodes_[connect_(i)] -> getAcceleration(1) * phi_(i);
        
        du_dx += nodes_[connect_(i)] -> getVelocity(0) * dphi_dx(0,i);
        du_dy += nodes_[connect_(i)] -> getVelocity(0) * dphi_dx(1,i);
        dv_dx += nodes_[connect_(i)] -> getVelocity(1) * dphi_dx(0,i);
        dv_dy += nodes_[connect_(i)] -> getVelocity(1) * dphi_dx(1,i);

        du_dxx += nodes_[connect_(i)] -> getVelocity(0) * ddphi_dx(0,0)(i);
        du_dyy += nodes_[connect_(i)] -> getVelocity(0) * ddphi_dx(1,1)(i);
        du_dxy += nodes_[connect_(i)] -> getVelocity(0) * ddphi_dx(0,1)(i);
        dv_dxx += nodes_[connect_(i)] -> getVelocity(1) * ddphi_dx(0,0)(i);
        dv_dyy += nodes_[connect_(i)] -> getVelocity(1) * ddphi_dx(1,1)(i);
        dv_dxy += nodes_[connect_(i)] -> getVelocity(1) * ddphi_dx(0,1)(i);

        p_ += nodes_[connect_(i)] -> getPressure() * phi_(i);
        
        pPrev_ += nodes_[connect_(i)] -> getPreviousPressure() * phi_(i);

        dp_dx += nodes_[connect_(i)] -> getPressure() * dphi_dx(0,i);
        dp_dy += nodes_[connect_(i)] -> getPressure() * dphi_dx(1,i);

        dp_dxx += nodes_[connect_(i)] -> getPressure() * ddphi_dx(0,0)(i);
        dp_dxy += nodes_[connect_(i)] -> getPressure() * ddphi_dx(0,1)(i);
        dp_dyx += nodes_[connect_(i)] -> getPressure() * ddphi_dx(1,0)(i);
        dp_dyy += nodes_[connect_(i)] -> getPressure() * ddphi_dx(1,1)(i);

        umesh_ += nodes_[connect_(i)] -> getMeshVelocity(0) * phi_(i);
        vmesh_ += nodes_[connect_(i)] -> getMeshVelocity(1) * phi_(i);
        
        dumesh_dx += nodes_[connect_(i)] -> getMeshVelocity(0) * dphi_dx(0,i);
        dumesh_dy += nodes_[connect_(i)] -> getMeshVelocity(0) * dphi_dx(1,i);
        dvmesh_dx += nodes_[connect_(i)] -> getMeshVelocity(1) * dphi_dx(0,i);
        dvmesh_dy += nodes_[connect_(i)] -> getMeshVelocity(1) * dphi_dx(1,i);
        // std::cout << "velocity " << phi_(i) << " " << nodes_[connect_(i)] -> getVelocity(0) << std::endl;

        lagMx_ += nodes_[connect_(i)] -> getLagrangeMultiplier(0) * phi_(i);
        lagMy_ += nodes_[connect_(i)] -> getLagrangeMultiplier(1) * phi_(i);

        dLx_dxx += nodes_[connect_(i)] -> getLagrangeMultiplier(0) * ddphi_dx(0,0)(i);
        dLx_dyy += nodes_[connect_(i)] -> getLagrangeMultiplier(0) * ddphi_dx(1,1)(i);
        dLx_dxy += nodes_[connect_(i)] -> getLagrangeMultiplier(0) * ddphi_dx(0,1)(i);
        dLy_dxx += nodes_[connect_(i)] -> getLagrangeMultiplier(1) * ddphi_dx(0,0)(i);
        dLy_dyy += nodes_[connect_(i)] -> getLagrangeMultiplier(1) * ddphi_dx(1,1)(i);
        dLy_dxy += nodes_[connect_(i)] -> getLagrangeMultiplier(1) * ddphi_dx(0,1)(i);        

        dLx_dx += nodes_[connect_(i)] -> getLagrangeMultiplier(0) * dphi_dx(0,i);
        dLx_dy += nodes_[connect_(i)] -> getLagrangeMultiplier(0) * dphi_dx(1,i);
        dLy_dx += nodes_[connect_(i)] -> getLagrangeMultiplier(1) * dphi_dx(0,i);
        dLy_dy += nodes_[connect_(i)] -> getLagrangeMultiplier(1) * dphi_dx(1,i);

        dax_dx += (nodes_[connect_(i)] -> getVelocity(0) - nodes_[connect_(i)] -> getPreviousVelocity(0)) * dphi_dx(0,i);
        dax_dy += (nodes_[connect_(i)] -> getVelocity(0) - nodes_[connect_(i)] -> getPreviousVelocity(0)) * dphi_dx(1,i);
        day_dx += (nodes_[connect_(i)] -> getVelocity(1) - nodes_[connect_(i)] -> getPreviousVelocity(1)) * dphi_dx(0,i);
        day_dy += (nodes_[connect_(i)] -> getVelocity(1) - nodes_[connect_(i)] -> getPreviousVelocity(1)) * dphi_dx(1,i);
    };  

    //std::cout << "VEL " << umesh_ << " " << vmesh_ << std::endl;

    ax_ = (u_ - uPrev_) / dTime_;
    ay_ = (v_ - vPrev_) / dTime_;

    return;
};

//------------------------------------------------------------------------------
//-------------INTERPOLATES VELOCITY, PRESSURE AND ITS DERIVATIVES--------------
//------------------------------------------------------------------------------
template<>
ublas::bounded_vector<double,2>
Element<2>::getBoundaryLoad(DimVector xsi) {

    setLocalNodes();

    ublas::bounded_vector<double,2>            load;
    DimMatrix                                  shearStress;
    typename QuadShapeFunction<2>::ValueDeriv  dphi;
    DimVector                                  t_vector, n_vector;
    ublas::identity_matrix<double>             ident(2);

    //Computes the shape functions        
    shapeQuad.evaluate(xsi,phi_);

    //Computes the jacobian matrix
    getJacobianMatrix(xsi);
    
    //Computes spatial derivatives
    getSpatialDerivatives(xsi);
    
    shapeQuad.evaluateGradient(xsi,dphi); 

    //Interpolates velocity and its derivatives values
    getVelAndDerivatives();
   
   
    shearStress(0,0) = 2. * visc_ * du_dx;
    shearStress(0,1) = visc_ * (du_dy + dv_dx);
    shearStress(1,0) = visc_ * (du_dy + dv_dx);
    shearStress(1,1) = 2. * visc_ * dv_dy;
    
    t_vector.clear();
    n_vector.clear();

    if (sideBoundary_ == 1){
        for (int i = 0; i < 6; i++){
            t_vector(0) -= dphi(1,i) * localNodes_(i,0);
            t_vector(1) -= dphi(1,i) * localNodes_(i,1);
        };        
    };

    if (sideBoundary_ == 2){
        for (int i = 0; i < 6; i++){
            t_vector(0) += dphi(0,i) * localNodes_(i,0);
            t_vector(1) += dphi(0,i) * localNodes_(i,1);
        };        
    };

    if (sideBoundary_ == 0){
        std::cout << "VERIFICAR VETOR NORMAL - getBoundaryLoad" << std::endl;
    };

    n_vector(0) =  t_vector(1) / norm_2(t_vector);
    n_vector(1) = -t_vector(0) / norm_2(t_vector);

    load.clear();
    load = -p_ * prod(ident,n_vector) + prod(shearStress,n_vector);

    //std::cout << "N Vector " << sideBoundary_ << " " <<  n_vector(0) << " " << n_vector(1) << " " << load(0) << " " << load(1) << " " << p_ << std::endl;
    
    return load;
};

//------------------------------------------------------------------------------
//-------------INTERPOLATES VELOCITY, PRESSURE AND ITS DERIVATIVES--------------
//------------------------------------------------------------------------------
template<>void Element<2>::getBoundaryIntegration() {

    
    ublas::bounded_vector<int, 3> nodesb_,corresp; 
    if(sideBoundary_ == 0){
        nodesb_(0) = connect_(1); corresp(0) = 1;
        nodesb_(1) = connect_(4); corresp(1) = 4;
        nodesb_(2) = connect_(2); corresp(2) = 2; 
        for (int i=0; i<2; i++){
            localNodesBoundary_(0,i) = localNodes_(1,i);
            localNodesBoundary_(1,i) = localNodes_(4,i);
            localNodesBoundary_(2,i) = localNodes_(2,i);
        };
    }else{
        if(sideBoundary_ == 1){
            nodesb_(0) = connect_(2); corresp(0) = 2;
            nodesb_(1) = connect_(5); corresp(1) = 5; 
            nodesb_(2) = connect_(0); corresp(2) = 0; 
            for (int i=0; i<2; i++){
                localNodesBoundary_(0,i) = localNodes_(2,i);
                localNodesBoundary_(1,i) = localNodes_(5,i);
                localNodesBoundary_(2,i) = localNodes_(0,i);
            };
        }else{
            nodesb_(0) = connect_(0); corresp(0) = 0;
            nodesb_(1) = connect_(3); corresp(1) = 3;
            nodesb_(2) = connect_(1); corresp(2) = 1;
            for (int i=0; i<2; i++){
                localNodesBoundary_(0,i) = localNodes_(0,i);
                localNodesBoundary_(1,i) = localNodes_(3,i);
                localNodesBoundary_(2,i) = localNodes_(1,i);
            };
        };        
    };

    BoundaryQuad           bQuad;     //Boundary Integration Quadrature
    std::pair<BoundaryQuad::PointCoord,BoundaryQuad::PointWeight> gaussQuad;
    typename QuadShapeFunction<2>::Coords xsi;
    gaussQuad = bQuad.GaussQuadrature();
    DimVector n_vector;
    
    int index = 0;
    for(typename BoundaryQuad::QuadratureListIt it = bQuad.begin(); 
        it != bQuad.end(); it++){
        
        double xsiB = gaussQuad.first(index);
        double weightB = gaussQuad.second(index);

        if(sideBoundary_ == 2){
            xsi(0) = (-xsiB + 1.) / 2.;
            xsi(1) = 0.;
        };
        if(sideBoundary_ == 1){
            xsi(1) = (xsiB + 1.) / 2.;
            xsi(0) = 0.;
        };
        if(sideBoundary_ == 0){
            xsi(0) = (xsiB + 1.) / 2.;
            xsi(1) = 1. - xsi(0);
        };

        //Computes the velocity shape functions
        shapeQuad.evaluate(xsi,phi_);
        
        //Computes the jacobian matrix
        getJacobianMatrix(xsi);

        //Computes spatial derivatives
        getSpatialDerivatives(xsi);

        //Interpolates velocity and its derivatives values
        getVelAndDerivatives();        

        phib_ = shapeBound.getShapeFunction(gaussQuad.first(index));
        dphib_ = shapeBound.getShapeFunctionDerivative(gaussQuad.first(index));

        double Tx=0.; double Ty = 0.;
        
        for (int i=0; i<3; i++){
            Tx += localNodesBoundary_(i,0) * dphib_(i);
            Ty += localNodesBoundary_(i,1) * dphib_(i);
        };

        double jacb_ = sqrt(Tx*Tx + Ty*Ty);


        n_vector(0) =  Ty / jacb_;
        n_vector(1) = -Tx / jacb_;

        for (int i=0; i<3; i++){
            if ((nodes_[nodesb_(i)] -> getConstrains(0) == 0) && (n_vector(0) < -0.9)){
                if (nodes_[nodesb_(i)] -> getConstrainValue(0) >0 )std::cout << "AQUI " <<  index_ << " " << n_vector(0) << " " << n_vector(1) << " " <<  nodes_[nodesb_(i)] -> getConstrainValue(0) * n_vector(0) * phi_(corresp(i)) * jacb_ * weightB << std::endl;
                rhsVector(2*corresp(i)  ) -= nodes_[nodesb_(i)] -> getConstrainValue(0) * n_vector(0) * phi_(corresp(i)) * jacb_ * weightB;
                rhsVector(2*corresp(i)+1) -= nodes_[nodesb_(i)] -> getConstrainValue(0) * n_vector(1) * phi_(corresp(i)) * jacb_ * weightB;
            }
        }
        //std::cout << "n_vector " << load_pressure(0) << " " << load_pressure(1) << std::endl;
        index++;

    };





    return;
}


//------------------------------------------------------------------------------
//-------------INTERPOLATES VELOCITY, PRESSURE AND ITS DERIVATIVES--------------
//------------------------------------------------------------------------------
template<>
void Element<2>::computeDragAndLiftForces() {

    setLocalNodes();
    
    ublas::bounded_vector<double, 3> nodesb_; 
    if(sideBoundary_ == 0){
        nodesb_(0) = connect_(1); 
        nodesb_(1) = connect_(4); 
        nodesb_(2) = connect_(2); 
        for (int i=0; i<2; i++){
            localNodesBoundary_(0,i) = localNodes_(1,i);
            localNodesBoundary_(1,i) = localNodes_(4,i);
            localNodesBoundary_(2,i) = localNodes_(2,i);
        };
    }else{
        if(sideBoundary_ == 1){
            nodesb_(0) = connect_(2); 
            nodesb_(1) = connect_(5); 
            nodesb_(2) = connect_(0); 
            for (int i=0; i<2; i++){
                localNodesBoundary_(0,i) = localNodes_(2,i);
                localNodesBoundary_(1,i) = localNodes_(5,i);
                localNodesBoundary_(2,i) = localNodes_(0,i);
            };
        }else{
            nodesb_(0) = connect_(0);
            nodesb_(1) = connect_(3);
            nodesb_(2) = connect_(1);
            for (int i=0; i<2; i++){
                localNodesBoundary_(0,i) = localNodes_(0,i);
                localNodesBoundary_(1,i) = localNodes_(3,i);
                localNodesBoundary_(2,i) = localNodes_(1,i);
            };
        };        
    };

    BoundaryQuad           bQuad;     //Boundary Integration Quadrature
    std::pair<BoundaryQuad::PointCoord,BoundaryQuad::PointWeight> gaussQuad;
    DimVector                                  n_vector;
    DimMatrix                                  shearStress;
    ublas::identity_matrix<double>             ident(2);
    ublas::bounded_vector<double,2>            load_friction;
    ublas::bounded_vector<double,2>            load_pressure;
    
    typename QuadShapeFunction<2>::Coords xsi;

    n_vector.clear();
    gaussQuad.first.clear();
    gaussQuad.second.clear();
    shearStress.clear();

    load_friction.clear();
    load_pressure.clear();

    gaussQuad = bQuad.GaussQuadrature();
    
    int index = 0;
    for(typename BoundaryQuad::QuadratureListIt it = bQuad.begin(); 
        it != bQuad.end(); it++){
        
        double xsiB = gaussQuad.first(index);
        double weightB = gaussQuad.second(index);

        if(sideBoundary_ == 2){
            xsi(0) = (-xsiB + 1.) / 2.;
            xsi(1) = 0.;
        };
        if(sideBoundary_ == 1){
            xsi(1) = (xsiB + 1.) / 2.;
            xsi(0) = 0.;
        };
        if(sideBoundary_ == 0){
            xsi(0) = (xsiB + 1.) / 2.;
            xsi(1) = 1. - xsi(0);
        };

        //Computes the velocity shape functions
        shapeQuad.evaluate(xsi,phi_);
        
        //Computes the jacobian matrix
        getJacobianMatrix(xsi);

        //Computes spatial derivatives
        getSpatialDerivatives(xsi);

        //Interpolates velocity and its derivatives values
        getVelAndDerivatives();        

        phib_ = shapeBound.getShapeFunction(gaussQuad.first(index));
        dphib_ = shapeBound.getShapeFunctionDerivative(gaussQuad.first(index));

        double Tx=0.; double Ty = 0.;
        
        for (int i=0; i<3; i++){
            Tx += localNodesBoundary_(i,0) * dphib_(i);
            Ty += localNodesBoundary_(i,1) * dphib_(i);
        };

        double jacb_ = sqrt(Tx*Tx + Ty*Ty);
        
        n_vector(0) =  Ty / jacb_;
        n_vector(1) = -Tx / jacb_;

        shearStress(0,0) = 2. * visc_ * du_dx;
        shearStress(0,1) = visc_ * (du_dy + dv_dx);
        shearStress(1,0) = visc_ * (du_dy + dv_dx);
        shearStress(1,1) = 2. * visc_ * dv_dy;
        load_pressure += -p_ * prod(ident,n_vector) * jacb_ * weightB;
        load_friction += prod(shearStress,n_vector) * jacb_ * weightB;
        //std::cout << "n_vector " << load_pressure(0) << " " << load_pressure(1) << std::endl;
        index++;
    };

    pressureDragForce = load_pressure(0);
    pressureLiftForce = load_pressure(1);
    
    frictionDragForce = load_friction(0);
    frictionLiftForce = load_friction(1);

    dragForce = pressureDragForce + frictionDragForce;
    liftForce = pressureLiftForce + frictionLiftForce;

  
    return;
};



//------------------------------------------------------------------------------
//-------------INTERPOLATES VELOCITY, PRESSURE AND ITS DERIVATIVES--------------
//------------------------------------------------------------------------------
template<>
double Element<2>::computeSeparationAngle() {

    setLocalNodes();

    ublas::bounded_vector<double, 3> nodesb_; 
    if(sideBoundary_ == 0){
        nodesb_(0) = connect_(1); 
        nodesb_(1) = connect_(4); 
        nodesb_(2) = connect_(2); 
        for (int i=0; i<2; i++){
            localNodesBoundary_(0,i) = localNodes_(1,i);
            localNodesBoundary_(1,i) = localNodes_(4,i);
            localNodesBoundary_(2,i) = localNodes_(2,i);
        };
    }else{
        if(sideBoundary_ == 1){
            nodesb_(0) = connect_(2); 
            nodesb_(1) = connect_(5); 
            nodesb_(2) = connect_(0); 
            for (int i=0; i<2; i++){
                localNodesBoundary_(0,i) = localNodes_(2,i);
                localNodesBoundary_(1,i) = localNodes_(5,i);
                localNodesBoundary_(2,i) = localNodes_(0,i);
            };
        }else{
            nodesb_(0) = connect_(0);
            nodesb_(1) = connect_(3);
            nodesb_(2) = connect_(1);
            for (int i=0; i<2; i++){
                localNodesBoundary_(0,i) = localNodes_(0,i);
                localNodesBoundary_(1,i) = localNodes_(3,i);
                localNodesBoundary_(2,i) = localNodes_(1,i);
            };
        };        
    };

    BoundaryQuad           bQuad;     //Boundary Integration Quadrature
    std::pair<BoundaryQuad::PointCoord,BoundaryQuad::PointWeight> gaussQuad;
    DimVector                                  n_vector, t_vector;
    DimMatrix                                  shearStress;
    ublas::identity_matrix<double>             ident(2);
    ublas::bounded_vector<double,2>            load_friction;
    ublas::bounded_vector<double,2>            load_pressure;
    ublas::bounded_vector<double,3>            shear_nodal;
    
    typename QuadShapeFunction<2>::Coords xsi;

    n_vector.clear();
    gaussQuad.first.clear();
    gaussQuad.second.clear();
    shearStress.clear();

    load_friction.clear();
    load_pressure.clear();

    gaussQuad = bQuad.GaussQuadrature();
    
    int index = 0;
    for(int it = 0; it < 3; it++){
        
        double xsiB;// = gaussQuad.first(index);
        double weightB = 1.;//gaussQuad.second(index);

        if (it == 0) xsiB = -1.;
        if (it == 1) xsiB =  0.;
        if (it == 2) xsiB =  1.;

        if(sideBoundary_ == 2){
            xsi(0) = (-xsiB + 1.) / 2.;
            xsi(1) = 0.;
        };
        if(sideBoundary_ == 1){
            xsi(1) = (xsiB + 1.) / 2.;
            xsi(0) = 0.;
        };
        if(sideBoundary_ == 0){
            xsi(0) = (xsiB + 1.) / 2.;
            xsi(1) = 1. - xsi(0);
        };

        //Computes the velocity shape functions
        shapeQuad.evaluate(xsi,phi_);
        
        //Computes the jacobian matrix
        getJacobianMatrix(xsi);

        //Computes spatial derivatives
        getSpatialDerivatives(xsi);

        //Interpolates velocity and its derivatives values
        getVelAndDerivatives();        

        phib_ = shapeBound.getShapeFunction(xsiB);
        dphib_ = shapeBound.getShapeFunctionDerivative(xsiB);

        double Tx = 0.; double Ty = 0.;
        
        for (int i=0; i<3; i++){
            Tx += localNodesBoundary_(i,0) * dphib_(i);
            Ty += localNodesBoundary_(i,1) * dphib_(i);
        };

        t_vector(0) = Tx;
        t_vector(1) = Ty;
        
        double jacb_ = sqrt(Tx*Tx + Ty*Ty);
        
        n_vector(0) =  Ty / jacb_;
        n_vector(1) = -Tx / jacb_;

        shearStress(0,0) = 2. * visc_ * du_dx;
        shearStress(0,1) = visc_ * (du_dy + dv_dx);
        shearStress(1,0) = visc_ * (du_dy + dv_dx);
        shearStress(1,1) = 2. * visc_ * dv_dy;
        
        load_friction = prod(shearStress,n_vector) * jacb_ * weightB;
        
        shear_nodal(it) = inner_prod(load_friction,t_vector);
        //std::cout << "n_vector " << load_pressure(0) << " " << load_pressure(1) << std::endl;
        index++;
    };

    //std::cout << "Shear_nodal " << shear_nodal(0) << " " << shear_nodal(1) << " " << shear_nodal(2) << std::endl; 

    double theta = 0.;

    if ((shear_nodal(0)/shear_nodal(1) < 0) || (shear_nodal(1)/shear_nodal(2) < 0)){
        double a = 0.5 * (shear_nodal(0) + shear_nodal(2) - 2*shear_nodal(1));
        double b = 0.5 * (shear_nodal(2) - shear_nodal(0));
        double c = shear_nodal(1);

        //Bhaskara
        double delta = b * b - 4 * a * c;

        double x1 = (-b + sqrt(delta)) / (2*a);
        double x2 = (-b - sqrt(delta)) / (2*a);

        double xsiB;
        if((x1 > -1.00001) && (x1 < 1.00001)) xsiB = x1;
        if((x2 > -1.00001) && (x2 < 1.00001)) xsiB = x2;

        if(sideBoundary_ == 2){
            xsi(0) = (-xsiB + 1.) / 2.;
            xsi(1) = 0.;
        };
        if(sideBoundary_ == 1){
            xsi(1) = (xsiB + 1.) / 2.;
            xsi(0) = 0.;
        };
        if(sideBoundary_ == 0){
            xsi(0) = (xsiB + 1.) / 2.;
            xsi(1) = 1. - xsi(0);
        };

        //Computes the velocity shape functions
        shapeQuad.evaluate(xsi,phi_);
        
        x_ = 0.;
        y_ = 0.;

        for (int i = 0; i < 6; ++i){
           x_ += localNodes_(i,0) * phi_(i);
           y_ += localNodes_(i,1) * phi_(i);
        }

        theta = acos(-x_/sqrt(x_*x_+y_*y_)) * 180 / pi;
        //if(rank == 0) std::cout << "AQUI SEPARATION ANGLE " << index_ << " " << theta << std::endl;

    }

  
    return theta;
};

//------------------------------------------------------------------------------
//------------------COMPUTES THE SUPG STABILIZATION PARAMETER-------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::getParameterSUPG() {

     DimVector   r, s, sl, rl;

    tSUPG_ = 0.;
    tSUGN1_ = 0.;
    tSUGN2_ = 0.;
    tSUGN3_ = 0.;
    hRGN_ = 0.;
    double hUGN_ = 0.;
    
    r.clear(); s.clear(); rl.clear(); sl.clear();

    double u__ = 0.;
    double v__ = 0.;
    double lx__ = 0.;
    double ly__ = 0.;

    for (int i = 0; i < 6; i++){
        double ua = nodes_[connect_(i)] -> getVelocity(0);
        double va = nodes_[connect_(i)] -> getVelocity(1);

        double uma = nodes_[connect_(i)] -> getMeshVelocity(0);
        double vma = nodes_[connect_(i)] -> getMeshVelocity(1);

        double lxa = nodes_[connect_(i)] -> getLagrangeMultiplier(0);
        double lya = nodes_[connect_(i)] -> getLagrangeMultiplier(1);
            
        ua -= uma;
        va -= vma;
        
        u__ += ua * phi_(i);
        v__ += va * phi_(i);

        lx__ += lxa * phi_(i);
        ly__ += lya * phi_(i);
    };

    double uNorm = sqrt(u__ * u__ + v__ * v__);
    double lNorm = sqrt(lx__ * lx__ + ly__ * ly__);

    if(uNorm > 1.e-10){
        s(0) = u__ / uNorm;
        s(1) = v__ / uNorm;
    }else{
        s(0) = 1. / sqrt(2.);
        s(1) = 1. / sqrt(2.);
    };
    if(lNorm > 1.e-10){
        sl(0) = lx__ / lNorm;
        sl(1) = ly__ / lNorm;
    }else{
        sl(0) = 1. / sqrt(2.);
        sl(1) = 1. / sqrt(2.);
    };


    for (int i = 0; i < 6; i++){
        double ua = nodes_[connect_(i)] -> getVelocity(0);
        double va = nodes_[connect_(i)] -> getVelocity(1);

        double uma = nodes_[connect_(i)] -> getMeshVelocity(0);
        double vma = nodes_[connect_(i)] -> getMeshVelocity(1);

        double lxa = nodes_[connect_(i)] -> getLagrangeMultiplier(0);
        double lya = nodes_[connect_(i)] -> getLagrangeMultiplier(1);

        ua -= uma;
        va -= vma;
        
        r(0) += sqrt(ua * ua + va * va) * dphi_dx(0,i);
        r(1) += sqrt(ua * ua + va * va) * dphi_dx(1,i);

        rl(0) += sqrt(lxa * lxa + lya * lya) * dphi_dx(0,i);
        rl(1) += sqrt(lxa * lxa + lya * lya) * dphi_dx(1,i);
    };

    double rNorm = norm_2(r);
    double rlNorm = norm_2(rl);

    if (rNorm >= 1.e-10){
        r /= rNorm;
    }else{
        r(0) = 1. / sqrt(2.);
        r(1) = 1. / sqrt(2.);
    };

    if (rlNorm >= 1.e-10){
        rl /= rlNorm;
    }else{
        rl(0) = 1. / sqrt(2.);
        rl(1) = 1. / sqrt(2.);
    };
    
    double hUGNL_ = 0.;
    double hRGNL_ = 0.;
    for (int i = 0; i < 6; i++){
        hRGN_ += fabs(r(0) * dphi_dx(0,i) + r(1) * dphi_dx(1,i));
        hUGN_ += fabs(s(0) * dphi_dx(0,i) + s(1) * dphi_dx(1,i));        

        hRGNL_ += fabs(rl(0) * dphi_dx(0,i) + rl(1) * dphi_dx(1,i));
        hUGNL_ += fabs(sl(0) * dphi_dx(0,i) + sl(1) * dphi_dx(1,i));        
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

   
    if (fabs(tSUGN1_) <= 1.e-10) tSUGN1_ = 1.e-10;
    if (fabs(tSUGN3_) <= 1.e-10) tSUGN3_ = 1.e-10;

    if (fabs(tSUGN1L_) <= 1.e-10) tSUGN1L_ = 1.e-10;
    if (fabs(tSUGN3L_) <= 1.e-10) tSUGN3L_ = 1.e-10;

    //if(model == false) std::cout << "SUPG " << tSUGN1_ << " " << tSUGN3_ << std::endl;
    //Computing tSUPG parameter
    tSUPG_ = 1. / sqrt(1. / (tSUGN1_ * tSUGN1_) + 
                       1. / (tSUGN2_ * tSUGN2_) + 
                       1. / (tSUGN3_ * tSUGN3_));
    //tSUPG_ = 0.;

    tARLQ_ = -1. / sqrt(1. / (tSUGN1L_ * tSUGN1L_) + 
                       1. / (tSUGN2_ * tSUGN2_) + 
                       1. / (tSUGN3L_ * tSUGN3L_));


    //if (tSUPG_ > 10) tSUPG_ = 0.;

    tPSPG_ = 1*tSUPG_;
    // if (fabs(lagMx_) > 0){
    //  //   std::cout << "aqe" << std::endl;
    //     tARLQ_ = djac_ * sqrt(u_ * u_ + v_ * v_) / sqrt(lagMx_ * lagMx_ + lagMy_ * lagMy_);//-tSUPG_*1;
    // }else{
    tARLQ_ = -1 * tSUPG_;
    //}
    //tARLQ_ = 0.;

    //tARLQ_ = 0.;
    //tSUPG_ = 0.;
    tLSIC_ = tSUPG_ * uNorm * uNorm;
 

     //tARLQ_ = 0.;
    // tSUPG_ = 0.;

    return;
};

//------------------------------------------------------------------------------
//----------------------ELEMENT DIFFUSION/VISCOSITY MATRIX----------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::getElemMatrix(int index){
    
    //tSUPG_ = 0.;
    // tPSPG_ = 0.;
    //tLSIC_ = 0.;

    double una_ = timeScheme_ * u_ + (1. - timeScheme_) * uPrev_;
    double vna_ = timeScheme_ * v_ + (1. - timeScheme_) * vPrev_;

    for (int i = 0; i < 6; i++){
        for (int j = 0; j < 6; j++){

            double acelx = (una_ - uPrev_) / dTime_;
            double acely = (vna_ - vPrev_) / dTime_;
            double wSUPGi = (una_ - umesh_) * dphi_dx(0,i) +
                (v_ - vmesh_) * dphi_dx(1,i);
            double wSUPGj = (una_ - umesh_) * dphi_dx(0,j) +      
                (v_ - vmesh_) * dphi_dx(1,j);
            
            //Mass matrix (use for both directions)
            double mM =  phi_(i) * phi_(j) * dens_ + 
                wSUPGi * phi_(j) * tSUPG_ * dens_;
            double sMxx = dphi_dx(0,i) * phi_(j) * acelx * tSUPG_ * dens_;
            double sMxy = dphi_dx(1,i) * phi_(j) * acelx * tSUPG_ * dens_;
            double sMyx = dphi_dx(0,i) * phi_(j) * acely * tSUPG_ * dens_;
            double sMyy = dphi_dx(1,i) * phi_(j) * acely * tSUPG_ * dens_;
                       
            //Difusion matrix (viscosity)
            double Kxx = 2. * dphi_dx(0,i) * dphi_dx(0,j) * visc_ + 
                dphi_dx(1,i) * dphi_dx(1,j) * visc_;
            double Kxy = dphi_dx(1,i) * dphi_dx(0,j) * visc_;
            double Kyx = dphi_dx(0,i) * dphi_dx(1,j) * visc_;
            double Kyy = 2. * dphi_dx(1,i) * dphi_dx(1,j) * visc_ + 
                dphi_dx(0,i) * dphi_dx(0,j) * visc_;
            
            //Convection matrix
            double Cxx = (dphi_dx(0,j) * (una_ - umesh_) + 
                          dphi_dx(1,j) * (vna_ - vmesh_)) * phi_(i) * dens_
                + wSUPGi * wSUPGj * tSUPG_ * dens_;
            
            double Cyy  = Cxx;
            
            double Cuu = phi_(i) * du_dx * phi_(j) * dens_ +
                tSUPG_ * wSUPGi * du_dx * phi_(j) * dens_ + 
                tSUPG_ * dphi_dx(0,i) * phi_(j) * 
                (u_ * du_dx + v_ * du_dy) * dens_;
            
            double Cuv = phi_(i) * du_dy * phi_(j) * dens_ +
                tSUPG_ * wSUPGi * du_dy * phi_(j) * dens_ + 
                tSUPG_ * dphi_dx(1,i) * phi_(j) * 
                (u_ * du_dx + v_ * du_dy) * dens_;
            
            double Cvu = phi_(i) * dv_dx * phi_(j) * dens_ +
                tSUPG_ * wSUPGi * dv_dx * phi_(j) * dens_ + 
                tSUPG_ * dphi_dx(0,i) * phi_(j) * 
                (u_ * dv_dx + v_ * dv_dy) * dens_;
            
            double Cvv = phi_(i) * dv_dy * phi_(j) * dens_ +
                tSUPG_ * wSUPGi * dv_dy * phi_(j) * dens_ + 
                tSUPG_ * dphi_dx(1,i) * phi_(j) * 
                (u_ * dv_dx + v_ * dv_dy) * dens_;
            
            double Quu = dphi_dx(0,i) * dp_dx * phi_(j) * tSUPG_;
            double Quv = dphi_dx(1,i) * dp_dx * phi_(j) * tSUPG_;
            double Qvu = dphi_dx(0,i) * dp_dy * phi_(j) * tSUPG_;
            double Qvv = dphi_dx(1,i) * dp_dy * phi_(j) * tSUPG_;
            
            double KLSxx = dphi_dx(0,i) * dphi_dx(0,j) * tLSIC_ * dens_;
            double KLSxy = dphi_dx(0,i) * dphi_dx(1,j) * tLSIC_ * dens_;
            double KLSyx = dphi_dx(1,i) * dphi_dx(0,j) * tLSIC_ * dens_;
            double KLSyy = dphi_dx(1,i) * dphi_dx(1,j) * tLSIC_ * dens_;

            double Sxx = + wSUPGi * (2. * ddphi_dx(0,0)(j) + 
                                     ddphi_dx(1,1)(j)) * tSUPG_ * visc_;
            double Sxy = + wSUPGi * ddphi_dx(0,1)(j)
                * tSUPG_ * visc_;
            double Syx = + wSUPGi * ddphi_dx(0,1)(j)
                * tSUPG_ * visc_;
            double Syy = + wSUPGi * (2. * ddphi_dx(1,1)(j) + 
                                     ddphi_dx(0,0)(j)) * tSUPG_ * visc_;
            
            
            jacobianNRMatrix(2*i  ,2*j  ) += (mM + timeScheme_ * dTime_ * 
                                              (sMxx + Kxx + KLSxx / timeScheme_
                                               + Cxx + Cuu + Quu + Sxx))
                * weight_ * djac_ * intPointWeightFunction(index);

            jacobianNRMatrix(2*i+1,2*j+1) += (mM + timeScheme_ * dTime_ *
                                              (sMyy + Kyy + KLSyy / timeScheme_
                                               + Cyy + Cvv + Qvv + Syy))
                * weight_ * djac_ * intPointWeightFunction(index);

            jacobianNRMatrix(2*i  ,2*j+1) += (timeScheme_ * dTime_ * 
                (Kxy + sMxy + Cuv + Quv + Sxy) + KLSxy * dTime_)
                * weight_ * djac_ * intPointWeightFunction(index);

            jacobianNRMatrix(2*i+1,2*j  ) += (timeScheme_ * dTime_ *
                (Kyx + sMyx + Cvu + Qvu + Syx) + KLSyx * dTime_)
                * weight_ * djac_ * intPointWeightFunction(index); 

            
            //SINAL DA PARCELA QUE MULTIPLICA O TSUPG ESTA
            //COM SINAL TROCADO NA FORMULAÇAO DO TEZDUYAR
            //multipy pressure direction x
            double QSUPGx = - (dphi_dx(0,i) * phi_(j)) + 
                ((dphi_dx(0,i) * (una_ - umesh_) - 
                  dphi_dx(1,i) * (vna_ - vmesh_)) * 
                 dphi_dx(0,j) * tSUPG_);
            //multiply pressure direction y
            double QSUPGy = - (dphi_dx(1,i) * phi_(j)) +
                ((dphi_dx(0,i) * (una_ - umesh_) - 
                  dphi_dx(1,i) * (vna_ - vmesh_)) *   
                 dphi_dx(1,j) * tSUPG_);
            //multiply velocity direction x
            double Qx = dphi_dx(0,i) * phi_(j);
            //multiply velocity direction y
            double Qy = dphi_dx(1,i) * phi_(j);
            
            
            jacobianNRMatrix(12+j,2*i  ) += Qx * dTime_ * weight_ * djac_ 
                * intPointWeightFunction(index);
            jacobianNRMatrix(12+j,2*i+1) += Qy * dTime_ * weight_ * djac_ 
                * intPointWeightFunction(index);
            jacobianNRMatrix(2*i  ,12+j) += QSUPGx * dTime_ * weight_ * djac_
                * intPointWeightFunction(index);
            jacobianNRMatrix(2*i+1,12+j) += QSUPGy * dTime_ * weight_ * djac_
                * intPointWeightFunction(index);

            double Hx = dphi_dx(0,i) * phi_(j) * tPSPG_;
            double Hy = dphi_dx(1,i) * phi_(j) * tPSPG_;
            
            double Gx = dphi_dx(0,i) * wSUPGj * tPSPG_;
            double Gy = dphi_dx(1,i) * wSUPGj * tPSPG_;
            
            double Tx = -(dphi_dx(0,j) * (2. * ddphi_dx(0,0)(i) + 
                                         ddphi_dx(1,1)(i)) + 
                         dphi_dx(1,j) * ddphi_dx(0,1)(i))
                * tPSPG_ * visc_ / dens_;
            double Ty = -(dphi_dx(1,j) * (2. * ddphi_dx(1,1)(i) + 
                                         ddphi_dx(0,0)(i)) +
                         dphi_dx(0,j) * ddphi_dx(0,1)(i))
                * tPSPG_ * visc_ / dens_;

            double Guu = (dphi_dx(0,i) * du_dx * phi_(j) + 
                          dphi_dx(1,i) * dv_dx * phi_(j)) * tPSPG_;
            double Gvv = (dphi_dx(0,i) * du_dy * phi_(j) + 
                          dphi_dx(1,i) * dv_dy * phi_(j)) * tPSPG_;
            
            double Q = (dphi_dx(0,i) * dphi_dx(0,j) + 
                        dphi_dx(1,i) * dphi_dx(1,j)) * tPSPG_ / (dens_);
            
            //Hx = 0; Hy = 0; Gx = 0; Gy = 0; Guu = 0; Gvv = 0;


            jacobianNRMatrix(12+j,2*i  ) += (Hx + (Gx + Guu + Tx) * 
                                              timeScheme_ * dTime_)
                * weight_ * djac_ * intPointWeightFunction(index);
            jacobianNRMatrix(12+j,2*i+1) += (Hy + (Gy + Gvv + Ty) * 
                                              timeScheme_ * dTime_)
                * weight_ * djac_ * intPointWeightFunction(index);
            jacobianNRMatrix(12+j,12+i) += Q * dTime_ * weight_ * djac_
                * intPointWeightFunction(index);
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
void Element<2>::setBoundaryConditions(){

    for (int i = 0; i < 6; i++){
         double w = nodes_[connect_(i)] -> getWeightFunction();

        if ((nodes_[connect_(i)] -> getConstrains(0) == 1) ||
            (nodes_[connect_(i)] -> getConstrains(0) == 3) ||
            (w < 0.011))  {
            for (int j = 0; j < 18; j++){
                jacobianNRMatrix(2*i  ,j) = 0.;
                jacobianNRMatrix(j,2*i  ) = 0.;
            };
            jacobianNRMatrix(2*i  , 2*i  ) = 1.;
            rhsVector(2*i  ) = 0.;
        };

        if ((nodes_[connect_(i)] -> getConstrains(1) == 1) ||
            (nodes_[connect_(i)] -> getConstrains(1) == 3) ||
            (w<0.011)) {
            for (int j = 0; j < 18; j++){
                jacobianNRMatrix(2*i+1,j) = 0.;
                jacobianNRMatrix(j,2*i+1) = 0.;
            };
            jacobianNRMatrix(2*i+1, 2*i+1) = 1.;
            rhsVector(2*i+1) =  0.;
        };

        // double w = nodes_[connect_(i)] -> getWeightFunction();
        // if (w < 0.012){
        //     rhsVector(2*i  ) = 0.;
        //     rhsVector(2*i+1) = 0.;
        // }


    };




    // typename Nodes::VecLocD x;

    // for (int i = 0; i < 6; i++){
    //     x = nodes_[connect_(i)] -> getCoordinates();
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
    //     x = nodes_[connect_(i)] -> getCoordinates();
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
void Element<2>::setBoundaryConditionsLaplace(){

    for (int i = 0; i < 6; i++){

        if (nodes_[connect_(i)] -> getConstrainsLaplace(0) == 1) {
            for (int j = 0; j < 18; j++){
                jacobianNRMatrix(2*i  ,j) = 0.;
                //jacobianNRMatrix(j,2*i  ) = 0.;
            };
            jacobianNRMatrix(2*i  , 2*i  ) = 1.;
            //rhsVector(2*i  ) = 0.;
        };

        if (nodes_[connect_(i)] -> getConstrainsLaplace(1) == 1) {
            for (int j = 0; j < 18; j++){
                jacobianNRMatrix(2*i+1,j) = 0.;
                //jacobianNRMatrix(j,2*i+1) = 0.;
            };
            jacobianNRMatrix(2*i+1, 2*i+1) = 1.;
            //rhsVector(2*i+1) =  0.;
        };
    };
    
    return;
};


//------------------------------------------------------------------------------
//---------------APPLY BOUNDARY CONDITIONS TO LAGRANGE MULTIPLIERS--------------
//------------------------------------------------------------------------------
template<>
void Element<2>::setBoundaryConditionsLagrangeMultipliers(){

    rhsVector.clear();
    rhsVectorLM.clear();
    
    LocalVector U_;

    U_.clear();

    for (int i=0; i<6; i++){
        U_(2*i  ) = nodes_[connect_(i)] -> getVelocity(0);
        U_(2*i+1) = nodes_[connect_(i)] -> getVelocity(1);
    };

    //noalias(rhsVector) = - prod(jacobianNRMatrix,U_); 
    
    for (int i = 0; i < 6; i++){
        if (nodes_[connect_(i)] -> getConstrains(0) == 1) {
            for (int j = 0; j < 18; j++){
                //jacobianNRMatrix(2*i  ,j) = 0.;
                jacobianNRMatrix(j,2*i  ) = 0.;
            };
            //jacobianNRMatrix(12+i,12+i) = 1.;
            rhsVector(2*i  ) = 0.0;
        };

        if (nodes_[connect_(i)] -> getConstrains(1) == 1) {
            // std::cout << "aqui" << std::endl;
            for (int j = 0; j < 18; j++){
                //jacobianNRMatrix(2*i+1,j) = 0.;
                jacobianNRMatrix(j,2*i+1) = 0.;
            };
            jacobianNRMatrix(12+i,12+i) = 1.;
            rhsVector(2*i+1) = 0.0;
        };
    };

    return;
};

//------------------------------------------------------------------------------
//-----------------------------RESIDUAL - RHS VECTOR----------------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::getResidualVector(int index){

    double una_ = timeScheme_ * u_ + (1. - timeScheme_) * uPrev_;
    double vna_ = timeScheme_ * v_ + (1. - timeScheme_) * vPrev_;

    double da_dx = 0.;
    double da_dy = 0.;
    
    for (int i=0; i<6; i++){
        da_dx += nodes_[connect_(i)] -> getWeightFunction() * dphi_dx(0,i);
        da_dy += nodes_[connect_(i)] -> getWeightFunction() * dphi_dx(1,i);
    };
    
    for (int i = 0; i < 6; i++){

        // double rMx = dens_ * ((u_ - uPrev_) / dTime_ + 
        //                       (u_ - umesh_) * du_dx)
        //     - dp_dx + visc_ * (ddu_dx + ddu_dy);

        // double rMy = dens_ * ((v_ - vPrev_) / dTime_ + 
        //                       (v_ - vmesh_) * dv_dx)
        //     - dp_dx + visc_ * (ddv_dx + ddv_dy);
        
        double mx = phi_(i) * (u_ - uPrev_) * dens_ + 
            ((una_ - umesh_) * dphi_dx(0,i) + (vna_ - vmesh_) * dphi_dx(1,i)) * 
            (u_ - uPrev_) * tSUPG_ * dens_;
        double my = phi_(i) * (v_ - vPrev_) * dens_ + 
            ((una_ - umesh_) * dphi_dx(0,i) + (vna_ - vmesh_) * dphi_dx(1,i)) * 
            (v_ - vPrev_) * tSUPG_ * dens_;
        
        double Kx = 2. * dphi_dx(0,i) * du_dx * visc_ + 
            dphi_dx(1,i) * du_dy * visc_ + 
            dphi_dx(1,i) * dv_dx * visc_;
        double Ky= dphi_dx(0,i) * du_dy * visc_ + 
            2. * dphi_dx(1,i) * dv_dy * visc_ + 
            dphi_dx(0,i) * dv_dx * visc_;            

        double KLSx = dphi_dx(0,i) * (du_dx + dv_dy) * tLSIC_ * dens_;
        double KLSy = dphi_dx(1,i) * (du_dx + dv_dy) * tLSIC_ * dens_;

        double Cx = (du_dx * (una_ - umesh_) + 
                     du_dy * (vna_ - vmesh_)) * phi_(i) * dens_ +
            ((una_ - umesh_) * dphi_dx(0,i) + (vna_ - vmesh_) * dphi_dx(1,i)) *
            ((una_ - umesh_) * du_dx + (vna_ - vmesh_) * du_dy) * tSUPG_ *dens_;
        double Cy = (dv_dx * (una_ - umesh_) + 
                     dv_dy * (vna_ - vmesh_)) * phi_(i) * dens_ +
            ((una_ - umesh_) * dphi_dx(0,i) + (vna_ - vmesh_) * dphi_dx(1,i)) *
            ((una_ - umesh_) * dv_dx + (vna_ - vmesh_) * dv_dy) * tSUPG_ *dens_;
 
        double Px = - (dphi_dx(0,i) * p_) + ((dphi_dx(0,i) * (una_ - umesh_) +
                                              dphi_dx(1,i) * (vna_ - vmesh_))
                                             * dp_dx * tSUPG_);
        double Py = - (dphi_dx(1,i) * p_) + ((dphi_dx(0,i) * (una_ - umesh_) +
                                              dphi_dx(1,i) * (vna_ - vmesh_))
                                             * dp_dy * tSUPG_);
           
        double Q = ((du_dx + dv_dy) * phi_(i)) +
            (dphi_dx(0,i) * dp_dx + dphi_dx(1,i) * dp_dy) * tPSPG_ / dens_ +
            dphi_dx(0,i) * (u_ - uPrev_) / dTime_ * tPSPG_ +
            dphi_dx(1,i) * (v_ - vPrev_) / dTime_ * tPSPG_ +
            dphi_dx(0,i) * ((una_ - umesh_) * du_dx +
                            (vna_ - vmesh_) * du_dy) * tPSPG_ +
            dphi_dx(1,i) * ((una_ - umesh_) * dv_dx +
                            (vna_ - vmesh_) * dv_dy) * tPSPG_;

        double T = -(dphi_dx(0,i) * (2. * du_dxx + du_dyy + dv_dxy) +
                    dphi_dx(1,i) * (2. * dv_dyy + dv_dxx + du_dxy))
            * tPSPG_ * visc_ / dens_;

        double dAx = 0.;
        double dAy = 0.;
        if (model){
            //Esse
            // dAx = (umesh_ * da_dx + vmesh_ * da_dy) * u_ * phi_(i);
            // dAy = (umesh_ * da_dx + vmesh_ * da_dy) * v_ * phi_(i);



            dAx = (umesh_ * u_ + vmesh_ * v_) * da_dx * phi_(i);
            dAy = (umesh_ * u_ + vmesh_ * v_) * da_dy * phi_(i);
            //std::cout << "AQUI " << da_dx << " " << dAx << " " << da_dy << std::endl;
        } else {
            mx -= u_ * (intPointWeightFunction(index) - 
                        intPointWeightFunctionPrev(index)) / dTime_ * phi_(i);
            my -= v_ * (intPointWeightFunction(index) - 
                        intPointWeightFunctionPrev(index)) / dTime_ * phi_(i);
              // std::cout << "AQUI " << (intPointWeightFunction(index) - 
              //                   intPointWeightFunctionPrev(index)) / dTime_ << " " << v_ << std::endl;
        };

        double Sx = - ((2. * du_dxx + du_dyy + dv_dxy) * 
                       ((una_ - umesh_) * dphi_dx(0,i) + 
                        (vna_ - vmesh_) * dphi_dx(1,i)))
            * tSUPG_ * visc_;
        
        double Sy = - ((2. * dv_dyy + dv_dxx + du_dxy) * 
                     ((una_ - umesh_) * dphi_dx(0,i) + 
                      (vna_ - vmesh_) * dphi_dx(1,i)))
            * tSUPG_ * visc_;
        
        // if (index_ == 10) std::cout << "Nx " << std::scientific << du_dxx << std::endl;

        rhsVector(2*i  ) += (-mx + (-Kx - Px - Cx - dAx - Sx) * dTime_ 
                             - KLSx * dTime_)
            * weight_ * djac_ * intPointWeightFunction(index);
        rhsVector(2*i+1) += (-my + (-Ky - Py - Cy - dAy - Sy) * dTime_ 
                             - KLSy * dTime_)
            * weight_ * djac_ * intPointWeightFunction(index);
        rhsVector(12+i) += (-Q - T) * dTime_ 
            * weight_ * djac_ * intPointWeightFunction(index);
                             
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
void Element<2>::getResidualVectorLaplace(){

    rhsVector.clear();
    
    LocalVector U_;
    typename Nodes::VecLocD x_up;

    U_.clear();

    for (int i=0; i<6; i++){

        x_up = nodes_[connect_(i)] -> getUpdatedCoordinates();        
        
        if (nodes_[connect_(i)] -> getConstrainsLaplace(0) == 1){
            U_(2*i  ) = x_up(0) - localNodes_(i,0);
        };
        if (nodes_[connect_(i)] -> getConstrainsLaplace(1) == 1){
            U_(2*i+1) = x_up(1) - localNodes_(i,1);
        };

        //if(connect_(i) == 29) std::cout << "Element 29 " << x_up(0) << " " << x_up(1) <<std::endl;
    };
    
    noalias(rhsVector) += U_;//- prod(laplMatrix,U_); 

    return;
};

//------------------------------------------------------------------------------
//---------------------------ELEMENT LAPLACIAN MATRIX---------------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::getElemLaplMatrix(){

     for (int i = 0; i < 6; i++){
        for (int j = 0; j < 6; j++){        
            laplMatrix(2*i  ,2*j  ) += (dphi_dx(0,i) * dphi_dx(0,j) +
                                        dphi_dx(1,i) * dphi_dx(1,j)) 
                * weight_ * djac_ * meshMovingParameter;
            laplMatrix(2*i+1,2*j+1) += (dphi_dx(0,i) * dphi_dx(0,j) +
                                        dphi_dx(1,i) * dphi_dx(1,j)) 
                * weight_ * djac_ * meshMovingParameter;
        };
    };
     
    return;
};

//------------------------------------------------------------------------------
//-------------------------COMPUTES NODAL GRADIENT VALUE------------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::computeNodalGradient(){
    
    double dp_dx, dp_dy;
    typename QuadShapeFunction<2>::Coords xsi;
    double potential[6];

    for (int i=0; i<6; i++){
        potential[i] = nodes_[connect_(i)] -> getPotential();
    };
    
    //Node 0
    xsi(0) = 0.;
    xsi(1) = 0.;
        
    getJacobianMatrix(xsi);
    getSpatialDerivatives(xsi);

    dp_dx = 0.; dp_dy = 0.;
    
    for (int i=0; i<6; i++){
        dp_dx += dphi_dx(0,i) * potential[i];
        dp_dy += dphi_dx(1,i) * potential[i];
    };

    nodes_[connect_(0)] -> setGradientComponent( dp_dy, 0);
    nodes_[connect_(0)] -> setGradientComponent(-dp_dx, 1);

    gradient_(0,0) = dp_dy;
    gradient_(0,1) =-dp_dx;

    //Node 1
    xsi(0) = 1.;
    xsi(1) = 0.;

    getJacobianMatrix(xsi);
    getSpatialDerivatives(xsi);

    dp_dx = 0.; dp_dy = 0.;
    
    for (int i=0; i<6; i++){
        dp_dx += dphi_dx(0,i) * potential[i];
        dp_dy += dphi_dx(1,i) * potential[i];
    };

    nodes_[connect_(1)] -> setGradientComponent( dp_dy, 0);
    nodes_[connect_(1)] -> setGradientComponent(-dp_dx, 1);      

    gradient_(1,0) = dp_dy;
    gradient_(1,1) =-dp_dx;

    //Node 2
    xsi(0) = 0.;
    xsi(1) = 1.;
        
    getJacobianMatrix(xsi);
    getSpatialDerivatives(xsi);

    dp_dx = 0.; dp_dy = 0.;
    
    for (int i=0; i<6; i++){
        dp_dx += dphi_dx(0,i) * potential[i];
        dp_dy += dphi_dx(1,i) * potential[i];
    };

    nodes_[connect_(2)] -> setGradientComponent( dp_dy, 0);
    nodes_[connect_(2)] -> setGradientComponent(-dp_dx, 1);

    gradient_(2,0) = dp_dy;
    gradient_(2,1) =-dp_dx;

    //Node 3
    xsi(0) = .5;
    xsi(1) = 0.;
        
    getJacobianMatrix(xsi);
    getSpatialDerivatives(xsi);

    dp_dx = 0.; dp_dy = 0.;
    
    for (int i=0; i<6; i++){
        dp_dx += dphi_dx(0,i) * potential[i];
        dp_dy += dphi_dx(1,i) * potential[i];
    };

    nodes_[connect_(3)] -> setGradientComponent( dp_dy, 0);
    nodes_[connect_(3)] -> setGradientComponent(-dp_dx, 1);

    gradient_(3,0) = dp_dy;
    gradient_(3,1) =-dp_dx;

    //Node 4
    xsi(0) = .5;
    xsi(1) = .5;
        
    getJacobianMatrix(xsi);
    getSpatialDerivatives(xsi);

    dp_dx = 0.; dp_dy = 0.;
    
    for (int i=0; i<6; i++){
        dp_dx += dphi_dx(0,i) * potential[i];
        dp_dy += dphi_dx(1,i) * potential[i];
    };

    nodes_[connect_(4)] -> setGradientComponent( dp_dy, 0);
    nodes_[connect_(4)] -> setGradientComponent(-dp_dx, 1);

    gradient_(4,0) = dp_dy;
    gradient_(4,1) =-dp_dx;

    //Node 5
    xsi(0) = 0.;
    xsi(1) = .5;
        
    getJacobianMatrix(xsi);
    getSpatialDerivatives(xsi);

    dp_dx = 0.; dp_dy = 0.;
    
    for (int i=0; i<6; i++){
        dp_dx += dphi_dx(0,i) * potential[i];
        dp_dy += dphi_dx(1,i) * potential[i];
    };

    nodes_[connect_(5)] -> setGradientComponent( dp_dy, 0);
    nodes_[connect_(5)] -> setGradientComponent(-dp_dx, 1);
    
    gradient_(5,0) = dp_dy;
    gradient_(5,1) =-dp_dx;

    return;
};

//------------------------------------------------------------------------------
//---------COMPUTES VALUES OF POTENTIAL AND GRADIENT IN A SPECIFIC POINT--------
//------------------------------------------------------------------------------
template<>
std::pair<double,ublas::bounded_vector<double,2> > Element<2>::
                            getPotentialResultsInPoint(VecLoc xsi){
    
    shapeQuad.evaluate(xsi,phi_);
    getJacobianMatrix(xsi);
    getSpatialDerivatives(xsi);

    double potential[6];
    
    for (int i=0; i<6; i++){
        potential[i] = nodes_[connect_(i)] -> getPressure();
    };
    
    VecLoc dp_dx; 
    double pot = 0.;
    
    dp_dx.clear();
    
    for (int i=0; i<6; i++){
        dp_dx(0) += dphi_dx(1,i) * potential[i];
        dp_dx(1) -= dphi_dx(0,i) * potential[i];

        pot += potential[i] * phi_(i);
    };

    return std::make_pair(pot,dp_dx);

};

//------------------------------------------------------------------------------
//-----------------------COMPUTES THE VELOCITY DIVERGENT------------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::computeVelocityDivergent(){
    
    ublas::bounded_matrix <double, 2, 6>      int_points;
    
    int_points(0,0) = 0.; int_points(1,0) = 0.;
    int_points(0,1) = 1.; int_points(1,1) = 0.;
    int_points(0,2) = 0.; int_points(1,2) = 1.;
    int_points(0,3) = .5; int_points(1,3) = 0.;
    int_points(0,4) = .5; int_points(1,4) = .5;
    int_points(0,5) = 0.; int_points(1,5) = .5;
    

    typename QuadShapeFunction<2>::Coords xsi;
    
    for (int i=0; i<6; i++){
        
        xsi(0) = int_points(0,i);
        xsi(1) = int_points(1,i);
        
        //Computes the velocity shape functions
        shapeQuad.evaluate(xsi,phi_);
        
        //Computes the jacobian matrix
        getJacobianMatrix(xsi);
        
        //Computes spatial derivatives
        getSpatialDerivatives(xsi);
        
        //Interpolates velocity and its derivatives values
        getVelAndDerivatives();
        
        if(fabs(nodes_[connect_(i)] -> getVelocityDivergent()) <= 
           fabs(du_dx+dv_dy)){
            nodes_[connect_(i)] -> setVelocityDivergent(du_dx + dv_dy);        
        };
    };
    
    return;
};

//------------------------------------------------------------------------------
//------------------------COMPUTES THE VORTICITY FIELD--------------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::computeVorticity(){
    
    ublas::bounded_matrix <double, 6, 6>      least_squares;
    
    typename QuadShapeFunction<2>::Coords xsi;
    ublas::bounded_vector <double, 6>         nodal_values,results;
    
    int in = 0;
    nodal_values.clear();
    least_squares.clear();
    results.clear();

    for(typename NormalQuad::QuadratureListIt it = nQuad.begin(); 
        it != nQuad.end(); it++){
        
        //Defines the integration points adimentional coordinates
        xsi(0) = nQuad.PointList(in,0);
        xsi(1) = nQuad.PointList(in,1);
                
        //Computes the velocity shape functions
        shapeQuad.evaluate(xsi,phi_);
        
        //Computes the jacobian matrix
        getJacobianMatrix(xsi);
        
        //Computes spatial derivatives
        getSpatialDerivatives(xsi);
        
        //Interpolates velocity and its derivatives values
        getVelAndDerivatives();

        for (int i=0; i<6; i++){
            nodal_values(i) += (-du_dy + dv_dx) * phi_(i);
            //nodal_values(i) += du_dxx * phi_(i);
        };
                        
        in++;
    };    

    least_squares.clear();
    
    least_squares(0,0) =  4.825;
    least_squares(0,1) =  1.025;
    least_squares(0,2) =  1.025;
    least_squares(0,3) = -0.275480769230769;  
    least_squares(0,4) =  0.712019230769232;
    least_squares(0,5) = -0.275480769230769;

    least_squares(1,0) =  least_squares(0,1);
    least_squares(1,1) =  4.825;
    least_squares(1,2) =  1.025;
    least_squares(1,3) = -0.275480769230769;  
    least_squares(1,4) = -0.275480769230769;
    least_squares(1,5) =  0.712019230769232;    

    least_squares(2,0) =  least_squares(0,2);
    least_squares(2,1) =  least_squares(1,2);
    least_squares(2,2) =  4.825;
    least_squares(2,3) =  0.712019230769232;  
    least_squares(2,4) = -0.275480769230769;
    least_squares(2,5) = -0.275480769230769;

    least_squares(3,0) =  least_squares(0,3);
    least_squares(3,1) =  least_squares(1,3);
    least_squares(3,2) =  least_squares(2,3);
    least_squares(3,3) =  1.304890902366860;  
    least_squares(3,4) = -0.432609097633136;
    least_squares(3,5) = -0.432609097633136;

    least_squares(4,0) =  least_squares(0,4);
    least_squares(4,1) =  least_squares(1,4);
    least_squares(4,2) =  least_squares(2,4);
    least_squares(4,3) =  least_squares(3,4);  
    least_squares(4,4) =  1.304890902366860;
    least_squares(4,5) = -0.432609097633136;

    least_squares(5,0) =  least_squares(0,5);
    least_squares(5,1) =  least_squares(1,5);
    least_squares(5,2) =  least_squares(2,5);
    least_squares(5,3) =  least_squares(3,5);  
    least_squares(5,4) =  least_squares(4,5);
    least_squares(5,5) =  1.304890902366860;

    results = prod(least_squares,nodal_values);

    for (int i=0; i<6; i++){
        nodes_[connect_(i)] -> incrementVorticity(results(i));
    };
    
    // ublas::bounded_matrix <double, 2, 6>      int_points;
    
    // int_points(0,0) = 0.; int_points(1,0) = 0.;
    // int_points(0,1) = 1.; int_points(1,1) = 0.;
    // int_points(0,2) = 0.; int_points(1,2) = 1.;
    // int_points(0,3) = .5; int_points(1,3) = 0.;
    // int_points(0,4) = .5; int_points(1,4) = .5;
    // int_points(0,5) = 0.; int_points(1,5) = .5;
    

    // typename QuadShapeFunction<2>::Coords xsi;
    
    // for (int i=0; i<6; i++){
        
    //     xsi(0) = int_points(0,i);
    //     xsi(1) = int_points(1,i);
        
    //     //Computes the velocity shape functions
    //     shapeQuad.evaluate(xsi,phi_);
        
    //     //Computes the jacobian matrix
    //     getJacobianMatrix(xsi);
        
    //     //Computes spatial derivatives
    //     getSpatialDerivatives(xsi);
        
    //     //Interpolates velocity and its derivatives values
    //     getVelAndDerivatives();
        
    //     nodes_[connect_(i)] -> incrementVorticity(-du_dy + dv_dx);        
    // };
    
    return;




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
template<>
void Element<2>::getTransientNavierStokes(){

    typename QuadShapeFunction<2>::Coords xsi;
    int index = 0;

    jacobianNRMatrix.clear();
    rhsVector.clear();
    setLocalNodes();
    if (arlequin_){setIntegPointWeightFunction();}


    for(typename NormalQuad::QuadratureListIt it = nQuad.begin(); 
        it != nQuad.end(); it++){

        //Defines the integration points adimentional coordinates
        xsi(0) = nQuad.PointList(index,0);
        xsi(1) = nQuad.PointList(index,1);

        //Computes the velocity shape functions
        shapeQuad.evaluate(xsi,phi_);

        //Returns the quadrature integration weight
        weight_ = nQuad.WeightList(index);

        //Computes the jacobian matrix
        getJacobianMatrix(xsi);

        //Computes spatial derivatives
        getSpatialDerivatives(xsi);

        //Interpolates velocity and its derivatives values
        getVelAndDerivatives();

        //Compute Stabilization Parameters
        getParameterSUPG();

        //Computes the element diffusion/viscosity matrix
        getElemMatrix(index);

        //Computes the RHS vector
        getResidualVector(index); 

        index++;        
    };  
    
    //getBoundaryIntegration();

    //Apply boundary conditions
    setBoundaryConditions();

  

    return;
};

//------------------------------------------------------------------------------
//----------------------------STEADY LAPLACE PROBEM-----------------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::getSteadyLaplace(){

    typename QuadShapeFunction<2>::Coords xsi;
    int index = 0;

    laplMatrix.clear();
    jacobianNRMatrix.clear();

    setLocalNodes();   


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

    
    for(typename NormalQuad::QuadratureListIt it = nQuad.begin(); 
        it != nQuad.end(); it++){
        
        //Defines the integration points adimentional coordinates
        xsi(0) = nQuad.PointList(index,0);
        xsi(1) = nQuad.PointList(index,1);

        //Computes the velocity shape functions
        shapeQuad.evaluate(xsi,phi_);

        //Returns the quadrature integration weight
        weight_ = nQuad.WeightList(index);

        //Computes the jacobian matrix
        getJacobianMatrix(xsi);

        //Computes spatial derivatives
        getSpatialDerivatives(xsi);




        // ublas::bounded_matrix<double, 3,18> bMatrix,aux;
        // bMatrix.clear(); aux.clear();

        // for (int i = 0; i < 6; i++){                
        //     bMatrix(0,2*i  ) = dphi_dx(0,i);
        //     bMatrix(2,2*i  ) = dphi_dx(1,i);

        //     bMatrix(1,2*i+1) = dphi_dx(1,i);
        //     bMatrix(2,2*i+1) = dphi_dx(0,i);
        // };
        

        // aux = prod(hooke,bMatrix);        

        // laplMatrix += prod(trans(bMatrix),aux) * djac_ * weight_;




        getElemLaplMatrix();

        index++;        
    };  
    
    //Computes the Newton's method jacobian
    jacobianNRMatrix = laplMatrix;
    
    //Computes the RHS vector
    getResidualVectorLaplace();

    //Apply boundary conditions
    setBoundaryConditionsLaplace();

    return;
};

//------------------------------------------------------------------------------
//----------------------------STEADY LAPLACE PROBEM-----------------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::getLagrangeMultipliersSameMesh(){

    typename QuadShapeFunction<2>::Coords xsi;
    int index = 0;

    double una_ = timeScheme_ * u_ + (1. - timeScheme_) * uPrev_;
    double vna_ = timeScheme_ * v_ + (1. - timeScheme_) * vPrev_;

    diffMatrix.clear();
    jacobianNRMatrix.clear();
    rhsVector.clear();
    lagrMultMatrix.clear();
    laplMatrix.clear();
    arlequinStab.clear();
    arlequinStabVector.clear();
    
    for(typename NormalQuad::QuadratureListIt it = nQuad.begin(); 
        it != nQuad.end(); it++){
        
        //Defines the integration points adimentional coordinates
        xsi(0) = nQuad.PointList(index,0);
        xsi(1) = nQuad.PointList(index,1);

        //Computes the velocity shape functions
        shapeQuad.evaluate(xsi,phi_);

        //Returns the quadrature integration weight
        weight_ = nQuad.WeightList(index);

        //Computes the jacobian matrix
        getJacobianMatrix(xsi);

        getSpatialDerivatives(xsi);
        
        getVelAndDerivatives();     

        // H1 COUPLING OPERATOR        
        ublas::bounded_matrix<double, 3,18> Bmatrix, Baux;
        ublas::bounded_matrix<double, 18,3> BmatrixT; 
        ublas::bounded_matrix<double, 3,3>  Maux;

        Bmatrix.clear(); BmatrixT.clear(); Maux.clear(); Baux.clear();

        for (int i = 0; i < 6; i++){
            Bmatrix(0,2*i  ) = dphi_dx(0,i);
            Bmatrix(2,2*i  ) = dphi_dx(1,i);       
            
            BmatrixT(2*i  ,0) = dphi_dx(0,i);
            BmatrixT(2*i  ,2) = dphi_dx(1,i);        
            
            Bmatrix(1,2*i+1) = dphi_dx(1,i);
            Bmatrix(2,2*i+1) = dphi_dx(0,i);
            
            BmatrixT(2*i+1,1) = dphi_dx(1,i);
            BmatrixT(2*i+1,2) = dphi_dx(0,i);
        };
        
        Maux(0,0) =  2.;
        Maux(1,1) =  2.;
        Maux(2,2) =  1.;
        
        Baux = prod(Maux,Bmatrix);
        
        diffMatrix += prod(BmatrixT,Baux) * weight_ * djac_ * k2;
        
        for (int i = 0; i < 6; i++){
            for (int j = 0; j < 6; j++){        

                // L2 COUPLING OPERATOR
                lagrMultMatrix(2*i  ,2*j  ) += 
                    (phi_(i) * phi_(j))
                    * weight_ * djac_ * k1;
                lagrMultMatrix(2*i+1,2*j+1) += 
                    (phi_(i) * phi_(j))
                    * weight_ * djac_ * k1;

                //SUPG STABILIZATION TERM
                double LM = 0.;
                // LM = - ((una_ - umesh_) * dphi_dx(0,i) + (vna_ - vmesh_) * dphi_dx(1,i))
                //     * phi_(j) * tSUPG_;
            
                jacobianNRMatrix(2*i  ,2*j  ) += (LM)
                    * weight_ * djac_;
                jacobianNRMatrix(2*i+1,2*j+1) += (LM)
                    * weight_ * djac_;

                //PSPG STABILIZATION TERM
                double Lx = 0.;
                double Ly = 0.;
                double Tx = 0.;
                double Ty = 0.;

                // Lx = dphi_dx(0,i) * phi_(j) * tPSPG_ / dens_ * k1;
                // Ly = dphi_dx(1,i) * phi_(j) * tPSPG_ / dens_ * k1;

                // Tx = (dphi_dx(0,j) * (2. * ddphi_dx(0,0)(i) + 
                //       ddphi_dx(1,1)(i)) + dphi_dx(1,j) * ddphi_dx(0,1)(i))
                //       * tPSPG_ / dens_ * k2;
                // Ty = (dphi_dx(1,j) * (2. * ddphi_dx(1,1)(i) + 
                //                       ddphi_dx(0,0)(i)) + dphi_dx(0,j) * ddphi_dx(0,1)(i))
                //         * tPSPG_ / dens_ * k2;

                jacobianNRMatrix(2*i  ,12+j) += (Lx + Tx)// * intPointWeightFunction(index)
                    * weight_ * djac_;
                jacobianNRMatrix(2*i+1,12+j) += (Ly + Ty)// * intPointWeightFunction(index)
                    * weight_ * djac_;


                //ARLEQUIN STABILIZATION TERMS
                double AM = 0.;
                double Lpx = 0.; double Lpy = 0.;
                double LC = 0.; double LL = 0.;

                // AM = phi_(i) * phi_(j) * tARLQ_ * intPointWeightFunction(index);

                // LL = -2 * phi_(i) * phi_(j) * tARLQ_ / dens_;
                LL = (dphi_dx(0,i) * dphi_dx(0,j) + dphi_dx(1,i) * dphi_dx(1,j)) * tARLQ_ / dens_;

                // Lpx = -(dphi_dx(0,i) * dp_dxx + dphi_dx(1,i) * dp_dxy) * intPointWeightFunction(index)
                //      * tARLQ_ / dens_;
                // Lpy = -(dphi_dx(0,i) * dp_dxy + dphi_dx(1,i) * dp_dyy) * intPointWeightFunction(index)
                //      * tARLQ_ / dens_;

                // LC = phi_(i) * ((una_ - umesh_) * dphi_dx(0,i) + (vna_ - vmesh_) * dphi_dx(1,i)) * phi_(j) * tARLQ_ * intPointWeightFunction(index);

                arlequinStab(2*i  ,2*j  ) += (LL + Lpx) 
                    * weight_ * djac_;
                arlequinStab(2*i+1,2*j+1) += (LL + Lpy) 
                    * weight_ * djac_;
                arlequinStab(12+i,12+j) += 0 * weight_ * djac_;

                // LC = -(dphi_dx(0,i)*(du_dx*dphi_dx(0,j) + dv_dx*dphi_dx(1,j)) +
                //        dphi_dx(1,i)*(du_dy*dphi_dx(0,j) + dv_dy*dphi_dx(1,j)) + 
                //        dphi_dx(0,i)*(u_*ddphi_dx(0,0)(j) + v_*ddphi_dx(0,1)(j)) + 
                //        dphi_dx(1,i)*(u_*ddphi_dx(1,0)(j) + v_*ddphi_dx(1,1)(j))) * tARLQ_ * intPointWeightFunction(index);

                // Lpx = 0.; Lpy = 0.;

                // Lpx = (dphi_dx(0,i) * ddphi_dx(0,0)(j) + 
                //        dphi_dx(1,i) * ddphi_dx(0,1)(j)) * tARLQ_ * intPointWeightFunction(index);
                // Lpy = (dphi_dx(0,i) * ddphi_dx(1,0)(j) + 
                //        dphi_dx(1,i) * ddphi_dx(1,1)(j)) * tARLQ_ * intPointWeightFunction(index);




                laplMatrix(2*i  ,2*j  ) += LC * weight_ * djac_;
                laplMatrix(2*i+1,2*j+1) += LC * weight_ * djac_;

                laplMatrix(2*i  ,12+j) += Lpx * weight_ * djac_;
                laplMatrix(2*i+1,12+j) += Lpy * weight_ * djac_;
                laplMatrix(12+j,2*i  ) += Lpx * weight_ * djac_;
                laplMatrix(12+j,2*i+1) += Lpy * weight_ * djac_;

            };

            //SUPG AND PSPG STABILIZATION TERMS
            double LMx = 0.;
            double LMy = 0.;
            double LMp = 0.;
            
            LMx = (((una_ - umesh_) * dphi_dx(0,i) + (vna_ - vmesh_) * dphi_dx(1,i))
                   * lagMx_) * tSUPG_;
            LMy = (((una_ - umesh_) * dphi_dx(0,i) + (vna_ - vmesh_) * dphi_dx(1,i))
                   * lagMy_) * tSUPG_;
            // LMp = (dphi_dx(0,i) * lagMx_ + dphi_dx(1,i) * lagMy_)
            //         * tPSPG_ / dens_ * k1 + 
            //       (dphi_dx(0,i) * (2. * dLx_dxx + dLx_dxy + dLx_dyy) +
            //        dphi_dx(1,i) * (2. * dLy_dyy + dLy_dxy + dLy_dxx)) * tPSPG_ / dens_ * k2;

            rhsVector(2*i  ) += (LMx) * weight_ * djac_;
            rhsVector(2*i+1) += (LMy) * weight_ * djac_;
            rhsVector(12+i) += (LMp) * weight_ * djac_;

            //ARLEQUIN STABILIZATION TERMS
            double Amx = 0.; double Amy = 0.;
            double LCx = 0.; double LCy = 0.;
            double LPx = 0.; double LPy = 0.;
            double LLx = 0.; double LLy = 0.;

            // Amx = -phi_(i) * (una_ - uPrev_) / dTime_ *
            //         intPointWeightFunction(index) * tARLQ_;
            // Amy = -phi_(i) * (vna_ - vPrev_) / dTime_ *
            //         intPointWeightFunction(index) * tARLQ_;


            //if (iTimeStep < 10){
                LLx = -(dphi_dx(0,i) * (dLx_dx/intPointWeightFunction(index)) + dphi_dx(1,i) * (dLx_dy/intPointWeightFunction(index))) * tARLQ_ / dens_;
                LLy = -(dphi_dx(0,i) * (dLy_dx/intPointWeightFunction(index)) + dphi_dx(1,i) * (dLy_dy/intPointWeightFunction(index))) * tARLQ_ / dens_;
           //}else{
                // LLx = -(dphi_dx(0,i) * (dLx_dx/intPointWeightFunction(index) - du_dx*du_dx - dv_dx*du_dy) + 
                //         dphi_dx(1,i) * (dLx_dy/intPointWeightFunction(index) - du_dy*du_dx - dv_dy*du_dy)) * tARLQ_ / dens_;
                // LLy = -(dphi_dx(0,i) * (dLy_dx/intPointWeightFunction(index) - du_dx*dv_dx - dv_dx*dv_dy) + 
                //         dphi_dx(1,i) * (dLy_dy/intPointWeightFunction(index) - du_dy*dv_dx - dv_dy*dv_dy)) * tARLQ_ / dens_;  
                // LLx = -(dphi_dx(0,i) * (dLx_dx/intPointWeightFunction(index) - dp_dxx) + 
                //         dphi_dx(1,i) * (dLx_dy/intPointWeightFunction(index) - dp_dxy)) * tARLQ_ / dens_;
                // LLy = -(dphi_dx(0,i) * (dLy_dx/intPointWeightFunction(index) - dp_dxy) + 
                //         dphi_dx(1,i) * (dLy_dy/intPointWeightFunction(index) - dp_dyy)) * tARLQ_ / dens_;  

                // LLx += - (dphi_dx(0,i)*dax_dx + dphi_dx(1,i)*dax_dy) * tARLQ_ / dens_;
                // LLx += - (dphi_dx(0,i)*day_dx + dphi_dx(1,i)*day_dy) * tARLQ_ / dens_;

            //}
            
            

            arlequinStabVector(2*i  ) += (LLx) * weight_ * djac_ * intPointWeightFunction(index);
            arlequinStabVector(2*i+1) += (LLy) * weight_ * djac_ * intPointWeightFunction(index);





        };         
        
        index++;        
    }; 

    // jacobianNRMatrix.clear();
    // rhsVector.clear();
    // jacobianNRMatrix += diffMatrix;


    lagrMultMatrix += diffMatrix;

    //lagrMultMatrix.clear();


    // ublas::identity_matrix<double>             ident(18);
    // lagrMultMatrix = ident;

    // lagrMultMatrix(12,12) = 0.;
    // lagrMultMatrix(13,13) = 0.;
    // lagrMultMatrix(14,14) = 0.;

    // //Computes the Newton's method jacobian
    //jacobianNRMatrix = lagrMultMatrix;

    // Set boundary conditions to the Lagrange Multipliers
    //setBoundaryConditionsLagrangeMultipliers();

    // rhsVectorLM.clear();  
    // LocalVector U_;
    // U_.clear();

    // for (int i=0; i<6; i++){
    //     if (nodes_[connect_(i)] -> getConstrains(0) == 0){
    //         U_(2*i  ) = nodes_[connect_(i)] -> getLagrangeMultiplier(0);
    //     }
    //     if (nodes_[connect_(i)] -> getConstrains(1) == 0){
    //         U_(2*i+1) = nodes_[connect_(i)] -> getLagrangeMultiplier(1);
    //     }
    //     // if(index_ == 32) std::cout << "Vel " << i << " " << U_(2*i+1) 
    //     //         << " " << nodes_[connect_(i)] -> getConstrains(1) <<  std::endl;
    // };

    // noalias(rhsVectorLM) = - prod(lagrMultMatrix,U_); 

    // for (int i = 0; i < 6; i++){
    //     if (nodes_[connect_(i)] -> getConstrains(0) == 1) {
    //         for (int j = 0; j < 18; j++){
    //             //lagrMultMatrix(2*i  ,j) = 0.;
    //             lagrMultMatrix(j,2*i  ) = 0.;
    //         };
    //         //lagrMultMatrix(12+i,12+i) = 1.;
    //         rhsVectorLM(2*i  ) = 0.0;
    //     };

    //     if (nodes_[connect_(i)] -> getConstrains(1) == 1) {
    //         // std::cout << "aqui" << std::endl;
    //         for (int j = 0; j < 18; j++){
    //             //lagrMultMatrix(2*i+1,j) = 0.;
    //             lagrMultMatrix(j,2*i+1) = 0.;
    //         };
    //         //lagrMultMatrix(12+i,12+i) = 1.;
    //         rhsVectorLM(2*i+1) = 0.0;
    //     };
    // };



    for (int i = 0; i < 6; i++){
        if (nodes_[connect_(i)] -> getConstrains(0) == 1) {
            for (int j = 0; j < 18; j++){
                arlequinStab(2*i  ,j) = 0.;
                arlequinStab(j,2*i  ) = 0.;
            };
            //lagrMultMatrix(12+i,12+i) = 1.;
            arlequinStabVector(2*i  ) = 0.0;
        };

        if (nodes_[connect_(i)] -> getConstrains(1) == 1) {
            // std::cout << "aqui" << std::endl;
            for (int j = 0; j < 18; j++){
                arlequinStab(2*i+1,j) = 0.;
                arlequinStab(j,2*i+1) = 0.;
            };
            //lagrMultMatrix(12+i,12+i) = 1.;
            arlequinStabVector(2*i+1) = 0.0;
        };
    };




    //!!!!! O erro é que não pode somar aqui porque senão vai adicionar isso na equação dos multiplicadores também e só deve adicionar na parte da equação do momentum
    // lagrMultMatrix += jacobianNRMatrix;
    // rhsVectorLM += rhsVectorLM;


    return;
};

//------------------------------------------------------------------------------
//----------------------------STEADY LAPLACE PROBEM-----------------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::getLagrangeMultipliersDifferentMesh(int ielem,double tPSPG2_,ublas::bounded_vector<double, 6> press, ublas::bounded_vector<double, 6> velx, ublas::bounded_vector<double, 6> vely){

    typename QuadShapeFunction<2>::Coords xsi,xsi_intp;
    int index = 0;
    //tARLQ_ = -tPSPG2_;

    typename QuadShapeFunction<2>::Values phiLM_; 
    
    double una_ = timeScheme_ * u_ + (1. - timeScheme_) * uPrev_;
    double vna_ = timeScheme_ * v_ + (1. - timeScheme_) * vPrev_;
    
    diffMatrix.clear();
    jacobianNRMatrix.clear();
    rhsVector.clear();
    lagrMultMatrix.clear();
    arlequinStab.clear();
    arlequinStabVector.clear();
    laplMatrix.clear();
       //std::cout << "PSPG Fine " << ielem << " " << tPSPG_ << std::endl;
    for(typename SpecialQuad::QuadratureListIt it = sQuad.begin(); 
        it != sQuad.end(); it++){
        
        if ((intPointCorrespElem(index) == ielem)){

            //Defines the integration points adimentional coordinates
            xsi(0) = sQuad.PointList(index,0);
            xsi(1) = sQuad.PointList(index,1);
            
            //Computes the velocity shape functions
            shapeQuad.evaluate(xsi,phi_);
            
            xsi_intp(0) = intPointCorrespXsi(index,0);
            xsi_intp(1) = intPointCorrespXsi(index,1);

            //Computes the coarse mesh shape functions
            shapeQuad.evaluate(xsi_intp,phiLM_);
            
            //Returns the quadrature integration weight
            weight_ = sQuad.WeightList(index);
            
            //Computes the jacobian matrix
            getJacobianMatrix(xsi_intp);
                        
            getSpatialDerivatives(xsi_intp);

            dphiL_dx = dphi_dx;
            ddphiL_dx = ddphi_dx;

            getJacobianMatrix(xsi);
            getSpatialDerivatives(xsi);

            u_ = 0.;
            v_ = 0.;

            du_dx = 0.;
            du_dy = 0.;
            dv_dx = 0.;
            dv_dy = 0.;

            du_dxx = 0.;
            du_dyy = 0.;
            du_dxy = 0.;
            dv_dxx = 0.;
            dv_dyy = 0.;
            dv_dxy = 0.;

            dp_dxx = 0.;
            dp_dxy = 0.;
            dp_dyx = 0.;
            dp_dyy = 0.;

    //Interpolates the velocity components and its spatial derivatives
            for (int i = 0; i < 6; i++){
                u_ += velx(i) * phiLM_(i);
                v_ += vely(i) * phiLM_(i);

                du_dx += velx(i) * dphiL_dx(0,i);
                du_dy += velx(i) * dphiL_dx(1,i);
                dv_dx += vely(i) * dphiL_dx(0,i);
                dv_dy += vely(i) * dphiL_dx(1,i);

                du_dxx += velx(i) * ddphiL_dx(0,0)(i);
                du_dyy += velx(i) * ddphiL_dx(1,1)(i);
                du_dxy += velx(i) * ddphiL_dx(0,1)(i);
                dv_dxx += vely(i) * ddphiL_dx(0,0)(i);
                dv_dyy += vely(i) * ddphiL_dx(1,1)(i);
                dv_dxy += vely(i) * ddphiL_dx(0,1)(i);

                dp_dxx += press(i) * ddphiL_dx(0,0)(i);
                dp_dxy += press(i) * ddphiL_dx(0,1)(i);
                dp_dyx += press(i) * ddphiL_dx(1,0)(i);
                dp_dyy += press(i) * ddphiL_dx(1,1)(i);
            };  

            ublas::bounded_matrix<double, 3,18> Bmatrix, Baux;
            ublas::bounded_matrix<double, 18,3> BmatrixT; 
            ublas::bounded_matrix<double, 3,3>  Maux;
            
            //H1 COUPLING OPERATOR
            Bmatrix.clear(); BmatrixT.clear(); Maux.clear(); Baux.clear();
            
            for (int i = 0; i < 6; i++){
                Bmatrix(0,2*i  ) = dphi_dx(0,i);
                Bmatrix(2,2*i  ) = dphi_dx(1,i);       
                
                BmatrixT(2*i  ,0) = dphiL_dx(0,i);
                BmatrixT(2*i  ,2) = dphiL_dx(1,i);        
                
                Bmatrix(1,2*i+1) = dphi_dx(1,i);
                Bmatrix(2,2*i+1) = dphi_dx(0,i);
                
                BmatrixT(2*i+1,1) = dphiL_dx(1,i);
                BmatrixT(2*i+1,2) = dphiL_dx(0,i);
            };
            
            Maux(0,0) =  2.;
            Maux(1,1) =  2.;
            Maux(2,2) =  1.;
            
            Baux = prod(Maux,Bmatrix);
            
            diffMatrix += prod(BmatrixT,Baux) * weight_ * djac_ * k2;


            p_ = 0.;
            dp_dxx = 0.;
            dp_dxy = 0.;
            dp_dyx = 0.;
            dp_dyy = 0.;

            // for (int i = 0; i < 6; i++){
            //     p_ += phiLM_(i) * press(i);

            //     dp_dxx += press(i) * ddphiL_dx(0,0)(i);
            //     dp_dxy += press(i) * ddphiL_dx(0,1)(i);
            //     dp_dyx += press(i) * ddphiL_dx(1,0)(i);
            //     dp_dyy += press(i) * ddphiL_dx(1,1)(i);
            // }


            
            for (int i = 0; i < 6; i++){
                for (int j = 0; j < 6; j++){     
                    //L2 COUPLING OPERATOR   
                    lagrMultMatrix(2*i  ,2*j  ) += (phiLM_(i) * phi_(j)) * weight_ * djac_ * k1;
                    lagrMultMatrix(2*i+1,2*j+1) += (phiLM_(i) * phi_(j)) * weight_ * djac_ * k1;

                    //SUPG STABILIZATION TERM
                    double LM = 0.;
                    // LM = ((una_ - umesh_) * dphi_dx(0,j) + (vna_ - vmesh_) * dphi_dx(0,j))
                    //     * phiLM_(i) * tPSPG2_;
                    
                    jacobianNRMatrix(2*i  ,2*j  ) += (timeScheme_ * dTime_ * LM)
                        * weight_ * djac_;
                    jacobianNRMatrix(2*i+1,2*j+1) += (timeScheme_ * dTime_ * LM)
                        * weight_ * djac_;
                    
                    //PSPG STABILIZATION TERM
                    double Lx = 0.; double Ly = 0.;
                    
                    // Lx = -dphiL_dx(0,i) * phi_(j) * tPSPG_ / dens_ * k1;
                    // Ly = -dphiL_dx(1,i) * phi_(j) * tPSPG_ / dens_ * k1;
                    
                    jacobianNRMatrix(2*i  ,12+j) += (Lx) * timeScheme_ * dTime_
                        * weight_ * djac_;
                    jacobianNRMatrix(2*i+1,12+j) += (Ly) * timeScheme_ * dTime_
                        * weight_ * djac_;  


                    double AM = 0.;
                    double Lpx = 0.; double Lpy = 0.;
                    double LC = 0.; double LL = 0.;

                    // AM = -phiLM_(i) * phi_(j) * 
                    // intPointWeightFunction(index) * tARLQ_;

                    // LL = phiLM_(i) * phi_(j) * tARLQ_ / dens_;
                    LL = (dphi_dx(0,i) * dphi_dx(0,j) + dphi_dx(1,i) * dphi_dx(1,j)) * tARLQ_ / dens_;

                    // Lpx = phi_(i) * dphi_dx(0,i) * intPointWeightFunction(index)
                    //     * tARLQ_ / dens_;
                    // Lpy = phi_(i) * dphi_dx(1,i) * intPointWeightFunction(index)
                    //     * tARLQ_ / dens_;

                    // LC = -phi_(j) * ((una_ - umesh_) * dphi_dx(0,j) + (vna_ - vmesh_) * dphi_dx(1,j)) * phiLM_(i) * tARLQ_ * intPointWeightFunction(index);

                    arlequinStab(2*i  ,2*j  ) += (LL) * weight_ * djac_;
                    arlequinStab(2*i+1,2*j+1) += (LL) * weight_ * djac_;
                    arlequinStab(12+i,12+j) += 0 * weight_ * djac_;

                    // LC = (dphi_dx(0,i)*(du_dx*dphiL_dx(0,j) + dv_dx*dphiL_dx(1,j)) +
                    //       dphi_dx(1,i)*(du_dy*dphiL_dx(0,j) + dv_dy*dphiL_dx(1,j)) + 
                    //       dphi_dx(0,i)*(u_*ddphiL_dx(0,0)(j) + v_*ddphiL_dx(0,1)(j)) + 
                    //       dphi_dx(1,i)*(u_*ddphiL_dx(1,0)(j) + v_*ddphiL_dx(1,1)(j))) * tARLQ_ * (1-intPointWeightFunction(index));

                    // laplMatrix(2*i  ,2*j  ) += LC * weight_ * djac_;
                    // laplMatrix(2*i+1,2*j+1) += LC * weight_ * djac_;

                     Lpx = 0.; Lpy = 0.;

                    // Lpx = (dphi_dx(0,i) * ddphi_dx(0,0)(j) + 
                    //        dphi_dx(1,i) * ddphi_dx(0,1)(j)) * tARLQ_ * intPointWeightFunction(index);
                    // Lpy = (dphi_dx(0,i) * ddphi_dx(1,0)(j) + 
                    //        dphi_dx(1,i) * ddphi_dx(1,1)(j)) * tARLQ_ * intPointWeightFunction(index);

                    laplMatrix(2*i  ,2*j  ) += LC * weight_ * djac_;
                    laplMatrix(2*i+1,2*j+1) += LC * weight_ * djac_;

                    laplMatrix(2*i  ,12+j) += Lpx * weight_ * djac_;
                    laplMatrix(2*i+1,12+j) += Lpy * weight_ * djac_;
                    laplMatrix(12+j,2*i  ) += Lpx * weight_ * djac_;
                    laplMatrix(12+j,2*i+1) += Lpy * weight_ * djac_;

                };

                //SUPG AND PSPG STABILIZATION TERMS
                double LMx = 0.; double LMy = 0.; double LMp = 0.;
             
                LMx = ((una_ - umesh_) * dphiL_dx(0,i) + 
                       (vna_ - vmesh_) * dphiL_dx(1,i)) * lagMx_ * -tPSPG2_;
                LMy = ((una_ - umesh_) * dphiL_dx(0,i) + 
                       (vna_ - vmesh_) * dphiL_dx(1,i)) * lagMy_ * -tPSPG2_;
                // // LMp = (dphiL_dx(0,i) * lagMx_ + dphiL_dx(1,i) * lagMy_)
                //     * tPSPG2_ / dens_ * k1 +
                //       (dphiL_dx(0,i) * (2. * dLx_dxx + dLx_dxy + dLx_dyy) +
                //        dphiL_dx(1,i) * (2. * dLy_dyy + dLy_dxy + dLy_dxx)) * tPSPG2_ / dens_ * k2;

                rhsVector(2*i  ) += (LMx) * weight_ * djac_;
                rhsVector(2*i+1) += (LMy) * weight_ * djac_;
                rhsVector(12+i) += (LMp) * weight_ * djac_;

                //ARLEQUIN STABILIZATION TERMS
                double Amx = 0.; double Amy = 0.;
                double LCx = 0.; double LCy = 0.;
                double LPx = 0.; double LPy = 0.;
                double LLx = 0.; double LLy = 0.;
                
               
                //if (iTimeStep < 10){
                    LLx = -(dphi_dx(0,i) * dLx_dx/(1-intPointWeightFunctionSpecial(index)) + dphi_dx(1,i) * dLx_dy/(1-intPointWeightFunctionSpecial(index))) * tARLQ_ / dens_;
                    LLy = -(dphi_dx(0,i) * dLy_dx/(1-intPointWeightFunctionSpecial(index)) + dphi_dx(1,i) * dLy_dy/(1-intPointWeightFunctionSpecial(index))) * tARLQ_ / dens_;
                //}else{
                    // LLx = -(dphi_dx(0,i) * (dLx_dx/(1-intPointWeightFunction(index)) + du_dx*du_dx + dv_dx*du_dy) + 
                    //         dphi_dx(1,i) * (dLx_dy/(1-intPointWeightFunction(index)) + du_dy*du_dx + dv_dy*du_dy)) * tARLQ_ / dens_;
                    // LLy = -(dphi_dx(0,i) * (dLy_dx/(1-intPointWeightFunction(index)) + du_dx*dv_dx + dv_dx*dv_dy) + 
                    //         dphi_dx(1,i) * (dLy_dy/(1-intPointWeightFunction(index)) + du_dy*dv_dx + dv_dy*dv_dy)) * tARLQ_ / dens_;
                    // LLx = -(dphi_dx(0,i) * (dLx_dx/(1-intPointWeightFunction(index)) + dp_dxx) + 
                    //         dphi_dx(1,i) * (dLx_dy/(1-intPointWeightFunction(index)) + dp_dxy)) * tARLQ_ / dens_;
                    // LLy = -(dphi_dx(0,i) * (dLy_dx/(1-intPointWeightFunction(index)) + dp_dxy) + 
                    //         dphi_dx(1,i) * (dLy_dy/(1-intPointWeightFunction(index)) + dp_dyy)) * tARLQ_ / dens_;
                //};
     

                arlequinStabVector(2*i  ) += (LLx) * weight_ * djac_ * (1-intPointWeightFunctionSpecial(index)) ;
                arlequinStabVector(2*i+1) += (LLy) * weight_ * djac_ * (1-intPointWeightFunctionSpecial(index)) ;
            };
        };
        index++;        
    };  
    // jacobianNRMatrix.clear();
    // rhsVector.clear();
    

    lagrMultMatrix += diffMatrix;

    //lagrMultMatrix.clear();
    //std::cout << "AQUI " << std::endl;


    // rhsVectorLM.clear();  
    // LocalVector U_;
    // U_.clear();

    // for (int i=0; i<6; i++){
    //     if (nodes_[connect_(i)] -> getConstrains(0) == 0){
    //         U_(2*i  ) = nodes_[connect_(i)] -> getLagrangeMultiplier(0);
    //     }
    //     if (nodes_[connect_(i)] -> getConstrains(1) == 0){
    //         U_(2*i+1) = nodes_[connect_(i)] -> getLagrangeMultiplier(1);
    //     }
    //     // if(index_ == 32) std::cout << "Vel " << i << " " << U_(2*i+1) 
    //     //         << " " << nodes_[connect_(i)] -> getConstrains(1) <<  std::endl;
    // };

    // noalias(rhsVectorLM) = -prod(lagrMultMatrix,U_); 

    // for (int i = 0; i < 6; i++){
    //     if (nodes_[connect_(i)] -> getConstrains(0) == 1) {
    //         for (int j = 0; j < 18; j++){
    //             //lagrMultMatrix(2*i  ,j) = 0.;
    //             lagrMultMatrix(j,2*i  ) = 0.;
    //         };
    //         //lagrMultMatrix(12+i,12+i) = 1.;
    //         rhsVectorLM(2*i  ) = 0.0;
    //     };

    //     if (nodes_[connect_(i)] -> getConstrains(1) == 1) {
    //         // std::cout << "aqui" << std::endl;
    //         for (int j = 0; j < 18; j++){
    //             //lagrMultMatrix(2*i+1,j) = 0.;
    //             lagrMultMatrix(j,2*i+1) = 0.;
    //         };
    //         //lagrMultMatrix(12+i,12+i) = 1.;
    //         rhsVectorLM(2*i+1) = 0.0;
    //     };
    // };

    for (int i = 0; i < 6; i++){
        if (nodes_[connect_(i)] -> getConstrains(0) == 1) {
            for (int j = 0; j < 18; j++){
                arlequinStab(2*i  ,j) = 0.;
                arlequinStab(j,2*i  ) = 0.;
            };
            //lagrMultMatrix(12+i,12+i) = 1.;
            arlequinStabVector(2*i  ) = 0.0;
            std::cout << "AQUI2 " << std::endl;
        };

        if (nodes_[connect_(i)] -> getConstrains(1) == 1) {
            // std::cout << "aqui" << std::endl;
            for (int j = 0; j < 18; j++){
                arlequinStab(2*i+1,j) = 0.;
                arlequinStab(j,2*i+1) = 0.;
            };
            //lagrMultMatrix(12+i,12+i) = 1.;
            arlequinStabVector(2*i+1) = 0.0;
            std::cout << "AQUI3 " << std::endl;
        };
    };



    return;
};


//------------------------------------------------------------------------------
//----------------------------STEADY LAPLACE PROBEM-----------------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::getLMStabilizationSameMesh(){

    //  typename QuadShapeFunction<2>::Coords xsi;
    // int index = 0;

    // diffMatrix.clear();
    // jacobianNRMatrix.clear();
    // rhsVector.clear();
    // lagrMultMatrix.clear();
    
    // for(typename NormalQuad::QuadratureListIt it = nQuad.begin(); 
    //     it != nQuad.end(); it++){
        
    //     //Defines the integration points adimentional coordinates
    //     xsi(0) = nQuad.PointList(index,0);
    //     xsi(1) = nQuad.PointList(index,1);

    //     //Computes the velocity shape functions
    //     shapeQuad.evaluate(xsi,phi_);

    //     //Returns the quadrature integration weight
    //     weight_ = nQuad.WeightList(index);

    //     //Computes the jacobian matrix
    //     getJacobianMatrix(xsi);

    //     getSpatialDerivatives(xsi);
        
    //     getVelAndDerivatives();     
        
        
    //     for (int i = 0; i < 6; i++){
    //         for (int j = 0; j < 6; j++){        
                
    //         //Mass matrix (use for both directions)
    //         double mM =  phi_(i) * phi_(j) * dens_ * tARLQ_;
                         
    //         //Convection matrix
    //         double Cxx = phi_(i) * ((u_ - umesh_) * dphi_dx(0,j) + 
    //                                 (v_ - vmesh_) * dphi_dx(1,j)) 
    //             * dens_ * tARLQ_;
            
    //         double Cyy  = Cxx;

            
    //         double Cuu = phi_(i) * (du_dx + dv_dx) * phi_(j) * dens_ * tARLQ_;
                        
    //         double Cvv = phi_(i) * (du_dy + dv_dy) * phi_(j) * dens_ * tARLQ_;
            
    //         // jacobianNRMatrix(2*i  ,2*j  ) += (mM + timeScheme_ * dTime_ * 
    //         //                                   (Cxx + Cuu))
    //         //     * weight_ * djac_;

    //         // jacobianNRMatrix(2*i+1,2*j+1) += (mM + timeScheme_ * dTime_ *
    //         //                                   (Cyy + Cvv))
    //         //     * weight_ * djac_;

    //         };

    //         double una_ = timeScheme_ * u_ + (1. - timeScheme_) * uPrev_;
    //         double vna_ = timeScheme_ * v_ + (1. - timeScheme_) * vPrev_;
            
    //         double mx = phi_(i) * (u_ - uPrev_) * tARLQ_ * dens_;
    //         double my = phi_(i) * (v_ - vPrev_) * tARLQ_ * dens_;
        
    //         double Cx = phi_(i) * ((una_ - umesh_) * du_dx + 
    //                                (vna_ - vmesh_) * du_dy) * tARLQ_ *dens_ + 
    //             phi_(i) * (du_dx + dv_dx) * u_ * dens_ * tARLQ_;
    //         double Cy = phi_(i) * ((una_ - umesh_) * dv_dx + 
    //                                (vna_ - vmesh_) * dv_dy) * tARLQ_ *dens_ + 
    //             phi_(i) * (du_dy + dv_dy) * v_ * dens_ * tARLQ_;

            
    //         // rhsVector(2*i  ) += (-mx - Cx * dTime_ * timeScheme_) * 
    //         //     weight_ * djac_;
    //         // rhsVector(2*i+1) += (-my - Cy * dTime_ * timeScheme_) * 
    //         //     weight_ * djac_;
    //     };         
        
    //     index++;        
    // }; 

    return;
};


//------------------------------------------------------------------------------
//----------------------ELEMENT DIFFUSION/VISCOSITY MATRIX----------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::getStabCoarse(int ielem){

    // typename QuadShapeFunction<2>::Coords xsi, xsi_intp;
    // int index = 0;

    // typename QuadShapeFunction<2>::Values phiLM_; 
   
    // diffMatrix.clear();
    // jacobianNRMatrix.clear();
    // rhsVector.clear();
    // lagrMultMatrix.clear();
    
    // for(typename SpecialQuad::QuadratureListIt it = sQuad.begin(); 
    //     it != sQuad.end(); it++){
        
    //     if ((intPointCorrespElem(index) == ielem)){

    //         //Defines the integration points adimentional coordinates
    //         xsi(0) = sQuad.PointList(index,0);
    //         xsi(1) = sQuad.PointList(index,1);
            
    //         //Computes the velocity shape functions
    //         shapeQuad.evaluate(xsi,phi_);
            
    //         xsi_intp(0) = intPointCorrespXsi(index,0);
    //         xsi_intp(1) = intPointCorrespXsi(index,1);

    //         //Computes the coarse mesh shape functions
    //         shapeQuad.evaluate(xsi_intp,phiLM_);
            
    //         //Returns the quadrature integration weight
    //         weight_ = sQuad.WeightList(index);
            
    //         //Computes the jacobian matrix
    //         getJacobianMatrix(xsi);
                        
    //         getSpatialDerivatives(xsi);

    //         dphiL_dx = dphi_dx;

    //         getSpatialDerivatives(xsi_intp);

    //         double LM = 0.;

    //         for (int i = 0; i < 6; i++){
    //             for (int j = 0; j < 6; j++){
    //                 LM = ((u_ - umesh_) * phi_(i) + (v_ - vmesh_) * phi_(i))
    //                     * phiLM_(j) * tARLQ_;
            
    //                 jacobianNRMatrix(2*i  ,2*j  ) += timeScheme_ * dTime_ * LM
    //                     * weight_ * djac_ * intPointWeightFunction(index);
                    
    //                 jacobianNRMatrix(2*i+1,2*j+1) += timeScheme_ * dTime_ * LM
    //                     * weight_ * djac_ * intPointWeightFunction(index);
                    
          
    //                 double Lx = 0.;
    //                 double Ly = 0.;
                    
    //                 if (glueZone){
    //                     if (model){
    //                         //Fine
    //                         Lx = -dphi_dx(0,i) * phi_(j) * tARLQ_ / dens_;
    //                         Ly = -dphi_dx(1,i) * phi_(j) * tARLQ_ / dens_;
    //                     }else{
    //                         //Coarse
    //                         Lx = dphi_dx(0,i) * phi_(j) * tARLQ_ / dens_ * 0.;
    //                         Ly = dphi_dx(1,i) * phi_(j) * tARLQ_ / dens_ * 0.;
    //                     };
    //                 };
                    
                    
    //                 jacobianNRMatrix(12+j,2*i  ) +=  ((Lx) * 
    //                                                   timeScheme_ * dTime_) * 
    //                     weight_ * djac_;// * intPointWeightFunction(index);
    //                 jacobianNRMatrix(12+j,2*i+1) +=  ((Ly) * 
    //                                                   timeScheme_ * dTime_)
    //                     * weight_ * djac_;// * intPointWeightFunction(index);
    //             };

    //             // double LMx = 0.; double LMy = 0.;
    //             // LMx = ((u_ - umesh_) * phi_(i) + (v_ - vmesh_) * phi_(i)) 
    //             //     * lagMx_ * tARLQ_;
    //             // LMx = ((u_ - umesh_) * phi_(i) + (v_ - vmesh_) * phi_(i)) 
    //             //     * lagMy_ * tARLQ_;   
                

    //         };
    //     };
    //     index++;
    // };
    
    return;
};










#endif

