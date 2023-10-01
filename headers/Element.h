
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
#include "ProblemParameters.h"
#include "PanicButton.h"

#include "Boundary.h"
#include "CompMesh.h"
#include "IntegrationQuadrature.h"
#include "IntegrationQuadrature11.h"
#include "DomainIntegration.h"
// #include "PartitionedQuadrature.hpp"

/// Defines the fluid element object and all the element information
class Element{
private:
    VecInt        fConnect; //Velocity mesh connectivity 
    int64_t       fIndex;             //Element index
    CompMesh *fMesh;
    VecDouble     xK, XK;
    int           fSideInBoundary;
    std::vector<int64_t> fNeighborElements;
    
public:
    VecDouble intPointWeightFunction;
    VecDouble intPointWeightFunctionPrev;

    MatrixDouble fIntPointCoordinates;

    bool          FSIInterface;    
    int DIM, DEG;
public:
    Element() = default;

    /// fluid element constructor
    /// @param int element index @param Connectivity element connectivity
    /// @param vector<Nodes> 
    Element(int64_t index, VecInt &connect, CompMesh* mesh){
        
        fMesh = mesh;
        fConnect.resize(fMesh->NElNodes());
        fIndex = index;
        for (int i = fMesh->NElNodes(); i--; ) fConnect[i] = connect[i];
        DIM = fMesh->Dimension();
        DEG = fMesh->GetDefaultOrder();

        FSIInterface = false;
        fSideInBoundary = -1;
        fNeighborElements.clear();

        IntegQuadrature nQuad(fMesh->Dimension(),fMesh->GetDefaultOrder());
        intPointWeightFunction.resize(nQuad.getNumberOfIntegrationPoints());
        intPointWeightFunctionPrev.resize(nQuad.getNumberOfIntegrationPoints());

        intPointWeightFunction.fill(1.);
        intPointWeightFunctionPrev.fill(1.);
    
        getIntegPointCoordinates();

    };

    CompMesh* Mesh() {return fMesh;}
    void SetMesh(CompMesh* mesh){fMesh = mesh;}

    //........................Element basic information.........................
    /// Clear all element variables
    void clearVariables();

    /// Sets the element connectivity
    /// @param int* element connectivity
    void setConnectivity(VecInt &connect){fConnect = connect;};

    /// Gets the element connectivity
    /// @return element connectivity
    VecInt &getConnectivity(){return fConnect;};

    /// Compute and store the spatial jacobian matrix
    /// @param bounded_vector integration point adimensional coordinates
    void getJacobianMatrix(VecDouble &xsi, MatrixDouble &ainv_, double &djac_, int index);

    /// Compute and store the shape function spatial derivatives
    /// @param bounded_vector integration point adimensional coordinates
    void getSpatialDerivatives(VecDouble &xsi, MatrixDouble &ainv_, MatrixDouble &dphi_dx);
    void getHighOrderSpatialDerivatives(VecDouble &xsi, MatrixDouble &ainv_, MatrixDouble &dphi_dx, MatrixDouble &dDphi_dx);

    void interpolateSolution(int &index, VecDouble &u_);
    void interpolateMeshVelocity(int &index, VecDouble &umesh_, VecDouble &umeshPrev_);
    void interpolateSolDerivatives(MatrixDouble &dphi_dx, MatrixDouble &du_dx);

    /// Compute and store the SUPG, PSPG and LSIC stabilization parameters
    void getParameterArlequin(int &index, double &tARLQ_, double &tSUPG_, double &tPSPG_, double &tLSIC_, MatrixDouble &dphi_dx);

    /// Gets the element jacobian determinant
    /// @return element jacobinan determinant
    double getJacobian(){
        VecDouble xsi(fMesh->Dimension());
        MatrixDouble ainv_(fMesh->Dimension(),fMesh->Dimension());
        MatrixDouble dphi_dx(fMesh->NElNodes(),fMesh->Dimension());
        
        ShapeFunction  shapeQuad(fMesh->Dimension(),fMesh->GetDefaultOrder());

        xsi[0] = 0.5;
        xsi[1] = 0.5;
        
        // std::cout << "AAA 1 "<< std::endl;

        double djac_ = 0.;
        //Computes the jacobian matrix
        getJacobianMatrix(xsi, ainv_, djac_, 0);
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

    /// Pushs back a term of the inverse incidence, i.e., an element which
    /// contains the node
    /// @param int element
    void pushNeighborElement(int el) {
        fNeighborElements.push_back(el);
    }

    /// Gets the number of elements which contains the node
    /// @return int number of elements which contains the node
    int getNumberOfNeighborElements(){
        return fNeighborElements.size();
    }

    /// Gets an specific member of the inverse incidence
    /// @param int index @return int element of the inverse incidence
    int64_t &getNeighborElement(int i){
        return fNeighborElements[i];
    }

    void clearInverseIncidence(){
        fNeighborElements.clear();
        fNeighborElements.shrink_to_fit();
    }

    int64_t &Index(){return fIndex;}

    void sortEraseNeighborElements(){
        std::sort(fNeighborElements.begin(), fNeighborElements.end());
        fNeighborElements.erase(std::unique(fNeighborElements.begin(),fNeighborElements.end()), fNeighborElements.end());
    }

    /// Sets the element side in boundary
    /// @param int side in boundary
    void setElemSideInBoundary(int side){fSideInBoundary = side;};

    /// Gets the element side in boundary
    /// @return side in boundary
    int &getElemSideInBoundary(){return fSideInBoundary;};

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
    /// Sets if the element belongs to the fluid structure interface
    void setFSIInterface(){FSIInterface = true;};
    bool &getFSIInterface(){return FSIInterface;};

    //......................Integration Points Information......................
    /// Gets the number of integration points of the special quadrature rule
    /// @retunr number of integration point of the special quadrature rule
    int getNumberOfIntegrationPoints(){IntegQuadrature sQuad(fMesh->Dimension(),fMesh->GetDefaultOrder()); return sQuad.getNumberOfIntegrationPoints();};

    /// Compute and store the integration points global coordinates
    void getIntegPointCoordinates();

    /// Gets the integration point global coordinates
    /// @param int integration point index @return integration point coordinates
    VecDouble getIntegPointCoordinatesValue(int index){
        VecDouble aux(2);
        aux[0] = fIntPointCoordinates(index,0); 
        aux[1] = fIntPointCoordinates(index,1); 
        return aux;
    };

    /// Sets the integration point energy weight function
    /// @param int integration point index 
    /// @param double energy weight function value
    void setIntegPointWeightFunction();

    /// Gets the integration point energy weight function
    /// @param int integration point index @return energy weight function value
    double &getIntegPointWeightFunction(int index)
    {return intPointWeightFunction[index];};

    //.......................Element vectors and matrices.......................
    /// Compute and store the element matrix for the incompressible flow problem
    /// @param int integration point index
    void getElemMatrix(int &index, MatrixDouble &dphi_dx, double &tSUPG_, double &tPSPG_, double &tLSIC_, double &weight_, double &djac_, MatrixDouble &jacobianNRMatrix);

    /// Compute and store the element matrix for the Laplace/Poisson problem
    void getElemLaplMatrix(double &weight_, double &djac_, MatrixDouble &dphi_dx, MatrixDouble &jacobianNRMatrix);
    void getElemElasticity2DMatrix(int &index, double &weight_, double &djac_, MatrixDouble &dphi_dx, MatrixDouble &jacobianNRMatrix);

    /// Sets the boundary conditions for the incompressible flow problem
    void setBoundaryConditions(MatrixDouble &jacobianNRMatrix, VecDouble &rhsVector);

    /// Sets the boundary conditions for the Laplace/Poisson problem
    void setBoundaryConditionsLaplace(MatrixDouble &jacobianNRMatrix, VecDouble &rhsVector);

    ///Compute and store the residual vector for the incompressible flow problem
    /// @param int integration point index
    void getResidualVector(int &index, MatrixDouble &dphi_dx, double &tSUPG_, double &tPSPG_, double &tLSIC_, double &weight_, double &djac_, VecDouble &rhsVector);
    void getResidualVectorElasticity2D(int &index, MatrixDouble &dphi_dx, double &weight_, double &djac_, VecDouble &rhsVector);

    /// Compute and store the residual vector for the Laplace/Poisson problem
    void getResidualVectorLaplace(VecDouble &rhsVector);

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
    void getLagrangeMultipliersSameMesh(MatrixDouble &lagrMultMatrix, VecDouble &lagrMultVector, VecDouble &rhsVector);
    void getLagrangeMultipliersSUPG_PSPG_SameMesh(MatrixDouble &jacobianNRMatrix, VecDouble &rhsVector);
    void getLagrangeMultipliersArlequinSameMesh(MatrixDouble &arlequinStab, MatrixDouble &laplMatrix, VecDouble &arlequinStabVector);

    /// Compute and store the Lagrange multiplier operator when integrationg
    /// the different mesh portion
    /// @param int element of the coarse mesh (used to verify which integration
    /// point belongs to the coarse mesh element)
    void getLagrangeMultipliersDifferentMesh(int &ielem, double &tPSPG2_,VecDouble &press, VecDouble &velx, VecDouble &vely,
                                             VecDouble &velxPrev, VecDouble &velyPrev, MatrixDouble &lagrMultMatrix, VecDouble &rhsVectorLM, VecDouble &rhsVector);
    void getLagrangeMultipliersSUPG_PSPG_DifferentMesh(int &ielem, double &tPSPG2_,VecDouble &press, VecDouble &velx, VecDouble &ely,
                                                       MatrixDouble &jacobianNRMatrix, VecDouble &rhsVector);
    void getLagrangeMultipliersArlequinDifferentMesh(int &ielem, double &tPSPG2_, VecDouble &press, VecDouble &velx, VecDouble &vely,
                                                     MatrixDouble &arlequinStab, MatrixDouble &laplMatrix, VecDouble &arlequinStabVector);

    //...............................Problem type...............................
    /// Compute the Transient Navier-Stokes problem matrices and vectors
    void getTransientNavierStokes(MatrixDouble &jacobianNRMatrix, VecDouble &rhsVector);
    
    /// Compute the Transient Navier-Stokes problem matrices and vectors
    void getPoisson(MatrixDouble &matrix, VecDouble &rhsVector);

    /// Compute the Steady Laplace problem matrices and vectors 
    /// (usually for the mesh moving step)
    void getSteadyLaplace(MatrixDouble &jacobianNRMatrix, VecDouble &rhsVector);
    void getSolidProblem(MatrixDouble &jacobianNRMatrix, VecDouble &rhsVector);
    void getElasticity2D(MatrixDouble &jacobianNRMatrix, VecDouble &rhsVector);

    void computeErrorPoisson(VecDouble &errors);





    void ComputeElContribution(MatrixDouble &Stiffness, VecDouble &Rhs);
    virtual void ApplyBC(MatrixDouble &Stiffness, VecDouble &Rhs){};
    virtual void ComputeStiffness(int &index, MatrixDouble &dphi_dx, double &weight_, double &djac_, MatrixDouble &Stiffness){};
    virtual void ComputeResidual(int &index, MatrixDouble &dphi_dx, double &weight_, double &djac_, VecDouble &Rhs){};
    virtual void ComputeError(VecDouble &errors){};
    


};


#endif

