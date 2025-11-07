//------------------------------------------------------------------------------
// 
//                   Jeferson W D Fernandes and Rodolfo A K Sanches
//                             University of Sao Paulo
//                           (C) 2017 All Rights Reserved
//
// <LicenseText>
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//----------------------------------ARLEQUIN------------------------------------
//------------------------------------------------------------------------------

#ifndef ARLEQUIN_H
#define ARLEQUIN_H

#include "Fluid.h"
#include "Glue.h"
#include <petscsnes.h>

/// Mounts the overlapping mesh problem for solving the incompressible flow problem

template<int DIM, int DEG>
class Arlequin{
public:
    /// Defines the class Fluid locally
    typedef Fluid<DIM,DEG>                 FluidMesh;

    /// Defines the class Element locally
    typedef typename FluidMesh::Elements   Elements;

    /// Defines the class Node locally
    typedef typename FluidMesh::Node       Nodes;

    /// Defines the class Boundary locally
    typedef typename FluidMesh::Boundaries Boundary;

    /// Defines the class SpecialQuad locally
    typedef typename Elements::SpecialQuad Quadrature;

    /// Defines the class Glue locally
    typedef Glue<DIM,DEG>                  GlueZone;

    typedef FluidParameters<DIM,DEG>       Parameters;

    FluidMesh coarseModel, fineModel;

    std::vector<Nodes *>     *nodesCoarse_;
    std::vector<Nodes *>     *nodesFine_;
    std::vector<Nodes *>     nodesLagrangeFine_;
    std::vector<Nodes *>     nodesLagrangeCoarse_;

    std::vector<Elements *>  elementsCoarse_;
    std::vector<Elements *>  elementsFine_;
    std::vector<GlueZone *>  glueZoneFine_;
    std::vector<GlueZone *>  glueZoneCoarse_;

    std::vector<Boundary *>  boundaryCoarse_;
    std::vector<Boundary *>  boundaryFine_;

    std::vector<int>         elementsGlueZoneFine_;
    std::vector<int>         nodesGlueZoneFine_;
    std::vector<int>         elementsGlueZoneCoarse_;
    std::vector<int>         nodesGlueZoneCoarse_;

    ProblemType fProbType;

    Mat               A, F;
    Vec               b, u, All;
    PetscErrorCode    ierr;
    PetscInt          Istart, Iend, Ii, Ione, iterations;
    KSP               ksp; 
    PC                pc;
    VecScatter        ctx;
    PetscScalar       val;
    PetscLogDouble bytes = 0;
    SNES              snes;

public:
    int numElemCoarse;
    int numElemFine;
    int numBoundElemCoarse;
    int numBoundElemFine;
    int numElemGlueZoneFine;
    int numElemGlueZoneCoarse;
    int numNodesCoarse;
    int numNodesFine;
    int numNodesGlueZoneFine;
    int numNodesGlueZoneCoarse;
    int numTimeSteps;
    double dTime;
    int rank;
    int iTimeStep;

    Parameters *parametersCoarse, *parametersFine;

    std::pair<idx_t*,idx_t*> domDecompCoarse;//Coarse Model Domain Decomposition
    std::pair<idx_t*,idx_t*> domDecompFine;  //Fine Model Domain Decomposition

    Quadrature quad;

    double pi = M_PI;

    double alpha_f;
    double alpha_m;
    double gamma;

    int nElNodes = (3+(DIM-2)*DEG)*(2+3*DEG+DEG*DEG)/6;
    int nLocDOF = -8*DIM -21*DEG + 15*DIM*DEG + 16;
    int nBdNodes = 3*(1-DEG)+DIM*(2*DEG-1);

    ArlequinStabType fArlequinStab = ArlequinStabType::ENoStab;

public:

    void setArlequinStabilization(ArlequinStabType stab){fArlequinStab = stab;}
    ArlequinStabType &getArlequinStabilization(){return fArlequinStab;}

    /// Sets the coarse and mesh models. It is considered that the fine model
    /// is completely immersed on the coarse model.
    /// @param Fluid coarse model @param Fluid fine model
    void setFluidModels(FluidMesh& coarse, FluidMesh& fine);

    /// Mounts and solve the incompressible flow problem with overlapping meshes
    /// using the Arlequin method whit the gluing zone defined in the fine model
    /// @param int maximum number of iterations of the Newton-Raphson's process
    /// @param double tolerance of the Newton-Raphson's process
    /// @param int type of problem to be solved: 1 - Stokes; 2 - Navier-Stokes.
    /// @param int 0 - Steady problem; 1 - Transient problem.
    int solveArlequinProblem(int iterNumber, double tolerance,
                             int problem_type, int time_dependency);

    /// Mounts and solve the incompressible flow problem with overlapping meshes
    /// using the Arlequin method with the gluing zone defined in the fine model
    /// and the fine model can be moved arbitrarily in an ALE description 
    /// framework
    /// @param int maximum number of iterations of the Newton-Raphson's process
    /// @param double tolerance of the Newton-Raphson's process
    /// @param int type of problem to be solved: 1 - Stokes; 2 - Navier-Stokes.
    /// @param int 0 - Steady problem; 1 - Transient problem.
    int solveArlequinProblemMoving(int iterNumber, double tolerance,
                                   int problem_type, int time_dependency);

    /// Mounts and solve the incompressible flow problem with overlapping meshes
    /// using the Arlequin method with the gluing zone defined in the fine model
    /// and the fine model can be moved arbitrarily in an ALE description 
    /// framework for the FSI interaction problem
    /// @param int maximum number of iterations of the Newton-Raphson's process
    /// @param double tolerance of the Newton-Raphson's process
    /// @param int type of problem to be solved: 1 - Stokes; 2 - Navier-Stokes.
    int solveFSIArlequin(int iterNumber, double tolerance,
                         int problem_type, int iTimeStep);

    /// Compute and store the element boxes for improving the correspondence 
    /// searching process
    void setElementBoxes();

    /// Defines the gluing (or coupling) zone
    void setCouplingZone();

    /// Compute and store the signaled distance function
    void setSignaledDistance();

    /// Sets the energy weight function 
    /// @param double reference value for computing the energy weight function
    void setWeightFunction(double val);

    /// Sets the nodal correspondence of the fine to the coarse models
    void setNodalCorrespondenceFine();

    /// Compute and print drag and lift coefficients
    void dragAndLiftCoefficients(std::ofstream& dragLift);

    /// Searchs the point correspondence in a fluid model
    /// @param VecLocD point
    /// @param vector<Nodes> vector of fluid model nodes
    /// @param vector<Elements> vector of fluid model elements
    /// @param int number of elements of the fluid model
    void searchNodeCorrespondence(VecDouble &x, std::vector<Nodes *> nodes,
                                  std::vector<Elements *> elements,
                                  int numElem, int &elCorr, VecDouble &xsiCorr, int elSearch);

    void setMatVecValuesFineModel(MatrixDouble &matrix, VecDouble &rhs, VecInt &connec);
    void setVecValuesFineModel(VecDouble &rhs, VecInt &connec);
    void setMatValuesFineModel(MatrixDouble &matrix, VecInt &connec);
    void setMatVecValuesFineModelPoisson(MatrixDouble &matrix, VecDouble &rhs, VecInt &connec);
    void setMatVecValuesCoarseModel(MatrixDouble &matrix, VecDouble &rhs, VecInt &connec);
    void setVecValuesCoarseModel(VecDouble &rhs, VecInt &connec);
    void setMatValuesCoarseModel(MatrixDouble &matrix, VecInt &connec);
    void setMatVecValuesCoarseModelPoisson(MatrixDouble &matrix, VecDouble &rhs, VecInt &connec);
    void setMatVecValuesLagMultFineFine(MatrixDouble &Ajac2, MatrixDouble &localMV_mat, 
                                        MatrixDouble &ArlequinA1, MatrixDouble &ArlequinA2, 
                                        VecDouble &Rhs2, VecDouble &rhsLagMult2,
                                        VecDouble &localMV_vec, VecDouble &RhsArlequin2,
                                        VecInt &connec, VecInt &connecL);
    void setVecValuesLagMultFineFine(VecDouble &Rhs2, VecDouble &rhsLagMult2,
                                     VecDouble &localMV_vec, VecDouble &RhsArlequin2,
                                     VecInt &connec, VecInt &connecL);
    void setMatValuesLagMultFineFine(MatrixDouble &Ajac2, MatrixDouble &localMV_mat, 
                                     MatrixDouble &ArlequinA1, MatrixDouble &ArlequinA2, 
                                     VecInt &connec, VecInt &connecL);

    void setMatVecValuesLagMultFineFinePoisson(MatrixDouble &Ajac2, MatrixDouble &localMV_mat, 
                                               MatrixDouble &ArlequinA1, MatrixDouble &ArlequinA2, 
                                               VecDouble &Rhs2, VecDouble &rhsLagMult2,
                                               VecDouble &localMV_vec, VecDouble &RhsArlequin2,
                                               VecInt &connec, VecInt &connecL);

    void setMatVecValuesLagMultFineCoarse(MatrixDouble &Ajac2, MatrixDouble &localMV_mat, 
                                          MatrixDouble &ArlequinA1, MatrixDouble &ArlequinA2, 
                                          VecDouble &Rhs2, VecDouble &rhsLagMult2,
                                          VecDouble &localMV_vec, VecDouble &RhsArlequin2,
                                          VecInt &connecC, VecInt &connecL);
    void setVecValuesLagMultFineCoarse(VecDouble &Rhs2, VecDouble &rhsLagMult2,
                                       VecDouble &localMV_vec, VecDouble &RhsArlequin2,
                                       VecInt &connecC, VecInt &connecL);
    void setMatValuesLagMultFineCoarse(MatrixDouble &Ajac2, MatrixDouble &localMV_mat, 
                                       MatrixDouble &ArlequinA1, MatrixDouble &ArlequinA2, 
                                       VecInt &connecC, VecInt &connecL);
    
    void setMatVecValuesLagMultFineCoarsePoisson(MatrixDouble &Ajac2, MatrixDouble &localMV_mat, 
                                                 MatrixDouble &ArlequinA1, MatrixDouble &ArlequinA2, 
                                                 VecDouble &Rhs2, VecDouble &rhsLagMult2,
                                                 VecDouble &localMV_vec, VecDouble &RhsArlequin2,
                                                 VecInt &connecC, VecInt &connecL);

    void assembleArlequinSystem();
    void assembleArlequinMatrix();
    void assembleArlequinVector();
    void assembleArlequinSystemPoisson();

    /// Print the results for Paraview post-processing
    /// @param int time step
    void printResultsCoarse(int step);
    void printResultsFine(int step);

    void computeErrorPoisson();

    void stabilizeArlequin(MatrixDouble &A0, MatrixDouble &A1, 
                           MatrixDouble &C0, MatrixDouble &C1,
                           MatrixDouble &E, VecDouble &b0, 
                           VecDouble &b1, double &tArlq0, double &tArlq1);

    void stabilizeArlequinNew(VecDouble &Ml1, VecDouble &t1, 
                             VecDouble &j1, VecDouble &k1, VecDouble &p1, double &tArleq);

    static PetscErrorCode FormFunction(SNES snes, Vec vecU,Vec vecB, void *ptr);
    static PetscErrorCode FormJacobian(SNES snes,Vec vecU,Mat matA, Mat matB, void *ptr);

};


#endif
