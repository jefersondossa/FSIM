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
#include "IntegrationQuadrature11.h"
#include "ElCouplingLocal.h"
#include <map>
#include <set>

/// Mounts the overlapping mesh problem for solving the incompressible flow problem
class Arlequin{
public:
    std::vector<Node *>     nodesLagrangeFine_;
    std::vector<Node *>     nodesLagrangeCoarse_;

    std::vector<Element *>     glueZoneFine_;
    std::vector<Element *>     glueZoneCoarse_;

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

    std::map<int64_t,int64_t> fNodeLocalToElementGlobal;
    std::map<int64_t,VecDouble> fNodeLocalToXsiGlobal;

    std::map<int64_t,VecDouble> fLocalIntPointToGlobalElement;
    std::map<int64_t,MatrixDouble> fLocalIntPointToGlobalXsi;

    std::map<int64_t,std::set<int64_t>> fGlobalElToLocalEl;

private:
    double fGlueZoneThickness = 0.5;
    double fArlequinEpsilon = 0.e-3;

    int numElemGlueZoneFine;
    int numElemGlueZoneCoarse;
    int numNodesGlueZoneFine;
    int numNodesGlueZoneCoarse;
    int rank;

    std::map<int64_t,double> fLocalSignaledDistance;
    std::map<int64_t,double> fGlobalSignaledDistance;
    

    ArlequinStabType fArlequinStab = ArlequinStabType::ENoStab;
    std::vector<CompMesh *> fMeshVector;

public:
    Arlequin() = default;

    Arlequin(std::vector<CompMesh *> &meshvec){
        fMeshVector = meshvec;
        fMeshVector.resize(3);
        fMeshVector[2] = new CompMesh(fMeshVector[0]->getProblemParameters());
        if (fMeshVector[2]->getProblemParameters().ProbType() == ProblemType::ENavierStokes || fMeshVector[2]->getProblemParameters().ProbType() == ProblemType::EStokes){
            fMeshVector[2]->SetNStateVariables(fMeshVector[2]->Dimension());
        }
    }

    std::vector<CompMesh *> &MeshVec(){return fMeshVector;}

    void setArlequinStabilization(ArlequinStabType stab){fArlequinStab = stab;}
    ArlequinStabType &getArlequinStabilization(){return fArlequinStab;}

    /// Sets the coarse and mesh models. It is considered that the fine model
    /// is completely immersed on the coarse model.
    /// @param Fluid coarse model @param Fluid fine model
    // void setFluidModels(FluidMesh& coarse, FluidMesh& fine);

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
    void SetElementBoxes();

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
    void searchNodeCorrespondence(VecDouble &x, CompMesh *cmesh,
                                  int &elCorr, VecDouble &xsiCorr, int elSearch);

    void setMatVecValuesFineModel(MatrixDouble &matrix, VecDouble &rhs, VecInt &connec);
    void setMatVecValuesFineModelPoisson(MatrixDouble &matrix, VecDouble &rhs, VecInt &connec);
    void setMatVecValuesFineModelElasticity(MatrixDouble &matrix, VecDouble &rhs, VecInt &connec);
    void setMatVecValuesCoarseModel(MatrixDouble &matrix, VecDouble &rhs, VecInt &connec);
    void setMatVecValuesCoarseModelPoisson(MatrixDouble &matrix, VecDouble &rhs, VecInt &connec);
    void setMatVecValuesCoarseModelElasticity(MatrixDouble &matrix, VecDouble &rhs, VecInt &connec);
    void setMatVecValuesLagMultFineFine(MatrixDouble &Ajac2, MatrixDouble &localMV_mat, 
                                        MatrixDouble &ArlequinA1, MatrixDouble &ArlequinA2, 
                                        VecDouble &Rhs2, VecDouble &rhsLagMult2,
                                        VecDouble &localMV_vec, VecDouble &RhsArlequin2,
                                        VecInt &connec, VecInt &connecL);
    void setMatVecValuesLagMultFineFinePoisson(MatrixDouble &Ajac2, MatrixDouble &localMV_mat, 
                                               MatrixDouble &ArlequinA1, MatrixDouble &ArlequinA2, 
                                               VecDouble &Rhs2, VecDouble &rhsLagMult2,
                                               VecDouble &localMV_vec, VecDouble &RhsArlequin2,
                                               VecInt &connec, VecInt &connecL);
    void setMatVecValuesLagMultFineFineElasticity(MatrixDouble &Ajac2, MatrixDouble &localMV_mat, 
                                               MatrixDouble &ArlequinA1, MatrixDouble &ArlequinA2, 
                                               VecDouble &Rhs2, VecDouble &rhsLagMult2,
                                               VecDouble &localMV_vec, VecDouble &RhsArlequin2,
                                               VecInt &connec, VecInt &connecL);

    void setMatVecValuesLagMultFineCoarse(MatrixDouble &Ajac2, MatrixDouble &localMV_mat, 
                                          MatrixDouble &ArlequinA1, MatrixDouble &ArlequinA2, 
                                          VecDouble &Rhs2, VecDouble &rhsLagMult2,
                                          VecDouble &localMV_vec, VecDouble &RhsArlequin2,
                                          VecInt &connecC, VecInt &connecL);
    void setMatVecValuesLagMultFineCoarsePoisson(MatrixDouble &Ajac2, MatrixDouble &localMV_mat, 
                                                 MatrixDouble &ArlequinA1, MatrixDouble &ArlequinA2, 
                                                 VecDouble &Rhs2, VecDouble &rhsLagMult2,
                                                 VecDouble &localMV_vec, VecDouble &RhsArlequin2,
                                                 VecInt &connecC, VecInt &connecL);
    void setMatVecValuesLagMultFineCoarseElasticity(MatrixDouble &Ajac2, MatrixDouble &localMV_mat, 
                                                 MatrixDouble &ArlequinA1, MatrixDouble &ArlequinA2, 
                                                 VecDouble &Rhs2, VecDouble &rhsLagMult2,
                                                 VecDouble &localMV_vec, VecDouble &RhsArlequin2,
                                                 VecInt &connecC, VecInt &connecL);

    void assembleArlequinSystem();
    void assembleCoarseModel();
    void assembleFineModel();
    void assembleCouplingOperator();

    void assembleArlequinSystemPoisson();

    /// Print the results for Paraview post-processing
    /// @param int time step
    void printResultsCoarse(int step);
    void printResultsFine(int step);
    
    void stabilizeArlequin(MatrixDouble &A0, MatrixDouble &A1, 
                           MatrixDouble &C0, MatrixDouble &C1,
                           MatrixDouble &E, VecDouble &b0, 
                           VecDouble &b1, double &tArlq0, double &tArlq1);

    void stabilizeArlequin(std::vector<MatrixDouble> &Stiffness, int64_t &element);

    void CreateGlobalCouplingElements();

    void SetUp(){
        SetElementBoxes();
        setSignaledDistance();
        //Construct the glue zone based on some defined criterion
        setCouplingZone();
        //Computes the Weight function for all the finite elements
        setWeightFunction(16.);
        PanicButton();
        // for (auto el:fMeshVector[2]->ElementVec()){
        //     el->getIntegPointCoordinates();
        // }
        setNodalCorrespondenceFine();
        CreateGlobalCouplingElements();
        // DeleteCoarseEls();
    };

    void DeleteCoarseEls();

    double GlobalWeightFunction(double dist);
    double LocalWeightFunction(double dist);

    void ComputeTriangleIntersections(int64_t iEl, std::set<int64_t> &elIntersected);

};


#endif
