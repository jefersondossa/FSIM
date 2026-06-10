#ifndef WEAKFORM_H
#define WEAKFORM_H

#include "DataTypes.h"
#include "IntPointData.h"
#include "Node.h"

enum class WeakFormType
{
    kNone = 0,
    kL2Projection = 1,
    kElasticTruss = 2,
    kElasticity2D = 3,
    kElasticity3D = 4,
    kElasticityPositional2D = 5,
    kLinearBeam = 6,
    kLinearFrame = 7,
    kMixedElasticity = 8,
    kPositionalFrame2D = 9,
    kPositionalTruss = 10,
    kTransientElasticity2D = 11,
    kTransientElasticTruss = 12,
    kTransientPositional2D = 13,
    kTransientPositionalFrame2D = 14,
    kTransientPositionalTruss = 15,
    kNavierStokes = 16,
    kTransientNavierStokes = 17,
    kPoisson = 18,
    kTransientPoisson = 19,
    kStokes = 20,
    kGlobalLocalEnrichment = 21,
    kMixedGlobalLocalEnrichment = 22,
    kInterpolatedBC = 23,
    kCouplingLocal = 24,
    kCouplingGlobal = 25
};



/// @brief Implements a base for all Weak Forms 
class WeakForm
{
protected:
    // The material id, i.e. physicaltag from gmsh
    int fMatId;
    // Number of state variables
    int fNState = 1;
    // Problem dimension
    int fDimension = 0;
    // Big number to enforce boundary conditions
    static double fBigNumber;
    // Lambda function returning the exact solution and its gradient
    std::function<void (const VecDouble &coord, VecDouble &u, MatrixDouble &gradU)> fExactSol = 0; 
    // Lambda function returning the forcing function
    std::function<void (const VecDouble &coord, VecDouble &force)> fForceFunction = 0; 
    // If elemental stiffness matrix can be reused accross multiple runs.
    bool fHasMemory = true;
    // Weak Form type
    WeakFormType fType = WeakFormType::kNone;
public:

    /// @brief Default constructor
    WeakForm() = default;
    WeakForm(int matid, int nstate) : fMatId(matid), fNState(nstate) {};

    /// @brief Default destructor
    virtual ~WeakForm() = default;

    /// @brief Returns the weak form stiffness matrix. It should never be called in this class
    /// @param index integration point index
    /// @param data integration point data
    /// @param Stiffness stiffness matrix
    virtual void ComputeStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness) {
        PanicButton();
    };
    virtual void ComputeStiffness(int &index, std::vector<IntPointData *> &data, MatrixDouble &Stiffness) {
        PanicButton();
    };
    
    /// @brief Overloads the weak form stiffness matrix computation in the case more than one contribution is provided
    /// @param index integration point index
    /// @param data integration point data
    /// @param Stiffness vector of stiffness matrices
    virtual void ComputeStiffness(int &index, IntPointData &data, std::vector<MatrixDouble> &Stiffness) {
        PanicButton();
    };
    
    /// @brief Returns the weak form residual vector. It should never be called in this class
    /// @param index integration point index
    /// @param data integration point data
    /// @param Rhs residual vector
    virtual void ComputeResidual(int &index, IntPointData &data, VecDouble &Rhs) {
        PanicButton();
    };
    virtual void ComputeResidual(int &index, std::vector<IntPointData *> &data, VecDouble &Rhs) {
        PanicButton();
    };

    virtual void ComputeResidual(int &index, IntPointData &data, VecDouble &Rhs, Tensor3D &Stress) {
        PanicButton();
    };

    /// @brief Overloads the weak form residual vector computation in the case more than one contribution is provided
    /// @param index integration point index
    /// @param data integration point data
    /// @param Rhs vector of residual vectors
    virtual void ComputeResidual(int &index, IntPointData &data, std::vector<VecDouble> &Rhs) {};

    
    /// @brief Computes the element error. The exact solution shoul be provided.
    /// @param data integration point data
    /// @param errors vector storing all errors
    virtual void ComputeError(IntPointData &data, VecDouble &errors) {};
    virtual void ComputeError(std::vector<IntPointData *> &data, VecDouble &errors) {};
    
    /// @brief Sets the exact solution
    /// @param exSol Lambda function to compute the exact solution and its gradient for a given coordinate
    void SetExactSolution(std::function<void (const VecDouble &coord, VecDouble &u, MatrixDouble &gradU)> exSol){
        fExactSol = exSol;
    }

    /// @brief Gets the exact solution Lambda function
    /// @return Exact solution Lambda function
    std::function<void (const VecDouble &coord, VecDouble &u, MatrixDouble &gradU)> &GetExactSolution(){
        return fExactSol;
    }

    /// @brief Sets the source term Lambda function
    /// @param ffunction Lambda function to compute the source term for given coordinate.
    void SetForcingFunction(std::function<void (const VecDouble &coord, VecDouble &force)> ffunction){
        fForceFunction = ffunction;
    }

    /// @brief Gets the source term Lambda function
    /// @return Source term Lambda function
    std::function<void (const VecDouble &coord, VecDouble &force)> &GetForcingFunction(){
        return fForceFunction;
    }

    /// @brief Gets the Weak form physical tag
    /// @return Physical tag
    int &Id() {return fMatId;}

    /// @brief Gets the number of state variables
    /// @return Number of State variables
    int &NState() {return fNState; }

    /// @brief 
    /// @return 
    int &Dimension() {return fDimension;}

    /// @brief Returns the variable index of a given solution variable
    /// @param name solution variable name
    /// @return solution variable's index
    virtual int VariableIndex(const std::string &name) const {
        PanicButton();
        return 0;
    };

    /// @brief Returns the number of solution variables of a given post processing variable
    /// @param var solution variable's index
    /// @return number of solution variables
    virtual int NSolutionVariables(int var) const {
        PanicButton();
        return 0;
    };

    /// @brief Post process the results for a given solution variable. It should never be called here, but in the derived weak form.
    /// @param data integration point data
    /// @param var solution variable's index
    /// @param Sol solution vector
    virtual void Solution(IntPointData &data, int var, VecDouble &Sol) {
        PanicButton();
    };
    virtual void Solution(std::vector<IntPointData *> &data, int var, VecDouble &Sol) {
        PanicButton();
    };

    virtual MatrixDouble &ConstitutiveMatrix(){
        PanicButton();
    }

    void SetNStateVariables(int nst){
        fNState = nst;
    }

    bool GetHasMemory() const {
        return fHasMemory;
    }

    void SetHasMemory(bool new_has_memory) {
        fHasMemory = new_has_memory;
    }

    WeakFormType Type() const {
        return fType;
    }
};


#endif
