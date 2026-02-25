#ifndef CONNECT_H
#define CONNECT_H

#include "DataTypes.h"

class Connect{
private:
    int fNState;
    int fNShape;
    int fOrder;
    int64_t fIndex;
    int64_t fSequenceNumber;

    /// @brief Solution
    VecDouble fSolution;
    /// @brief Solution for the previous time/load step
    VecDouble fPrevSolution;
    /// @brief Solution first time derivative
    VecDouble fDTimeSolution;
    /// @brief Solution second time derivative
    VecDouble fDDTimeSolution;
        /// @brief Nodal mesh velocity
    VecDouble fMeshVelocity;
    /// @brief Nodal Arlequin weight function
    double    fWeightFunction;
    /// @brief Nodal Arlequin previous time step weight funtion
    double    fPrevWeightFunction;

public:
    Connect(int nstate, int nshape, int order, int64_t index, int64_t seqnum);
    // void CreateConnections();

    void AllocateTimeDerivatives();

    int GetNStateVariables(){
        return fNState;
    };

    int GetNShapeFunctions(){
        return fNShape;
    };

    int GetOrder(){
        return fOrder;
    };

    int64_t GetIndex(){
        return fIndex;
    };

    int64_t GetSequenceNumber(){
        return fSequenceNumber;
    };

    /// @brief Sets the solution
    /// @param istate state variable index
    /// @param sol solution
    void SetSolution(int istate, double sol);

    /// @brief Increments the solution of a given state variable
    /// @param istate state variable index
    /// @param sol increment
    void IncrementSolution(int istate, double sol);

    /// @brief Returns a vector with all solutions
    /// @return Solution vector
    VecDouble &Solution() {
        return fSolution;
    }
    
    /// @brief Returns the previous time/load step solution 
    /// @return Previous time/load step solution vector
    VecDouble &PrevSolution() {
        return fPrevSolution;
    }

    /// @brief Returns the solution of a given state variable
    /// @param istate state variable index
    /// @return solution
    double GetSolution(int istate) const {
        return fSolution[istate];
    }

    /// @brief Returns the nodal solution of a given state variable in the previous time/load step
    /// @param istate state variable index
    /// @return solution
    double GetPreviousSolution(int istate) const {
        return fPrevSolution[istate];
    }
    void SetPreviousSolution(int istate, double val);

    /// @brief Returns the solution first time derivative of a given state variable
    /// @param istate state variable index
    /// @return solution first time derivative
    double GetDSolutionDTime(int istate) const {
        return fDTimeSolution[istate];
    }
    VecDouble &SolutionDTime(){
        return fDTimeSolution;
    }
    void SetDSolutionDTime(int istate, double val);

    /// @brief Returns the solution second time derivative of a given state variable
    /// @param istate state variable index
    /// @return solution second time derivative
    double GetDSolutionDDTime(int istate) const {
        return fDDTimeSolution[istate];
    }
    VecDouble &SolutionDDTime(){
        return fDDTimeSolution;
    }
    void SetDSolutionDDTime(int istate, double val);


    //.........................Mesh Velocity functions..........................
    /// Sets the node mesh velocity
    /// @param double* mesh velocity
    void setMeshVelocity(VecDouble &u);
    void setMeshVelocityComponent(int dir,double u);

    /// Gets the node mesh velocity
    /// @param int direction @return mesh velocity component
    double getMeshVelocity(int dir) const {
        return fMeshVelocity[dir];
    }
   
    //............................Arlequin functions............................
    /// Sets the nodal energy weight function value
    /// @param double weight function value
    void setWeightFunction(double val);

    /// Gets the nodal energy weight function value
    /// @return weight function value
    double getWeightFunction() const {
        return fWeightFunction;
    };
    double getPreviousWeightFunction() const {
        return fPrevWeightFunction;
    };
};

#endif // CONNECT_H