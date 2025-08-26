#ifndef PLASTICITYMODEL_H
#define PLASTICITYMODEL_H

#include "DataTypes.h"
#include "WeakForm.h"
#include "Tensor3D.h"

/// @brief Implements the base class for all plasticity models
class PlasticityModel : public WeakForm
{
protected: 
    // The real dimension, corresponding to the order of the stress/strain tensor
    int fRealDimension;

    // The elastic model
    WeakForm *fElasticModel;

    //Number of Stress components
    int fNStressComponents;
    
    // The Young modulus
    double fYoungModulus;
    double fPoissonRatio;
    double fShearModulus;
    double fBulkModulus;
    bool fPlaneStress = false;
    double fThickness = 1.;

    double fHardening;

    // Lambda function returning the stress-stran relation
    std::function<void (const double &plasticStrain, double &yieldStress, double &hardening)> fUniaxialYield = 0; 
public:
    //2nd order identity tensor
    MatrixDouble fId2xId2;
    
    //Deviatory 4th order identity tensor
    MatrixDouble fIdentity4Dev;


public:
    /// @brief Default constructor
    PlasticityModel() = default;

    /// @brief Plasticity model constructor
    /// @param elast Elasticity model
    PlasticityModel(WeakForm *elast);

    /// @brief Default destructor
    ~PlasticityModel(){};

    /// @brief Returns the weak form stiffness matrix. It should never be called in this class
    /// @param index integration point index
    /// @param data integration point data
    /// @param Stiffness stiffness matrix
    virtual void ComputeTangentStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness, Tensor3D &Stress) {
        PanicButton();
    };


    void ComputeResidual(int &index, IntPointData &data, VecDouble &Rhs, Tensor3D &Stress) override;

    virtual double PlasticMultiplier(int &index, IntPointData &data, Tensor3D &Stress){
        PanicButton();
    };

    /// @brief Computes the principal stress Tensor3D by means of eigenvalue problem
    /// @param Stress Stress tensor
    /// @param PrincipalS principal stresses
    void ComputePrincipalStress(MatrixDouble &Stress, VecDouble &PrincipalS);

    /// @brief Returns the real dimension
    /// @return real dimension
    int &RealDimension(){return fRealDimension;}

    /// @brief Converts a Tensor3D to Voigt notation
    /// @param Tensor3D tensor
    /// @param voigt vector in Voigt notation
    void TensorToVoigt(MatrixDouble &tensor, VecDouble &voigt);

    /// @brief Converts a Tensor3D in Voigt notation to Tensor3D notation
    /// @param Tensor3D tensor
    /// @param voigt vector in Voigt notation
    void VoigtToTensor(MatrixDouble &tensor, VecDouble &voigt);

    virtual double YieldFunction(int &index, IntPointData &data, Tensor3D &Stress){
        PanicButton();
    }

    WeakForm* ElasticModel(){return fElasticModel;}

    std::function<void (const double &plasticStrain, double &yieldStress, double &hardening)> &UniaxialYield(){
        return fUniaxialYield;
    }

    void SetUniaxialYieldFunction(std::function<void (const double &plasticStrain, double &yieldStress, double &hardening)> yield){
        fUniaxialYield = yield;
    }

    int NStressComponents(){return fNStressComponents;}

    virtual void UpdateStateVariables(int &index, IntPointData &data, Tensor3D &Stress){
        PanicButton();
    };

    bool &PlaneStress(){
        return fPlaneStress;
    }

    double &YoungModulus(){
        return fYoungModulus;
    }

    double &PoissonRatio(){
        return fPoissonRatio;
    }
    double &ShearModulus(){
        return fShearModulus;
    }
    double &BulkModulus(){
        return fBulkModulus;
    }

    /// @brief Returns the variable index of a given solution variable
    /// @param name solution variable name
    /// @return solution variable's index
    int VariableIndex(const std::string &name) const override;

    /// @brief Returns the number of solution variables of a given post processing variable
    /// @param var solution variable's index
    /// @return number of solution variables
    int NSolutionVariables(int var) const override;

    /// @brief Post process the results for a given solution variable. It should never be called here, but in the derived weak form.
    /// @param data integration point data
    /// @param var solution variable's index
    /// @param Sol solution vector
    void Solution(IntPointData &data, int var, VecDouble &Sol) override;
};



#endif