#ifndef PLASTICITYMODEL_H
#define PLASTICITYMODEL_H

#include "DataTypes.h"
#include "WeakForm.h"
#include "Tensor.h"

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
    //2nd order identity tensor
    MatrixDouble Id2xId2;
    
    //Deviatory 4th order identity tensor
    MatrixDouble fIdentity4Dev;

    // The Young modulus
    double fYoungModulus;
    double fPoissonRatio;
    double fShearModulus;
    double fBulkModulus;
    bool fPlaneStress = false;

    double fHardening;

    // Lambda function returning the stress-stran relation
    std::function<void (const double &plasticStrain, double &yieldStress, double &hardening)> fUniaxialYield = 0; 

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
    virtual void ComputeTangentStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness, Tensor &Stress) {
        PanicButton();
    };


    virtual void ComputeResidual(int &index, IntPointData &data, VecDouble &Rhs, Tensor &Stress);

    virtual double PlasticMultiplier(int &index, IntPointData &data, Tensor &Stress){
        PanicButton();
    };

    /// @brief Computes the principal stress tensor by means of eigenvalue problem
    /// @param Stress Stress tensor
    /// @param PrincipalS principal stresses
    void ComputePrincipalStress(MatrixDouble &Stress, VecDouble &PrincipalS);

    /// @brief Returns the real dimension
    /// @return real dimension
    int &RealDimension(){return fRealDimension;}

    /// @brief Converts a tensor to Voigt notation
    /// @param tensor tensor
    /// @param voigt vector in Voigt notation
    void TensorToVoigt(MatrixDouble &tensor, VecDouble &voigt);

    /// @brief Converts a tensor in Voigt notation to tensor notation
    /// @param tensor tensor
    /// @param voigt vector in Voigt notation
    void VoigtToTensor(MatrixDouble &tensor, VecDouble &voigt);

    virtual double YieldFunction(int &index, IntPointData &data, Tensor &Stress){
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

    virtual void UpdateStateVariables(int &index, IntPointData &data, Tensor &Stress){
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
};



#endif