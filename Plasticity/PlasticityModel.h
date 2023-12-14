#ifndef PLASTICITYMODEL_H
#define PLASTICITYMODEL_H

#include "DataTypes.h"
#include "WeakForm.h"

/// @brief Implements the base class for all plasticity models
class PlasticityModel : public WeakForm
{
protected: 
    // The real dimension, corresponding to the order of the stress/strain tensor
    int fRealDimension;

    // The elastic model
    WeakForm *fElasticModel;

public:
    /// @brief Default constructor
    PlasticityModel() = default;

    /// @brief Plasticity model constructor
    /// @param elast Elasticity model
    PlasticityModel(WeakForm *elast);

    /// @brief Default destructor
    ~PlasticityModel(){};

    /// @brief Verify the plastic creterion and computes the plastic strain
    /// @param data integration point data
    /// @param plasticstrain tensor of plastic strain
    /// @param totalstrain tensor of total strain
    virtual void ComputePlasticStrain(IntPointData &data, MatrixDouble &plasticstrain, MatrixDouble &totalstrain) = 0;
    
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

    WeakForm* ElasticModel(){return fElasticModel;}
};



#endif