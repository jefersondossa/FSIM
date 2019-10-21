//------------------------------------------------------------------------------
// 
//                   Jeferson W D Fernandes and Rodolfo A K Sanches
//                             University of Sao Paulo
//                           (C) 2017 All Rights Reserved
//
// <LicenseText>
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//--------------------------FLUID CONSTANT PARAMETERS---------------------------
//------------------------------------------------------------------------------

#ifndef FLUID_PARAMETERS_H
#define FLUID_PARAMETERS_H

using namespace boost::numeric;

/// Defines the fluid boundary shape functions

template<int DIM>
class FluidParameters {
public:
    
    /// Sets the element viscosity
    /// @param double element viscosity
    void setViscosity(double& visc){viscosity = visc;}

    /// Sets the element density
    /// @param double element density
    void setDensity(double& dens){density = dens;}

    /// Sets the time step size
    /// @param double time step size
    void setTimeStep(double& dt){timeStepSize = dt;}

    /// Sets the time integration scheme
    /// @param double time integration scheme: 0.0 - Explicit forward Euler;
    /// 1.0 - Implicit backward Euler;
    /// 0.5 - Implicit Trapezoidal Rule.
    void setSpectralRadius(double& b){
        spectralRadius = b;
        alpha_f = 1. / (1. + spectralRadius);
        alpha_m = 0.5 * (3. - spectralRadius) / (1. + spectralRadius);
        gamma = 0.5 + alpha_m - alpha_f;
    }

    void setFieldForce(double* ff){fieldForce = ff;}

    void setArlequinOperatorConstants(double& k_1, double& k_2){k1 = k_1; k2 = k_2;}



    double& getTimeStep() {return timeStepSize;}
    double& getDensity() {return density;}
    double& getViscosity() {return viscosity;}
    double& getAlphaM() {return alpha_m;}
    double& getAlphaF() {return alpha_f;}
    double& getGamma() {return gamma;}
    double* &getFieldForce() {return fieldForce;}
    double& getArlequinK1() {return k1;}
    double& getArlequinK2() {return k2;}

private:
    double viscosity;
    double density;
    double timeStepSize;
    double spectralRadius;
    double alpha_m;
    double alpha_f;
    double gamma;
    double* fieldForce;
    double k1, k2;
    
};






#endif
