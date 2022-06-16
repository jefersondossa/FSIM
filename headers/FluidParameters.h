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

#include <math.h>
#include "DataTypes.h"

/// Defines the fluid boundary shape functions

template<int DIM, int DEG>
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

        // alpha_f = 1.;
        // alpha_m = 1.;
        // gamma = 1.;
    }

    /// Sets the domain force vector
    /// @param double* domain force vector
    void setFieldForce(double* ff){fieldForce.resize(3); for(int i=0; i<3; i++) fieldForce[i] = ff[i];}

    void setArlequinOperatorConstants(double& k_1, double& k_2){k1 = k_1; k2 = k_2;}
    void setTimeInstant(int& it){timeInstant = it;}
    /// Sets the undisturbed velocity field
    /// @param double* undisturbed velocity field
    void setVelocityInf(double* u){velocityInf.resize(3); for(int i=0; i<3; i++) velocityInf[i] = u[i];}
    /// Gets the undisturbed velocity field
    /// @return undisturbed velocity field
    VecDouble &getVelocityInf() {return velocityInf;}

    double& getTimeStep() {return timeStepSize;}
    double& getDensity() {return density;}
    double& getViscosity() {return viscosity;}
    double& getAlphaM() {return alpha_m;}
    double& getAlphaF() {return alpha_f;}
    double& getGamma() {return gamma;}
    
    /// Gets the domain force vector
    /// @return domain force vector
    VecDouble &getFieldForce() {return fieldForce;}

    double& getArlequinK1() {return k1;}
    double& getArlequinK2() {return k2;}
    double& getPi() {return pi;}
    int& getTimeInstant() {return timeInstant;}

private:
    double viscosity;
    double density;
    double timeStepSize;
    double spectralRadius;
    double alpha_m;
    double alpha_f;
    double gamma;
    VecDouble fieldForce;
    double k1, k2;
    double pi = M_PI;
    VecDouble velocityInf;

    int timeInstant;
};






#endif
