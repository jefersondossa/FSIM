//------------------------------------------------------------------------------
// 
//                   Jeferson W D Fernandes and Rodolfo A K Sanches
//                             University of Sao Paulo
//                           (C) 2017 All Rights Reserved
//
// <LicenseText>
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//--------------------------QUADRATIC SHAPE FUNCTION----------------------------
//------------------------------------------------------------------------------

#ifndef QUADSHAPEFUNCTION_H
#define QUADSHAPEFUNCTION_H

using namespace boost::numeric;

/// Defines the quadratic shape functions and its derivatives

template<int DIM>

class QuadShapeFunction {
public:
    
    /// Evaluates the shape function value
    /// @param Coords Adimensional coordinates 
    /// @param Values Shape function values
    void evaluate(double xi[], double phi[]) const;
    
    /// Evaluates the values of the shape funtion derivatives
    /// @param Coords Adimensional coordinates 
    /// @param ValueDeriv Shape function derivatives values
    void evaluateGradient(double xi[], double dphi[][6]) const;   

    /// Evaluates the values of the shape funtion second derivatives    
    /// @param Coords Adimensional coordinates 
    /// @param ValueDeriv Shape function second derivatives values
    void evaluateHessian(double xi[], double ***ddphi) const;   
};


//------------------------------------------------------------------------------
//-------------------------COMPUTE SHAPE FUNCTION VALUE-------------------------
//------------------------------------------------------------------------------
// Defines quadratic shape functions and its derivatives 
// for triangles and tetrahedrons
template<>
void QuadShapeFunction<2>::evaluate(double xi[], double phi[]) const {
    
    double xsi1 = xi[0];
    double xsi2 = xi[1];
    double xsi3 = 1. - xsi1 - xsi2;

    phi[0] = xsi3 * (2.0 * xsi3 - 1.0);
    phi[1] = xsi1 * (2.0 * xsi1 - 1.0);
    phi[2] = xsi2 * (2.0 * xsi2 - 1.0);
    phi[3] = 4.0 * xsi3 * xsi1;
    phi[4] = 4.0 * xsi1 * xsi2;
    phi[5] = 4.0 * xsi2 * xsi3;

     // element conectivity
     //     2
     //     54
     //     031
     
     return;
}

//------------------------------------------------------------------------------
//-------------------COMPUTE SHAPE FUNCTION DERIVATIVE VALUE--------------------
//------------------------------------------------------------------------------
template<>
void QuadShapeFunction<2>::evaluateGradient(double xi[], double dphi[][6]) const {

    const double xsi1 = xi[0];
    const double xsi2 = xi[1];
    const double xsi3 = 1. - xsi1 - xsi2;

    dphi[0][1] = 4. * xsi1 - 1.;
    dphi[1][1] = 0.;

    dphi[0][2] = 0.;
    dphi[1][2] = 4. * xsi2 - 1.;

    dphi[0][0] = -4. * xsi3 + 1.;
    dphi[1][0] = -4. * xsi3 + 1.;

    dphi[0][4] = 4. * xsi2;
    dphi[1][4] = 4. * xsi1;

    dphi[0][5] = -4. * xsi2;
    dphi[1][5] = 4. * (xsi3 - xsi2);

    dphi[0][3] = 4. * (xsi3 - xsi1);
    dphi[1][3] = -4. * xsi1;

    // element conectivity
    //     2
    //     54
    //     031
     
    return;
}

//------------------------------------------------------------------------------
//----------------COMPUTE SHAPE FUNCTION SECOND DERIVATIVE VALUE----------------
//------------------------------------------------------------------------------
template<>
void QuadShapeFunction<2>::evaluateHessian(double xi[], double ***ddphi) const {

     // ddphi(0,0)(0) = 4.;
     // ddphi(0,1)(0) = 4.;
     // ddphi(1,0)(0) = 4.;
     // ddphi(1,1)(0) = 4.;

     // ddphi(0,0)(1) = 4.;
     // ddphi(0,1)(1) = 0.;
     // ddphi(1,0)(1) = 0.;
     // ddphi(1,1)(1) = 0.;

     // ddphi(0,0)(2) = 0.;
     // ddphi(0,1)(2) = 0.;
     // ddphi(1,0)(2) = 0.;
     // ddphi(1,1)(2) = 4.;

     // ddphi(0,0)(3) = -8.;
     // ddphi(0,1)(3) = -4.;
     // ddphi(1,0)(3) = -4.;
     // ddphi(1,1)(3) = 0.;

     // ddphi(0,0)(4) = 0.;
     // ddphi(0,1)(4) = 4.;
     // ddphi(1,0)(4) = 4.;
     // ddphi(1,1)(4) = 0.;

     // ddphi(0,0)(5) = 0.;
     // ddphi(0,1)(5) = -4.;
     // ddphi(1,0)(5) = -4.;
     // ddphi(1,1)(5) = -8.;

     // element conectivity
     //     2
     //     54
     //     031
     
    return;
}

#endif
