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

/// Defines the quadratic shape functions and its derivatives

template<int DIM, int DEG>
class QuadShapeFunction{
public:
    //Number of element notes    
    static const int numIntpNodes = 4*DIM-2;

    /// Evaluates the shape function value
    /// @param double* Non-dimensional coordinates 
    /// @param double* Shape function values
    void evaluate(VecDouble &xi, double *phi) const;
    
    /// Evaluates the values of the shape funtion derivatives
    /// @param double* Non-dimensional coordinates 
    /// @param double** Shape function derivatives values
    void evaluateGradient(VecDouble &xi, double **dphi) const;   

    /// Evaluates the values of the shape funtion second derivatives    
    /// @param double** Shape function second derivatives values
    void evaluateHessian(double ***ddphi) const;     

    /// Gets the element nodal coordinates
    /// @param double** Non-dimensional coordinate
    void getCoordinates(double** &coord) const;
};


//------------------------------------------------------------------------------
//-------------------------COMPUTE SHAPE FUNCTION VALUE-------------------------
//------------------------------------------------------------------------------
template<>
void QuadShapeFunction<2,1>::evaluate(VecDouble &xi, double *phi) const {

    double xsi1 = xi[0];
    double xsi2 = xi[1];
    double xsi3 = 1. - xsi1 - xsi2;
    
    phi[0] = xsi3;
    phi[1] = xsi1;
    phi[2] = xsi2;
    
    // element conectivity
    //     2
    //     01

    return;
}

// Defines quadratic shape functions and its derivatives 
// for triangles and tetrahedrons
template<>
void QuadShapeFunction<2,2>::evaluate(VecDouble &xi, double *phi) const {
    
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

template<>
void QuadShapeFunction<3,1>::evaluate(VecDouble &xi, double *phi) const {

    double xsi1 = xi[0];
    double xsi2 = xi[1];
    double xsi3 = xi[2];

    phi[0] = 1.0 - xsi1 - xsi2 - xsi3;
    phi[1] = xsi1;
    phi[2] = xsi2;
    phi[3] = xsi3;

    return;
}

//------------------------------------------------------------------------------
//-------------------------COMPUTE SHAPE FUNCTION VALUE-------------------------
//------------------------------------------------------------------------------
// Defines quadratic shape functions and its derivatives 
// for triangles and tetrahedrons
template<>
void QuadShapeFunction<3,2>::evaluate(VecDouble &xi, double *phi) const {
    
    double xsi1 = xi[0];
    double xsi2 = xi[1];
    double xsi3 = xi[2];

    phi[3] = 2.0 * (xsi1 - 0.50) * xsi1;
    phi[1] = 2.0 * (xsi2 - 0.50) * xsi2;
    phi[2] = 2.0 * (xsi3 - 0.50) * xsi3;
    phi[0] = (2.0 - 2.0 * xsi3 - 2.0 * xsi2 - 2.0 * xsi1 - 1.0) * (1.0 - xsi1 - xsi2 - xsi3);
    phi[8] = 4.0 * xsi1 * xsi2;
    phi[9] = 4.0 * xsi1 * xsi3;
    phi[7] = 4.0 * xsi1 * (1.0 - xsi1 - xsi2 - xsi3);
    phi[6] = 4.0 * xsi3 * (1.0 - xsi1 - xsi2 - xsi3);
    phi[4] = 4.0 * xsi2 * (1.0 - xsi1 - xsi2 - xsi3);
    phi[5] = 4.0 * xsi2 * xsi3;

    // element conectivity
    //layer 1
    //     3
    //     8 9
    //     2 7 1
    //layer 2
    //     6
    //     54
    //layer 3
    //     0
    
     
     return;
}

//------------------------------------------------------------------------------
//-------------------COMPUTE SHAPE FUNCTION DERIVATIVE VALUE--------------------
//------------------------------------------------------------------------------
template<>
void QuadShapeFunction<2,1>::evaluateGradient(VecDouble &xi, double **dphi) const {

    dphi[0][0] = -1.;
    dphi[0][1] = -1.;
    
    dphi[1][0] = 1.;
    dphi[1][1] = 0.;
    
    dphi[2][0] = 0.;
    dphi[2][1] = 1.;
    
    // element conectivity
    //     2
    //     01

    return;
}

template<>
void QuadShapeFunction<2,2>::evaluateGradient(VecDouble &xi, double **dphi) const {

    const double xsi1 = xi[0];
    const double xsi2 = xi[1];
    const double xsi3 = 1. - xsi1 - xsi2;

    dphi[1][0] = 4. * xsi1 - 1.;
    dphi[1][1] = 0.;

    dphi[2][0] = 0.;
    dphi[2][1] = 4. * xsi2 - 1.;

    dphi[0][0] = -4. * xsi3 + 1.;
    dphi[0][1] = -4. * xsi3 + 1.;

    dphi[4][0] = 4. * xsi2;
    dphi[4][1] = 4. * xsi1;

    dphi[5][0] = -4. * xsi2;
    dphi[5][1] = 4. * (xsi3 - xsi2);

    dphi[3][0] = 4. * (xsi3 - xsi1);
    dphi[3][1] = -4. * xsi1;

    // element conectivity
    //     2
    //     54
    //     031
     
    return;
}

template<>
void QuadShapeFunction<3,1>::evaluateGradient(VecDouble &xi, double **dphi) const {

    dphi[0][0] = -1.0;
    dphi[0][1] = -1.0;
    dphi[0][2] = -1.0;
    
    dphi[1][0] = 1.0;
    dphi[1][1] = 0.0;
    dphi[1][2] = 0.0;

    dphi[2][0] = 0.0;
    dphi[2][1] = 1.0;
    dphi[2][2] = 0.0;

    dphi[3][0] = 0.0;
    dphi[3][1] = 0.0;
    dphi[3][2] = 1.0;


    return;
}

template<>
void QuadShapeFunction<3,2>::evaluateGradient(VecDouble &xi, double **dphi) const {

    const double xsi1 = xi[0];
    const double xsi2 = xi[1];
    const double xsi3 = xi[2];

    dphi[3][0] = 4. * xsi1 - 1.;
    dphi[3][1] = 0.;
    dphi[3][2] = 0.;

    dphi[1][0] = 0.;
    dphi[1][1] = 4. * xsi2 - 1.;
    dphi[1][2] = 0.;
     
    dphi[2][0] = 0.;
    dphi[2][1] = 0.;
    dphi[2][2] = 4. * xsi3 - 1.;

    dphi[0][0] = 4. * (xsi1 + xsi2 + xsi3) - 3.;
    dphi[0][1] = 4. * (xsi1 + xsi2 + xsi3) - 3.;
    dphi[0][2] = 4. * (xsi1 + xsi2 + xsi3) - 3.;

    dphi[8][0] = 4. * xsi2;
    dphi[8][1] = 4. * xsi1;
    dphi[8][2] = 0.;

    dphi[9][0] = 4. * xsi3;
    dphi[9][1] = 0.;
    dphi[9][2] = 4. * xsi1;

    dphi[7][0] = 4. * (1. - 2. * xsi1 - xsi2 - xsi3);
    dphi[7][1] = -4. * xsi1;
    dphi[7][2] = -4. * xsi1;

    dphi[5][0] = 0.;
    dphi[5][1] = 4. * xsi3;
    dphi[5][2] = 4. * xsi2;

    dphi[6][0] = -4. * xsi3;
    dphi[6][1] = -4. * xsi3;
    dphi[6][2] = 4. * (1. - 2. * xsi3 - xsi2 - xsi1);

    dphi[4][0] = -4. * xsi2;
    dphi[4][1] = 4. * (1. - 2. * xsi2 - xsi1 - xsi3);
    dphi[4][2] = -4. * xsi2;





    return;
}

//------------------------------------------------------------------------------
//----------------COMPUTE SHAPE FUNCTION SECOND DERIVATIVE VALUE----------------
//------------------------------------------------------------------------------
template<>
void QuadShapeFunction<2,1>::evaluateHessian(double ***ddphi) const {

    return;
}

template<>
void QuadShapeFunction<2,2>::evaluateHessian(double ***ddphi) const {

    ddphi[0][0][0] = 4.;
    ddphi[0][1][0] = 4.;
    ddphi[1][0][0] = 4.;
    ddphi[1][1][0] = 4.;

    ddphi[0][0][1] = 4.;
    ddphi[0][1][1] = 0.;
    ddphi[1][0][1] = 0.;
    ddphi[1][1][1] = 0.;

    ddphi[0][0][2] = 0.;
    ddphi[0][1][2] = 0.;
    ddphi[1][0][2] = 0.;
    ddphi[1][1][2] = 4.;

    ddphi[0][0][3] = -8.;
    ddphi[0][1][3] = -4.;
    ddphi[1][0][3] = -4.;
    ddphi[1][1][3] = 0.;

    ddphi[0][0][4] = 0.;
    ddphi[0][1][4] = 4.;
    ddphi[1][0][4] = 4.;
    ddphi[1][1][4] = 0.;

    ddphi[0][0][5] = 0.;
    ddphi[0][1][5] = -4.;
    ddphi[1][0][5] = -4.;
    ddphi[1][1][5] = -8.;

     // element conectivity
     //     2
     //     54
     //     031
     
    return;
}

template<>
void QuadShapeFunction<3,1>::evaluateHessian(double ***ddphi) const {

    return;
}



template<>
void QuadShapeFunction<3,2>::evaluateHessian(double ***ddphi) const {

    ddphi[0][0][0] = 4.;
    ddphi[0][1][0] = 0.;
    ddphi[0][2][0] = 0.;
    ddphi[1][0][0] = 0.;
    ddphi[1][1][0] = 0.;
    ddphi[1][2][0] = 0.;
    ddphi[2][0][0] = 0.;
    ddphi[2][1][0] = 0.;
    ddphi[2][2][0] = 0.;

    ddphi[0][0][1] = 0.;
    ddphi[0][1][1] = 0.;
    ddphi[0][2][1] = 0.;
    ddphi[1][0][1] = 0.;
    ddphi[1][1][1] = 4.;
    ddphi[1][2][1] = 0.;
    ddphi[2][0][1] = 0.;
    ddphi[2][1][1] = 0.;
    ddphi[2][2][1] = 0.;

    ddphi[0][0][2] = 0.;
    ddphi[0][1][2] = 0.;
    ddphi[0][2][2] = 0.;
    ddphi[1][0][2] = 0.;
    ddphi[1][1][2] = 0.;
    ddphi[1][2][2] = 0.;
    ddphi[2][0][2] = 0.;
    ddphi[2][1][2] = 0.;
    ddphi[2][2][2] = 4.;

    ddphi[0][0][3] = 4.;
    ddphi[0][1][3] = 4.;
    ddphi[0][2][3] = 4.;
    ddphi[1][0][3] = 4.;
    ddphi[1][1][3] = 4.;
    ddphi[1][2][3] = 4.;
    ddphi[2][0][3] = 4.;
    ddphi[2][1][3] = 4.;
    ddphi[2][2][3] = 4.;

    ddphi[0][0][4] = 0.;
    ddphi[0][1][4] = 4.;
    ddphi[0][2][4] = 0.;
    ddphi[1][0][4] = 4.;
    ddphi[1][1][4] = 0.;
    ddphi[1][2][4] = 0.;
    ddphi[2][0][4] = 0.;
    ddphi[2][1][4] = 0.;
    ddphi[2][2][4] = 0.;

    ddphi[0][0][5] = 0.;
    ddphi[0][1][5] = 0.;
    ddphi[0][2][5] = 4.;
    ddphi[1][0][5] = 0.;
    ddphi[1][1][5] = 0.;
    ddphi[1][2][5] = 0.;
    ddphi[2][0][5] = 4.;
    ddphi[2][1][5] = 0.;
    ddphi[2][2][5] = 0.;

    ddphi[0][0][6] = -8.;
    ddphi[0][1][6] = -4.;
    ddphi[0][2][6] = -4.;
    ddphi[1][0][6] = -4.;
    ddphi[1][1][6] = 0.;
    ddphi[1][2][6] = 0.;
    ddphi[2][0][6] = -4.;
    ddphi[2][1][6] = 0.;
    ddphi[2][2][6] = 0.;

    ddphi[0][0][7] = 0.;
    ddphi[0][1][7] = 0.;
    ddphi[0][2][7] = 0.;
    ddphi[1][0][7] = 0.;
    ddphi[1][1][7] = 0.;
    ddphi[1][2][7] = 4.;
    ddphi[2][0][7] = 0.;
    ddphi[2][1][7] = 4.;
    ddphi[2][2][7] = 0.;

    ddphi[0][0][8] = 0.;
    ddphi[0][1][8] = 0.;
    ddphi[0][2][8] = -4.;
    ddphi[1][0][8] = 0.;
    ddphi[1][1][8] = 0.;
    ddphi[1][2][8] = -4.;
    ddphi[2][0][8] = -4.;
    ddphi[2][1][8] = -4.;
    ddphi[2][2][8] = -8.;

    ddphi[0][0][9] = 0.;
    ddphi[0][1][9] = -4.;
    ddphi[0][2][9] = 0.;
    ddphi[1][0][9] = -4.;
    ddphi[1][1][9] = -8.;
    ddphi[1][2][9] = -4.;
    ddphi[2][0][9] = 0.;
    ddphi[2][1][9] = -4.;
    ddphi[2][2][9] = 0.;


    return;
}


template<>
void QuadShapeFunction<2,1>::getCoordinates(double** &coord) const {

    coord[0][0] = 0.0;
    coord[1][0] = 0.0;

    coord[0][1] = 1.0;
    coord[1][1] = 0.0;

    coord[0][2] = 0.0;
    coord[1][2] = 1.0;

    return;
}

template<>
void QuadShapeFunction<2,2>::getCoordinates(double** &coord) const {

    coord[0][0] = 0.0;
    coord[1][0] = 0.0;

    coord[0][1] = 1.0;
    coord[1][1] = 0.0;

    coord[0][2] = 0.0;
    coord[1][2] = 1.0;

    coord[0][3] = 0.5;
    coord[1][3] = 0.0;

    coord[0][4] = 0.5;
    coord[1][4] = 0.5;

    coord[0][5] = 0.0;
    coord[1][5] = 0.5;

    return;
}

template<>
void QuadShapeFunction<3,1>::getCoordinates(double** &coord) const {

    coord[0][0] = 0.0;
    coord[1][0] = 0.0;
    coord[2][0] = 0.0;

    coord[0][1] = 1.0;
    coord[1][1] = 0.0;
    coord[2][1] = 0.0;

    coord[0][2] = 0.0;
    coord[1][2] = 1.0;
    coord[2][2] = 0.0;

    coord[0][3] = 0.0;
    coord[1][3] = 0.0;
    coord[2][3] = 1.0;

    return;
}

template<>
void QuadShapeFunction<3,2>::getCoordinates(double** &coord) const {

    coord[0][0] = 0.0;
    coord[1][0] = 0.0;
    coord[2][0] = 0.0;

    coord[0][1] = 0.0;
    coord[1][1] = 1.0;
    coord[2][1] = 0.0;

    coord[0][2] = 0.0;
    coord[1][2] = 0.0;
    coord[2][2] = 1.0;

    coord[0][3] = 1.0;
    coord[1][3] = 0.0;
    coord[2][3] = 0.0;

    coord[0][4] = 0.0;
    coord[1][4] = 0.5;
    coord[2][4] = 0.0;

    coord[0][5] = 0.0;
    coord[1][5] = 0.5;
    coord[2][5] = 0.5;

    coord[0][6] = 0.0;
    coord[1][6] = 0.0;
    coord[2][6] = 0.5;

    coord[0][7] = 0.5;
    coord[1][7] = 0.0;
    coord[2][7] = 0.0;

    coord[0][8] = 0.5;
    coord[1][8] = 0.5;
    coord[2][8] = 0.0;

    coord[0][9] = 0.5;
    coord[1][9] = 0.0;
    coord[2][9] = 0.5;
    return;
}
#endif
